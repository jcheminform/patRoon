#' @include main.R
NULL

getBTBin <- function()
{
    ret <- path.expand(getOption("patRoon.path.BioTransformer", ""))
    if (is.null(ret) || !nzchar(ret) || !file.exists(ret))
        stop("Please set the 'biotransformer' option with a (correct) path to the BioTransformer JAR file. Example: options(patRoon.path.BioTransformer = \"C:/biotransformerjar/biotransformer-1-0-8.jar\")")
    
    if (!nzchar(Sys.which("java")))
        stop("Please make sure that java is installed and its location is correctly set in PATH.")
    
    return(ret)
}

initBTCommand <- function(SMILES, type, steps, extraOpts, btBin, logFile)
{
    outFile <- tempfile("btresults", fileext = ".csv")
    
    args <- c("-b", type,
              "-k", "pred",
              "-ismi", SMILES,
              "-s", as.character(steps),
              "-ocsv", outFile,
              extraOpts)
    
    return(list(command = "java", args = c("-jar", btBin, args), logFile = logFile, outFile = outFile, SMILES = SMILES))
}

processBTResults <- function(cmd)
{
    ret <- fread(cmd$outFile)
    
    # UNDONE: transform column names, more?
    
    return(ret)
}

predictTPsBioTransformer <- function(suspects, type = "env", steps = 2, extraOpts = NULL,
                                     logPath = file.path("log", "biotransformer"),
                                     maxProcAmount = getOption("patRoon.maxProcAmount"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDataFrame(suspects, any.missing = FALSE, min.rows = 1, add = ac)
    assertHasNames(suspects, c("name", "SMILES"), add = ac)
    checkmate::assertChoice(type, c("ecbased", "cyp450", "phaseII", "hgut", "superbio", "allHuman", "env"), add = ac)
    checkmate::assertCount(steps, positive = TRUE, add = ac)
    checkmate::assertCharacter(extraOpts, null.ok = TRUE, add = ac)
    assertMultiProcArgs(logPath, maxProcAmount, add = ac)
    checkmate::reportAssertions(ac)
    
    cacheDB <- openCacheDBScope()
    baseHash <- makeHash(type, steps, extraOpts)
    setHash <- makeHash(suspects, baseHash)
    hashes <- sapply(suspects$SMILES, makeHash, baseHash)
    cachedSet <- loadCacheSet("predictTPsBT", setHash, cacheDB)
    
    btBin <- getBTBin()
    btPath <- dirname(btBin) # BT has to be executed from its own directory
    
    if (!is.null(logPath))
    {
        mkdirp(logPath)
        logPath <- normalizePath(logPath) # need full path since we will temporarily change the working directory.
    }
    
    cachedResults <- sapply(hashes, function(hash)
    {
        ret <- NULL
        if (!is.null(cachedSet))
            ret <- cachedSet[[hash]]
        if (is.null(ret))
            ret <- loadCacheData("predictTPsBT", hash, cacheDB)
        return(ret)
    }, simplify = FALSE)
    names(cachedResults) <- suspects$name
    cachedResults <- cachedResults[!sapply(cachedResults, is.null)]

    doSuspects <- suspects[!SMILES %in% names(cachedResults)]
    cmdQueue <- mapply(doSuspects$name, doSuspects$SMILES, SIMPLIFY = FALSE, FUN = function(n, sm)
    {
        logf <- if (!is.null(logPath)) file.path(logPath, paste0("biotr-", n, ".txt")) else NULL
        initBTCommand(sm, type = type, steps = steps, extraOpts = extraOpts, btBin = btBin, logFile = logf)
    })
    names(cmdQueue) <- doSuspects$name
    
    results <- list()
    
    if (length(cmdQueue) > 0)
    {
        withr::with_dir(btPath, {
            results <- executeMultiProcess(cmdQueue, function(cmd)
            {
                res <- processBTResults(cmd)
                saveCacheData("predictTPsBT", res, hashes[[cmd$SMILES]], cacheDB)
                return(res)
            }, maxProcAmount = maxProcAmount)
        })
        
        if (is.null(cachedSet))
            saveCacheSet("predictTPsBT", hashes, setHash, cacheDB)
    }
    
    if (length(cachedResults) > 0)
    {
        results <- c(results, cachedResults)
        results <- results[intersect(suspects$names, names(results))] # re-order
    }
    
    # UNDONE: make nice S4 class?
    
    return(results)
}
