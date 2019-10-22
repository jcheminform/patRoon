#' @include main.R

makeGraph <- function(components, onlyLinked, titles)
{
    cInfo <- copy(components@componentInfo)
    
    # convert link IDs to numeric indices
    cInfo[, id := .I]
    cInfo[, linksIDs := lapply(links, match, table = names(components))]
    allLinks <- unique(unlist(cInfo$linksIDs))
    cInfo <- cInfo[lengths(links) > 0 | id %in% allLinks]
    
    edges <- rbindlist(mapply(cInfo$name, cInfo$linksIDs, FUN = function(n, l)
    {
        data.table::data.table(from = n, to = cInfo$name[match(unlist(l), cInfo$id)])
    }, SIMPLIFY = FALSE))
    
    graph <- igraph::simplify(igraph::graph_from_data_frame(edges, directed = FALSE))
    fc <- igraph::fastgreedy.community(graph)
    
    data <- visNetwork::toVisNetworkData(graph)
    nodes <- as.data.table(data$nodes)
    nodes[, group := fc$membership]
    
    if (!onlyLinked)
    {
        unNodes <- data.table(id = setdiff(names(components), cInfo$name), group = 0)
        unNodes[, label := id]
        nodes <- rbind(nodes, unNodes)
    }
    
    nodes[, shape := "circle"]
    nodes[, title := titles[match(id, names(components))]]
    
    visNetwork::visNetwork(nodes = nodes, edges = data$edges)
}
