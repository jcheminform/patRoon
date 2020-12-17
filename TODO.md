# Release

## general
- test negative subset indices
- refs to OpenBabel?
- convertMSFiles()
    - Agilent .d is also a directory?
    - Remove necessity to have different input/output formats? (at least OK for pwiz)
    - Support OpenMS vendor conversion? (eg thermo)

## AutoID

- ID level rules
    - add scorings for SIRIUS/DA
- interface
    - also convert TASQ?
    - newProject()
        - also allow suspect annotation with only peak lists? currently only selectable if formulas/compounds selected
    - annotateSuspects()
        - check why it's is sometimes slow
            - seems to be logging, disable by default? --> only slow with testthat?
    - filter():
        - cache?
    - don't assign level <1 if suspect is a target? or give the choice (or make filter?)
    - spec similarity:
        - port from TPs someday
- misc
    - prepareSuspectList(): export?
        - mainly to allow merging of lists, perhaps make util for that instead? Would also be handy for MF databases
            - could also fix column names, replace "-" with NAs etc
        - if yes, mention in ref docs for screenSuspects()
- expand reporting
    - eg include suspect name in EICs
    - mention suspect similarities/ranks etc for candidates (or somehow in compounds?)
    - optionally report with collapsed suspects
- update docs and handbook
    - mention that components should be done prior to onlyHits=T?


## docs
- improve instructions for MF and SIRIUS installation?


## sets
- methods to implement
    - consensus / comparison()
        - could first make consensus of setObjects and then make new set object from that
        - for compounds (and formulas?) need proper way to average scorings
        - or for later?
    - more sub class specific methods
        - compoundsSetMF sub-class (for settings slot)?
    - provide methods for non-implemented functionality
        - consensus()? (see above)
        - compoundViewer?
        - groupFeaturesXCMS3?
- misc
    - as.data.table() for formulas: average=T will now produce strange averaged ionized formula, for now simply remove this column...
        - also give a note in docs?
        - or maybe only remove is not all adducts are equal?
    - handle errors when object has <=1 set
        - groupFeaturesScreening()
        - mergeScreeningSetInfos()
    - fix empty MS(MS) peaklists if unavailable during merging sets
        - already fixed?
- merging setObjects
    - check if more has to be cached and may need status messages
    - filter features, annotation results on minimum abundance amongst different setObjects
        - (unlike setThreshold also take objects without results in to account)
        - add new filter()?
    - compound set consensus: scoreRanges should be re-determined from annotation results?
- setThreshold
    - filter() argument
    - keep/change default setThreshold?
        - or remove argument from generators?
- suspect screening
    - implement TASQ?
    - consensus?
        - or just new set filter described above?
    - annotation columns not in report, fine? (there are many) If yes document
    - as.data.table(fGroupsScrSets, collapseSuspects=NULL): omits sets column
        - still true? check
- neutralizing / ionization
    - mergeIons()
        - other name?
        - makes sense to not choose monoisotopic mass? nor for annotation at least
		- also method for fGroupsSets?
			- do it per component set
			- used to update adducts
			- needs to re-group afterwards
			    - make general reGroup() method, that uses stored settings from groupFeatures()?
	    - prefer adducts based on MS/MS? eg handy for Na/K adducts
	- remove adduct slots
	- unset()
	    - can now only work for 1 set?
	    - update all methods
	- formula/compounds: get adduct from gInfo if present
		- also for screening? could (optionally) look for matches with neutralized fGroup/susp masses
		- adduct arg still overrides?
	- cliqueMS components
	- component selection tool/function
	    - otherwise perhaps make a fGroup remover function to help subsetting
			- similarly as for feature remover in fGroups...
			- del()/delete()/rmResult()/delResult() generic? could also be for other classes
		- similarly: set() like method to change data, such as adduct annotations
- NEWS
    - [..., reAverage = FALSE] and implications of filtering when setting it to TRUE
- docs
    - filter() for features/fGroups: apply to neutral masses
    - CAMERA/RAMClustR/nontarget components: clearly mention it is simply a merge between sets
    - intclust is not a componentsSet
    - find nice way to re-use docs
    - mention that setObjects are _not_ filtered by setThreshold for formulas/compounds
    - mention that new consensus for formulas/compounds is made after filter() and addFormulaScoring()
        - this could mean really different results if subsetting on sets is done prior to filtering
        - put message() in sync function?
    - mention that set coverage/formula feature coverages do not consider sets/analyses without any results
        - put message() in sync function?
    - document for every object how consensus/merge is done
    - improve docs for areas (only affects when features=FALSE) and average (different behavior when features=TRUE/FALSE) for as.data.table() of featureGroups
    - update/check version nr mentioned in filter() for MSPeakLists
    - explain xlim/ylim behavior for annotations/mols for plotSpec()
    - update/add aliases
- tests
    - handle/test empty objects
    - test DA algorithms
    


## features
- feature optim:
    - docs
        - mention parameters default unless specified
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- suspect screening
    - rename patRoonData::targets?
    - rename groupFeaturesScreening?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- importFeaturesXCMS/importFeaturesXCMS3/importFeatureGroupsXCMS: get rid of anaInfo arg requirement? (or make import func?)
- comparison(): support xcms3? (needs missing support for missing raw data)
- Fix: blank filter with multiple replicate groups (and maybe others?)
- Check: units of plotChord() rt/mz graphs seems off
- plotEIC(): get mzWindow from features if possible to show more representative results (eg when OpenMS mz window is very small)
- remove features not in any group from fGroups


## MSPeakLists
- isotope tagging is lost after averaging
- collapse averagedPeakLists
- test avg params
- metadata() generic?


## compounds
- SIRIUS: use --auto-charge instead of manually fixing charge of fragments (or not? conflicting docs on what it does)
- test score normalization?
- timeouts for SIRIUS?
- do something about negative H explained fragments by MF?
- SusDat MF support
- MetFrag: auto-include suspect results if suspectListScore is selected?


## formulas
- customize/document ranking column order? (only do rank for sirius?)
- getFormInfoList(): take care of consensus results like getPrecursorFormScores()

## components
- RC: check spearmans correlation
- NT: minimum size argument, combine rows for multiple rGroups?


## reporting
- add more options to reportPlots argument of reportHTML()?


## Cleanup
- Reduce non-exported class only methods

## MP

- future MP
    - delayBetweenProc?


# Future

## General

- msPurity integration
- suspect screening: add MS/MS qualifiers
- fillPeaks for CAMERA (and RAMClustR?)
- support fastcluster for compounds clustering/int component clusters?
- algorithmObject() generic: for xset, xsa, rc, ...
- newProject(): fix multi line delete (when possible)
- more withr wrapping? (dev, par)
- improve default plotting for plotInt and cluster plot functions
- newProject()
    - concentration column for anaInfo
    - generate more detailed script with e.g. commented examples of subsetting, extraction etc
- support more of the new SIRIUS functionality
	- newProject(): import Bruker seq file?


## Features

- integrate OpenMS feature scoring and isotopes and PPS in general (also include filters?)
- parallel enviPick
- OpenMS MetaboliteAdductDecharger support?
- OpenMS: Support KD grouper?
- suspect screening: tag fGroups with suspect instead of converting fGroups object (and add filter to remove non-hits)
- Integration of mzMine features (package pending...), MS-DIAL and KPIC2, peakonly, SIRIUS?


## MSPeakLists

- DA
    - generateMSPeakListsDA: find precursor masses with larger window
    - tests
        - utils? EICs with export/vdiffr?
        - test MS peak lists deisotoping?
- metadata for Bruker peaklists?


## Formulas

- DBE calculation for SIRIUS?
- OM reporting
- as.data.table: option to average per replicate group?


## Compounds

- do something with sirius fingerprints? --> comparison?
- fix compoundViewer
- add new MF HD scorings and make sure default normalization equals that of MF web
- CFM-ID and MS-FINDER integration
- utility functions to make custom DBs for MetFrag and SIRIUS and support to use them with the latter


## components
- mass defect components
- CliqueMS
- split peak correlation and adduct etc annotation? would allow better non-target integration
- intclust
    - optionally take areas instead of intensities
    - cache results

