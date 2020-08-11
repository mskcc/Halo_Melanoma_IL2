#' Load global variables containing previously processed/calculated data 
#'
#' Given just the study configuration in list form, load into the global namespace
#' study data including parsed cell types, marker names, condition analysis list, 
#' sample annotation, annotated cell data, parsed study questions, neighborhood counts
#' and neighborhood distances. NOTE: this assumes all meta data have been processed,
#' cells have been annotated and metrics have been calculated prior to running
#' this function.
#'
#' @param config                  study configuration in list form (see docs for details)
#' @param cellTypes               logical; when TRUE, parse and expand cell types XLSX 
#'                                file and store in variable 'cellTypes'; default: TRUE
#' @param markerList              logical; when TRUE, pull Marker_name column from *Markers.xlsx 
#'                                file and store in variable 'markers'; default: TRUE  
#' @param analyses                logical; when TRUE, read stats conditions XLSX file into variable
#'                                'analysisList'; default: TRUE
#' @param conditions              logical; when TRUE, load conditions index into variable 'conds'; default: TRUE
#' @param sampleAnnotation        logical; when TRUE, load all sample annotations into variable 'sampAnn';
#'                                default: TRUE
#' @param annotatedCells          logical; when TRUE, load annotated cell data into variable 'annCells';
#'                                default: TRUE
#' @param questions               logical; when TRUE, parse stats questions XLSX file into list form and 
#'                                store in variable 'allQuestions'; default=FALSE
#' @param neighborhoodCounts      logical; when TRUE, read neighborhood counts file into variable 'nbhdCounts';
#'                                default: FALSE
#' @param neighborhoodDistances   logical; when TRUE, read neighborhood distances file into variable 
#'                                'nbhdDistances'; default: FALSE
#'
#' @return nothing
loadGlobalStudyData <- function(config,
                                all = FALSE,
                                cellTypes = TRUE,
                                markerList = TRUE,
                                sampleAnnotation = TRUE,
                                annotatedCells = TRUE,
                                analyses = FALSE,
                                conditions = FALSE,
                                questions = FALSE,
                                neighborhoodCounts = FALSE,
                                cellsInTumorNeighborhood = FALSE,
                                tmeSampleStatus = FALSE,
                                tmeCellStatus = FALSE){
    if(all){
        cellTypes <- markerList <- analyses <- conditions <- sampleAnnotation <- TRUE
        annotatedCells <- questions <- neighborhoodCounts <- cellsInTumorNeighborhood <- TRUE
        #tmeSampleStatus <- TRUE 
        tmeCellStatus <- TRUE
    }

    if(is.null(config$cell_dive_id)){
        config$cell_dive_id <- "All"
    }

    if(!is.null(config$debug) && (config$debug == "yes" || config$debug)){
        log_threshold(DEBUG)
    }

    if(cellTypes){
        log_info("Loading global variable: cellTypes")
        cellTypes <<- getCellTypes(getFiles(path = config$meta_dir, 
                                            files = config$meta_files, 
                                            pattern="_CellTypes.xlsx"))
    }

    if(markerList){
        log_info("Loading global variable: markers")
        markers <<- read.xlsx(getFiles(path = config$meta_dir,
                                       files = config$meta_files,
                                       pattern="_Markers.xlsx"), 
                              1, check.names = F) %>%
                    pull(Marker_name)
    }

    if(questions){
        log_info("Loading global variable: allQuestions")
        allQuestions <<- parseStatsQuestions(config$statistics_questions_file)
        if(!is.null(config$question)){
            allQuestions <<- allQuestions[names(allQuestions) == config$question]
        }
    }

    if(conditions || analyses){
        log_info("Loading global variable: conds")
        conds <<- getConditionsIndex(config$statistics_conditions_index,
                                     config$statistics_conditions_file,
                                     config$arrange_annotation)
    }

    if(analyses){
        log_info("Loading global variable: analysisList")
        analysisList <<- getStatsAnalyses(config$statistics_conditions_file, conds)
    }

    if(sampleAnnotation || annotatedCells){
        log_info("Loading global variable: sampAnn")
        sampAnn <<- loadStudyAnnotations(getFiles(path = config$meta_dir, pattern = ".xlsx"))$flat %>%
                                         mutate(Sample = FOV_ID) %>%
                                         filter(is.na(FOV_exclusion))
        if(!is.null(config$cell_dive_id) && config$cell_dive_id != "All"){
            sampAnn <<- sampAnn %>% filter(CellDive_ID == config$cell_dive_id)
        }
    }

    if(annotatedCells || neighborhoodCounts || tmeSampleStatus || tmeCellStatus){
        log_info("Loading global variable: annCells")
        annCells <<- loadAnnotatedCells(sampAnn, config$cell_data_dir, 
                                        bandDir = config$band_dir, 
                                        fovAreaDir = config$fov_area_dir, 
                                        cellDiveID = config$cell_dive_id)
    }

    if(neighborhoodCounts){
        log_info("Loading global variabel: nbhdCounts")
        neighborhoodDirs <- file.path(config$neighborhood_dir, c("C2","C3"))

        nbhdCounts <<- loadMacroNeighborhoodData(neighborhoodDirs, 
                                                 analysisList,
                                                 config$neighborhood_data_dir,
                                                 sampAnn, 
                                                 cellDiveID = config$cell_dive_id)

        nbhdCounts <<- formatNeighborhoodCounts(config$neighborhood_data_dir, annCells, cellDiveID = config$cell_dive_id)
    }

    if(cellsInTumorNeighborhood | tmeSampleStatus | tmeCellStatus){

        log_info("Loading global variable: tumorNbhd")
        nbhdDir <- file.path(config$neighborhood_dir, "C2")
        tFile   <- file.path(nbhdDir, dir(nbhdDir)[grepl("Tumor", dir(nbhdDir))]) 
        tumorNbhd <<- readRDS(tFile)

        if(!is.null(config$cell_dive_id)){
            fovs <- sampAnn %>% pull(FOV_ID) %>% unique()
            tumorNbhd <<- tumorNbhd %>% filter(FOV %in% fovs)
        }

        tumorNbhdCells <<- c()
        if(nrow(tumorNbhd) > 0){
            log_info("Loading global variable: tumorNbhdCells")
            tumorNbhdCells <<- tumorNbhd %>% pull(N.UUID) %>% unique()
        }
    }

    if(tmeSampleStatus){

        log_info("Loading TME assignments by sample & adding to global variable: tumorNbhd")
        files <- file.path(config$tme_by_sample_dir, dir(config$tme_by_sample_dir))

        for(fl in files){
            tme <- readRDS(fl)
            mrkrCol <- names(tme)[grepl("microEnv", names(tme))]
            mrkr <- gsub("microEnv\\.","",mrkrCol)
            annCells <<- annCells %>%
                         left_join(tme %>%
                                   select_at(c("UUID", mrkrCol)) %>%
                                   rename_at(vars(mrkrCol), list(~(. = paste0("sample_",mrkr)))),
                                   by = "UUID")
        }
    }

    if(tmeCellStatus){

        log_info("Loading TME assignments by cell & adding to global variable: tumorNbhd")
        files <- file.path(config$tme_by_cell_dir, dir(config$tme_by_cell_dir))

        for(fl in files){
            tme <- readRDS(fl)
            mrkrCol <- names(tme)[grepl("microEnv", names(tme))]
            mrkr <- gsub("microEnv\\.","",mrkrCol)
            annCells <<- annCells %>%
                         left_join(tme %>%
                                   select_at(c("UUID", mrkrCol)) %>%
                                   rename_at(vars(mrkrCol), list(~(. = paste0("cell_",mrkr)))),
                                   by = "UUID")
        }
    }

}

#' Generate a single tibble from several files
#' 
#' Provided a directory path or a vector of files, read all RDA files
#' and create a single tibble by binding the rows of all individual tibbles
#'
#' @param dataDir   full path to directory of RDA files; all files found in directory
#'                  will be included (default = NULL)
#' @param files     vector of paths to RDA files to be included
#' @param keepCols  vector of column names to keep; all others will be excluded from
#'                  final tibble
#' @return a single tibble 
tibbleFromMultipleFiles <- function(dataDir = NULL, files = NULL, keepCols=NULL){

    if(is.null(files) && !is.null(dataDir)){
        files <- file.path(dataDir, dir(dataDir))
    }
    tbl <- tibble()
    for(fl in files){
        flTbl <- readRDS(fl)
        if(!is.null(keepCols) && nrow(flTbl) > 0){
            flTbl <- flTbl %>% select_at(keepCols)
        }
        tbl <- tbl %>% bind_rows(flTbl %>% unique())
    }
    tbl %>% unique()
}

#' Load cell level data with all available cell and study annotations
#' 
#' Read pre-processed cell level data from one or more RDA files and join
#' additional information including sample annotation, band assignments
#' and areas
#' 
#' @param annDir       directory path containing one or more cell-level RDA files
#' @param bandDir      directory path containing one or more RDA files including 
#'                     interface band assignments, area and distance between cell
#'                     and its nearest point on the tumor interface
#' @param fovAreaDir   directory path containing one or more RDA file of total area 
#'                     of each FOV
#' @param sampAnn      tibble of all sample annotation data to be joined to cell data
#' @param cellDiveID   cellDiveID for which data should be loaded; default: All 
#' 
#' @return single tibble of all cell and study annotations
loadAnnotatedCells <- function(sampAnn, annDir, fovAreaDir = NULL, bandDir = NULL, 
                               cellDiveID = "all"){ 

    baFiles <- faFiles <- NULL

    annFiles <- file.path(annDir, dir(annDir))

    if(!is.null(bandDir)){
        baFiles  <- file.path(bandDir, dir(bandDir))
    }
    if(!is.null(fovAreaDir)){
        faFiles  <- file.path(fovAreaDir, dir(fovAreaDir))
    }

    if(!is.null(cellDiveID) && tolower(cellDiveID[1]) != "all"){
        cdidPat  <- paste0("^", cellDiveID, "_")
        annFiles <- annFiles[grepl(cdidPat, basename(annFiles))]
    }

    dat <- tibble()
    for(af in annFiles){
        log_debug(paste0("loading cells from file ",af))

        sampCells <- readRDS(af) %>%
                     left_join(sampAnn, by = intersect(names(.), names(sampAnn))) %>%
                     select(-MSK_ID_accession_num, -MSK_ID_block_num, -Replicate, -Microscope,
                            -Staining_batch, -Tumor_subtype, -Site, -dplyr::matches("detailed"),
                            -dplyr::matches("exclusion"))

        ## get band data for sample, if it exists
        cdids <- sampCells %>% pull(CellDive_ID) %>% unique()

        for(id in cdids){
            faFile <- NULL
            if(!is.null(faFiles)){ faFiles <- faFiles[grep(paste0("^", id, "_"), basename(faFiles))] }
            if(!is.null(faFile) && length(faFile) == 1 && file.exists(faFile)){
                log_debug(paste0("joining FOV areas from file: ", faFile))
                fa <- readRDS(faFile)
                sampCells <- sampCells %>%
                             left_join(fa, by = intersect(names(sampCells), names(fa)))
            }

            baFile <- NULL
            if(!is.null(baFiles)){ baFile <- baFiles[grep(paste0("^", id, "_"), basename(baFiles))] }
            if(!is.null(baFile) && length(baFile) == 1 && file.exists(baFile)){
                log_debug(paste0("joining cell band info from file: ", baFile))
                ba <- readRDS(baFile)
                sampCells <- sampCells %>%
                             left_join(ba, by = intersect(names(sampCells), names(ba)))
            }
        }

        dat <- dat %>% bind_rows(sampCells)
    }

    dat
}

#' Load or generate index of all fraction and density conditions to be analyzed
#'
#' If an XLSX index file already exists, read and return the index; if not, generate
#' an index based on the Fraction and Density sheets in a manually-created conditions
#' XLSX file and condition configuration
#' 
#' @param indexFile         XLSX file containing indexed conditions or file to which newly
#'                          generated index will be writeen
#' @param conditionsFile    XLSX file containing a sheet for Fraction conditions and one for
#'                          Density conditions (see docs for details)
#' @param conditionConfig   config in nested list form (key 'arrange_annotation' in study 
#'                          overview config file); see docs for details
#'
#' @return tibble of indexed conditons including a Cell State ID
getConditionsIndex <- function(indexFile, conditionsFile, conditionConfig){

    if(!fileDone(indexFile)){
        condIdx <- indexConditions(conditionsFile, conditionConfig)
        openxlsx::write.xlsx(condIdx, indexFile, row.names=F)
    } else {
        log_debug(paste0("Loading previously indexed conditions from file ",indexFile))
        condIdx <- read.xlsx(indexFile, 1, check.names=F) %>% as_tibble()
    }
    condIdx
}

#' Load a single Halo megatable from RDA file
#' 
#' Load all halo data from a single RDA file. Only necessary columns are kept
#' in order to minimize memory. Only positive marker data is kept and any cells NOT
#' positive for DAPI are removed. Optionally, all data marked EXCLUDE are also
#' removed. Finally, column names are adjusted for consistency throughout pipeline, 
#' specifically 'Sample' is changed to 'CellDive_ID' and 'SPOT' is changed to 'FOV_number'.
#' 
#' @param file               full path to RDA file
#' @param filterExclusions   logical; when set to TRUE, any row with text in column EXCLUDE
#'                           will be filtered out
#' @param controlMarker      name of marker used to identify individual, usable cells
#'
#' @return  filtered & formatted Halo data
loadHaloDataFile <- function(file, filterExclusions = FALSE, controlMarker = "DAPI"){

    dat <- readRDS(file)
    if(filterExclusions){
        dat <- dat %>% filter(EXCLUDE == "")
    }

    dat <- dat %>%
           filter(ValueType == "Positive") %>%
           select(Sample, SubSample, SPOT, UUID, Marker, dplyr::matches("XM|YM"), Value)

    ctrlNeg <- dat$UUID[which(dat$Marker == controlMarker && dat$Value == 0)]
    if(length(ctrlNeg) > 0){
        log_warn(paste0(length(ctrlNeg), " ", controlMarker, " negative cells being filtered out."))
        dat <- dat %>% filter(!UUID %in% ctrlNeg)
    }
    dat %>% select(CellDive_ID = Sample, FOV_number = SPOT, everything())

}

#' Load all halo data files into a single table
#' 
#' Load all halo data files into a single table
#'
#' @param dataDir            directory containing ONLY halo RDA files, all of which will be loaded
#' @param dataFiles          vector of specific file paths, possibly a subset of files in a directory
#' @param nThreads           number of threads
#' @param filterExclusions   logical; when set to TRUE, any row with text in column EXCLUDE
#'                           will be filtered out
#' @param controlMarker      name of marker used to identify individual, usable cells
#'
#' @return a single table of filtered & formatted Halo data
loadAllHaloData <- function(dataDir = NULL, dataFiles = NULL, nThreads = 1, 
                            filterExclusions = TRUE, controlMarker = "DAPI"){
    allDat <- NULL
    log_info("Loading all pre-processed halo data...")

    if(is.null(dataFiles) && !is.null(dataDir)){
        dataFiles <- file.path(dataDir, dir(dataDir))
    }

    cl <- makeCluster(nThreads, type = "FORK", outfile = "")
    clusterExport(cl, c('filterExclusions', 'controlMarker'), envir = environment())
    allDat <- parLapply(cl, dataFiles, function(df){
                    log_info(paste0("Reading file ",df))
                    loadHaloDataFile(df, 
                                     filterExclusions = filterExclusions, 
                                     controlMarker = controlMarker)
              })
    allDat <- bind_rows(allDat)

    stopCluster(cl)
    return(allDat)
}

#' Read, format and assemble all patient/sample/fov meta data
#' 
#' Read XLSX files for Patients, Samples and FOVs. Generate unique identifiers
#' for each sample in the form [Patient_ID]_[Sample_number] and for each FOV
#' in the form [Patient_ID]_[Sample_number]_[FOV_number]. Join and flatten all
#' data into a single complete table and create a separate table containing a 
#' map of all identifiers.
#'
#' @param metaFiles     vector of XLSX meta files, including at minimum *_Patients.xlsx,
#'                      *_Samples.xlsx and *_FOVs.xlsx; see docs for details on the formats
#'                      of these files
#' @param metaDataFile  RDA file that already contains or will contain the list of tables
#'                      returned from this function; if file already exists, data in the
#'                      file are what will be returned; otherwise, data will be compiled
#'                      and saved in this file
#'
#' @return list of five tables: Patients, Samples, FOVs, flat, IDs
loadStudyAnnotations <- function(metaFiles = NULL, metaDataFile = NULL){

    if(fileDone(metaDataFile)){
        log_info(paste0("Reading pre-processed meta data from file ",metaDataFile))
        return(readRDS(metaDataFile))
    }

    ptAnnFile     <- metaFiles[grep("Patients.xlsx",metaFiles)]
    sampleAnnFile <- metaFiles[grep("Samples.xlsx",metaFiles)]
    fovAnnFile    <- metaFiles[grep("FOVs.xlsx",metaFiles)]

    ## Read Patients, Samples, FOVs excel spreadsheets
    log_debug("Reading Patient annotations...")
    if(is.null(ptAnnFile) || length(ptAnnFile) != 1){
        log_error("No *_Patients.xlsx file found.")
    }
    pts <- as_tibble(xlsx::read.xlsx(ptAnnFile, 1, check.names = FALSE))

    log_debug("Reading Sample annotations...")
    if(is.null(sampleAnnFile) || length(sampleAnnFile) != 1){
        log_error("No *_Samples.xlsx file found.")
    }
    samples <- as_tibble(xlsx::read.xlsx(sampleAnnFile, 1, check.names = FALSE))

    log_debug("Reading FOV annotations...")
    if(is.null(fovAnnFile) || length(fovAnnFile) != 1){
        log_error("No *_FOVs.xlsx file found.")
    }
    fovs <- as_tibble(xlsx::read.xlsx(fovAnnFile, 1, check.names = FALSE))

    ##
    ### FORMAT IDENTIFIERS
    ##
    ptd <- paste0("%0",max(nchar(as.character(samples$Patient_ID))),"d")
    snd <- paste0("%0",max(nchar(as.character(samples$Sample_number))),"d")
    fnd  <- paste0("%0",max(nchar(as.character(fovs$FOV_number))),"d")

    ##
    ### JOIN ALL DATA INTO A SINGLE FLATTENED TABLE
    ##
    flat <- pts %>%
            full_join(samples, by = "Patient_ID") %>%
            full_join(fovs, by = "CellDive_ID") %>%
            mutate(Lesion_ID = paste(sprintf(ptd, Patient_ID),
                                     sprintf(snd, Sample_number),
                                     sep = "_"),
                   Sample_ID = Lesion_ID,
                   FOV_ID = paste(Lesion_ID,
                                  sprintf(fnd, FOV_number),
                                  sep = "_")) %>%
            filter(is.na(FOV_exclusion))

    ##
    ### CREATE A TABLE OF ONLY IDENTIFIERS
    ##
    idMap <- flat %>%
             select(Patient_ID, CellDive_ID, Sample_ID, Lesion_ID, FOV_ID,
                    Sample_number, FOV_number) %>%
             unique()

    ##
    ### RETURN A LIST OF ALL TABLES
    ##
    annot <- list(Patients = pts, Samples = samples, FOVs = fovs, flat = flat, IDs = idMap)

    return(annot)
}


#' Filter conditions index for all conditions of a single analysis type
#' 
#' Given a table of all conditions including column 'AnalysisType', filter
#' for conditions of a single type (e.g., spatial or general)
#' 
#' @param condIdx       table of conditions including column 'AnalysisType'
#' @param analysisType  character string matching a single value in column
#'                      'AnalysisType'
#'
#' @return  subset of condIdx
filterConditionsByAnalysisType <- function(condIdx, analysisType){
    condIdx %>%
    filter(AnalysisType == analysisType) %>%
    select_if(function(x){ !all(is.na(x)) })
}



#' Compile all analyses to be run
#' 
#' Read 'conditions' XLSX file, in which each tab contains a list of conditions to 
#' be analyzed using a specific calculation type. These types include general fractions
#' and densities in addition to spatial analyses 'neighborhood' fractions and 'neightobhood'
#' average counts
#' 
#' @param condFile    XLSX 'conditions' file described in the docs
#' @param condIdx     table of indexed conditions including Cell State IDs
#'
#' @return list of indexed and formatted conditions where each element contains
#'         conditions pertaining to one calculation type
getStatsAnalyses <- function(condFile, condIdx){ #, cellTypes, funcMarkers, funcCombos){
    fracConds  <- as_tibble(xlsx::read.xlsx(condFile, sheetName="Fraction", check.names=F))
    fracConds  <- fracConds %>%
                  left_join(condIdx %>%
                            filterConditionsByAnalysisType("general"),
                            by=intersect(names(fracConds), names(condIdx)))

    densConds  <- as_tibble(xlsx::read.xlsx(condFile, sheetName="Density", check.names=F))
    densConds  <- densConds %>%
                  left_join(condIdx %>%
                            filterConditionsByAnalysisType("general"),
                            by=intersect(names(densConds), names(condIdx)))

    navgcountsConds <- as_tibble(xlsx::read.xlsx(condFile, sheetName="Neighborhood_averagecounts", check.names=F)) %>%
                       mutate(`Center Population A` = Center,
                              `Neighborhood Population A` = Neighborhood) %>%
                       unite("Condition", c("Center Population A", "Neighborhood Population A"), sep="____", remove = FALSE)
    navgcountsConds <- navgcountsConds %>%
                       left_join(condIdx %>%
                                 filterConditionsByAnalysisType("spatial"),
                                 by = intersect(names(navgcountsConds), names(condIdx)))

    nfracsConds     <- as_tibble(xlsx::read.xlsx(condFile, sheetName="Neighborhood_fraction", check.names=F)) %>%
                       mutate(Condition = paste0(`Center Population A`, "____", `Neighborhood Population A`),
                              Population = paste0(`Center Population B`, "____", `Neighborhood Population B`))
    nfracsConds     <- nfracsConds %>%
                       left_join(condIdx %>%
                                 filterConditionsByAnalysisType("spatial"),
                                 by = intersect(names(nfracsConds), names(condIdx)))

    list(fractions  = fracConds,
         densities  = densConds,
         navgcounts = navgcountsConds,
         nfracs     = nfracsConds)
}


#' Convert XLSX table of study 'question' info to a list
#' 
#' Read XLSX file where each row contains info describing a question
#' and how to subset data in order to answer the question. Convert all
#' rows into a nested list
#'
#' @param  questionXLSXfile  XLSX file
#'
#' @return  question data in list form
parseStatsQuestions <- function(questionsXLSXfile){

    qList <- list()
    questions <- as_tibble(read.xlsx(questionsXLSXfile, 1, check.names=FALSE))

    for(x in 1:nrow(questions)){
        grp  <- as.list(questions[x,])
        if(!grp$QuestionNumber %in% names(qList)){
            q <- list(question  = grp$Question,
                      groupVar  = grp$ComparisonVariable,
                      qNum      = grp$QuestionNumber,
                      calcUnit  = grp$CalculationUnit,
                      `Group 1` = list(),
                      `Group 2` = list())
            qList[[grp$QuestionNumber]] <- q
        }
        filt <- grp[5:length(grp)]
        filt <- lapply(filt[!is.na(filt)], function(x) trimws(unlist(strsplit(as.character(x), ";;"))) )
        qList[[grp$QuestionNumber]][[paste0("Group ",grp$Group)]] <- filt
    }
    qList
}


#' Given a directory containing neighborhood data, load only those files
#' that are needed for the analyses to be run
#' 
#' @param ndhdDirs       main directory of neighborhood data including subdirectories 'C2'
#'                       and 'C3' which contain neighborhood files for Cell_type (C2) and
#'                       Subtype (C3) center cells
#' @param nbhdAnalyses   list with two items: 'nfracs' and 'navgcounts' describing all 
#'                       neighborhood analyses to be run (see docs for format of [StudyName]_Conditions.xlsx)
#' @param fmtNbhdDirs    directory of pre-compiled neighborhood data files (one per samples)
##' @param nbhdFile       file of pre-compiled neighborhood data or name of file to which compiled data
##'                       should be saved
#' @param fovs           character vector of FOV IDs to be included in returned result; default: "all"
#' @return complete table of all neighborhood data associated with analyses to be run
loadMacroNeighborhoodData <- function(nbhdDirs, nbhdAnalyses, fmtNbhdDirs, sampAnn, cellDiveID = "All"){

    nameMap <- c("MHCIIpos_macro" = "M1",
                 "MHCIIneg_macro" = "M2")

    fmtFiles <- file.path(fmtNbhdDirs, dir(fmtNbhdDirs))
    cdids <- sampAnn %>% pull(CellDive_ID) %>% unique

    if(tolower(cellDiveID) != "all"){
        cdids <- cellDiveID
        fmtFiles <- fmtFiles[grepl(paste0("^", cellDiveID, "_"), basename(fmtFiles))]
    }

    dat <- tibble()
    
    for(cdid in cdids){

        nbhdFile <- fmtFiles[grepl(paste0("^", cdid, "_"), basename(fmtFiles))]
        fovs <- sampAnn %>% 
                filter(CellDive_ID == cdid) %>%
                pull(FOV_ID)

        if(fileDone(nbhdFile)){
            log_info(paste0("Loading pre-computed neighborhood counts from file ",nbhdFile))
            dat <- dat %>%
                   bind_rows(readRDS(nbhdFile)) 
        } else {
 
            nbhdFiles <- c()
            for(nd in nbhdDirs){ nbhdFiles <- c(nbhdFiles, file.path(nd, dir(nd))) }
            centers <- unique(c(nbhdAnalyses$nfracs$`Center Population A`,
                                nbhdAnalyses$nfracs$`Center Population B`,
                                nbhdAnalyses$navgcounts$Center))
            mainCTs <- unique(unlist(lapply(strsplit(centers,","), function(x){ x[1] })))

            ########
            oldCenters <- centers
            for(nm in names(nameMap)){
                oldCenters <- gsub(paste0("^",nm), nameMap[[nm]], oldCenters)
            }
            ########

            nbhdDat <- tibble()
            for(ct in mainCTs){
                oldC <- ifelse(ct %in% names(nameMap), nameMap[[ct]], ct)

                log_debug("**************************************")
                log_debug(ct)
                log_debug("**************************************")

                nFile <- nbhdFiles[grepl(gsub("/","_",paste0("-",oldC,"____")),nbhdFiles)]
                nDat  <- readRDS(nFile) %>% filter(FOV %in% fovs)

                nbhds <- nbhdAnalyses$navgcounts %>% filter(Center %in% centers) %>% pull(Neighborhood) %>% unique()
                if(length(nbhds) == 0){ stop() }

                ## first get counts for main center cell type (c)
                for(x in 1:length(nbhds)){
                    nbhd <- nbhds[x]
                    oldNbhd <- nbhd
                    cls <- unlist(strsplit(nbhd,","))[1]
                    if(cls %in% names(nameMap)){
                        other <- unlist(strsplit(oldNbhd,","))[-1]
                        oldNbhd <- paste0(c(nameMap[[cls]], other), collapse=",")
                    }

                    log_debug(paste0(x," out of ",length(nbhds),": ", nbhd))
                    tmp <- filterForNeighborhood(nDat, oldNbhd)
                    log_debug(paste0("    dim(tmp): ",paste(dim(tmp), collapse=",")))
                    if(nrow(tmp) == 0){ next }
                    nbhdDat <- nbhdDat %>%
                               bind_rows(tmp %>%
                                         group_by(FOV,C.UUID) %>%
                                         summarize(N.Count = n()) %>%
                                         mutate(CenterCellType = ct, NeighborhoodCellType = nbhd))
                }

            }
            if(!is.null(nbhdDat) && nrow(nbhdDat) > 0){
                saveRDS(nbhdDat, nbhdFile)
            }
            dat <- dat %>% bind_rows(nbhdDat)    
        }
    }

    dat

}


tmeStatistics <- function(tmeFile, analysisList, nbhdAssignmentFile, assignmentLevel = "cell"){ 

    analyses <- analysisList$fractions %>%
                select(`Cell State ID`,
                       `Neighborhood Population A` = Condition,
                       `Neighborhood Population B` = Population) %>%
                mutate(`Center Population A` = "Tumor", `Center Population B` = "Tumor")

    tme <- readRDS(tmeFile)

    ### JOIN POS/NEG ENV ASSIGNMENTS of NEIGHBORHOOD CELLS
    prefix  <- paste0("C.", assignmentLevel, "_")
    tmp     <- readRDS(nbhdAssignmentFile)
    mrkrCol <- names(tmp)[grepl("microEnv", names(tmp))]
    groupVar  <- gsub("microEnv.", prefix, mrkrCol)
    tme <- tme %>%
           left_join(tmp %>%
                     select_at(c("UUID", mrkrCol)) %>%
                     rename(C.UUID = UUID, !!as.name(groupVar) := !!as.name(mrkrCol)),
                     by = "C.UUID")

    ## reduce at least some for now
    tme <- tme %>%
           rename(FOV_ID = FOV) %>%
           select(-Center, -Neighbor, -N.UUID, -Dij.micron, -C2) %>%
           group_by_at(names(.)) %>%
           summarize(Count = n()) %>%
           mutate(PositiveMarkers = N.FullPosMarkerStr, Classifiers = N.Classifiers) ## rename for filtering 
    tmeStats <- list()
    #for(x in 1:nrow(analyses)){

    tmeStats <- lapply(1:nrow(analyses), function(x){

        pop <- unlist(analyses[x,])
        log_debug(paste0("[", pop[3], "]    ", pop[1], "/", pop[2]))
        csid <- pop$`Cell State ID`

        fracs <- getTMEfractions(tme, csid, "Tumor", pop[1], pop[2])
        #fracs <- tibble()
        #for(p in c("A", "B")){
        #    col <- paste0("Neighborhood Population ", p)
        #    popDat <- tme %>% 
        #              filterDataForCondition(pop[[col]], markers) %>%
        #              mutate(`Cell State ID` = csid, 
        #                     `Center Population A` = pop$`Center Population A`,
        #                     `Center Population B` = pop$`Center Population A`,
        #                     !!as.name(col) := pop[[col]]) %>%
        #              group_by_at(names(.)[!grepl("Marker|Count|Classifiers", names(.))]) %>%
        #              summarize(!!as.name(paste0("N.pop.",p,"_Count")) := sum(Count))
        #    if(p == "A"){
        #        fracs <- popDat
        #    } else {
        #        fracs <- fracs %>% 
        #                 full_join(popDat, by = intersect(names(.), names(popDat)))
        #    
        #        naPop <- which(is.na(fracs$`Neighborhood Population A`))
        #        fracs$`Neighborhood Population A`[naPop]        <- pop[["Neighborhood Population A"]] 
        #        fracs$N.pop.A_Count[is.na(fracs$N.pop.A_Count)] <- 0
        #    }
        #}
        #calcUnit <- "C.UUID" 
        #fracs <- fracs %>%
        #         mutate(Fraction = N.pop.A_Count/N.pop.B_Count) %>%
        #         filter(!!as.name(groupVar) %in% c("Pos.Env","Neg.Env"))

        reportTMELogOdds(fracs, groupVar, calcUnit)
    })

}
