#' Set global variable 'cellTypes'
#'
#' Read and parse cell types XLSX file and save table in global table
#' 'cellTypes'. In case of error, catch and print error message, and warn 
#' that variable is not set.
#'
#' @param ctFile   XLSX file containing cell type definitions
#' 
#' @return nothing
setGlobalCellTypes <- function(ctFile){
    log_info("Setting global variable: 'cellTypes'")
    cellTypes <<- tryCatch({
                      getCellTypes(ctFile)
                    }, error = function(e){
                      log_error(e)
                      log_warn("Global variable 'cellTypes' not set.")
                  })
}

#' Set cell type definition table
#'
#' Read and parse cell types XLSX file and save table in tibble.
#' In case of error, catch and print error message, and warn that 
#' variable is not set.
#'
#' @param ctFile   XLSX file containing cell type definitions
#' 
#' @return nothing
setCellTypes <- function(ctFile){
    log_info("Setting cell type definition table...") 
    tryCatch({
        getCellTypes(ctFile)
      }, error = function(e){
        log_error(e)
        log_warn("Cell type definition table NOT set.")
    })
}

#' Set global variable 'markers'
#'
#' Read marker description file and save marker names in global vector 'markers'.
#' In case of error, catch and print error message, and warn that variable is not set.
#'
#' @param mFile   XLSX file containing all marker information, at least Marker_name
#'
#' @return nothing 
setGlobalMarkers <- function(mFile){
    log_info("Setting global variable: 'markers'")
    markers <<- tryCatch({
                    read.xlsx(mFile, 1, check.names = F) %>%
                    pull(Marker_name)
                  }, error = function(e){
                    log_error(e)
                    log_warn("Global variable 'markers' not set.")
                })
}

#' Set vector of all study markers 
#'
#' Read marker description file and save marker names in a vector. 
#' In case of error, catch and print error message, and warn that 
#' variable is not set.
#'
#' @param mFile   XLSX file containing all marker information, at least Marker_name
#'
#' @return nothing 
setMarkers <- function(mFile){
    log_info("Setting vector of all study markers...") 
    tryCatch({
        read.xlsx(mFile, 1, check.names = F) %>%
        pull(Marker_name)
      }, error = function(e){
        log_error(e)
        log_warn("Vector of study markers NOT set.")
     })
}

#' Set global variable 'sampAnn'
#' 
#' Read, parse and flatten all sample meta files and save in global 
#' table 'sampAnn'. If parameter 'id' is provided, filter annotation
#' on CellDive_ID for given value. In case of error, catch and print error 
#' message, and warn that variable is not set.
#'
#' @param   metaFiles   vector of all sample meta data files
#' @param   id          CellDive_ID of sample to keep; default = NULL (all samples kept)
#'
#' @return
setGlobalSampleAnnotation <- function(metaFiles, id = NULL){
    log_info("Setting global variable: 'sampAnn'")
    sampAnn <<- tryCatch({
                    sa <- loadStudyAnnotations(metaFiles)$flat %>%
                          mutate(Sample = FOV_ID) %>%
                          filter(is.na(FOV_exclusion))
                    if(!is.null(id) && tolower(id) != "All"){
                        sa <- sampAnn %>% filter(CellDive_ID == id)
                    }
                    sa
                  }, error = function(e){
                    log_error(e)
                    log_warn("Global variable 'sampAnn' not set.")
                })
}

#' Set flat table of all sample annotation
#' 
#' Read, parse and flatten all sample meta files and save in  
#' table 'sampAnn'. If parameter 'id' is provided, filter annotation
#' on CellDive_ID for given value. In case of error, catch and print error 
#' message, and warn that variable is not set.
#'
#' @param   metaFiles   vector of all sample meta data files
#' @param   id          CellDive_ID of sample to keep; default = NULL (all samples kept)
#'
#' @return
setSampleAnnotation <- function(metaFiles, id = NULL){
    log_info("Setting flattened table of sample annotation...")
    tryCatch({
        sa <- loadStudyAnnotations(metaFiles)$flat %>%
              mutate(Sample = FOV_ID) %>%
              filter(is.na(FOV_exclusion))
        if(!is.null(id) && tolower(id) != "All"){
             sa <- sampAnn %>% filter(CellDive_ID == id)
        }
        sa
      }, error = function(e){
        log_error(e)
        log_warn("Flattened table of sample annotation NOT set.")
    })
}


#' Set global variable 'allQuestions'
#'
#' Read, parse and reformat statistics questions file and save in global
#' list 'allQuestions'. In case of error, catch and print error 
#' message, and warn that variable is not set.
#' 
#' @param  qFile       XLSX file describing all study questions
#' @param  question    QuestionNumber to return; default = NULL (all questions)
#' 
#' @return nothing
setGlobalQuestions <- function(qFile, question = NULL){
    log_info("Setting global variable: 'allQuestions'")
    allQuestions <<- tryCatch({
                         qs <- parseStatsQuestions(qFile)
                         if(!is.null(question)){ qs <- qs[names(qs) == question] }
                         qs
                       }, error = function(e){
                         log_error(e)
                         log_warn("Global variable 'allQuestions' not set.")
                     })
}


#' Set list of study questions
#'
#' Read, parse and reformat statistics questions. In case of error, 
#' catch and print error message, and warn that variable is not set.
#' 
#' @param  qFile       XLSX file describing all study questions
#' @param  question    QuestionNumber to return; default = NULL (all questions)
#' 
#' @return nothing
setQuestions <- function(qFile, question = NULL){
    log_info("Setting list of study questions...")
    tryCatch({
        qs <- parseStatsQuestions(qFile)
        if(!is.null(question) && tolower(question) != 'all'){ qs <- qs[names(qs) == question] }
        qs
    }, error = function(e){
        log_error(e)
        log_warn("Study questions not set.")
    })
}


#' Set global variable 'conds'
#' 
#' Generate and index or read all previously indexed cell states/conditions
#' and save in global table 'conds'. In case of error, catch and print error
#' message and warn that variable is not set.
#'
#' @param cFile              XLSX file containing human-generated cell states/conditions
#' @param cIndex             XLSX file that does or will contain indexed version of cFile
#' @param arrangeAnnotation  list describing how conditions should be sorted
#'
#' @return nothing
setGlobalConds <- function(cFile, cIndex, arrangeAnnotation){
    log_info("Setting global variable: 'conds'")
    conds <<- tryCatch({
                  getConditionsIndex(cIndex, cFile, arrangeAnnotation)
                }, error = function(e){
                  log_error(e)
                  log_warn("Global variable 'conds' not set.")
              })
}


#' Set variable 'conds'
#' 
#' Generate and index or read all previously indexed cell states/conditions
#' and save in global table 'conds'. In case of error, catch and print error
#' message and warn that variable is not set.
#'
#' @param cFile              XLSX file containing human-generated cell states/conditions
#' @param cIndex             XLSX file that does or will contain indexed version of cFile
#' @param arrangeAnnotation  list describing how conditions should be sorted
#'
#' @return nothing
setConditions <- function(cFile, cIndex, arrangeAnnotation){
    log_info("Setting table of indexed conditions...")
    tryCatch({
        getConditionsIndex(cIndex, cFile, arrangeAnnotation)
      }, error = function(e){
        log_error(e)
      log_warn("Table of indexed conditions NOT set.")
    })
}


#' Set global variable 'analysisList'
#'
#' Divide indexed conditions into list according to calculation type
#' and save list in global list 'analysisList'. 
#'
#' @param cFile   XLSX file containing human-generated cell states/conditions
#' @param cTbl    table containing all indexed conditions
#'
#' @return nothing
setGlobalAnalyses <- function(cFile, cTbl){
    log_info("Setting global variable: 'analysisList'")
    analysisList <<- tryCatch({
                         getStatsAnalyses(cFile, cTbl)
                       }, error = function(e){
                         log_error(e)
                         log_warn("Global variable 'analysisList' not set.")
                     })
}

#' Set variable 'analysisList'
#'
#' Divide indexed conditions into list according to calculation type
#' and save list in global list 'analysisList'. 
#'
#' @param cFile   XLSX file containing human-generated cell states/conditions
#' @param cTbl    table containing all indexed conditions
#'
#' @return nothing
setAnalyses <- function(cFile, cTbl){
    log_info("Setting list of analyses by calculation type...") 
    tryCatch({
        getStatsAnalyses(cFile, cTbl)
      }, error = function(e){
        log_error(e)
        log_warn("List of all analyses NOT set.")
    })
}


#' Set global variable 'annCells'
#' 
#' Load all annotated cells tables and save in global table 'annCells'. 
#' In case of error, catch and print error message and warn that variable is
#' not set.
#' 
#' @param cellDir     directory containing all annotated cells tables
#' @param fovAreaDir  directory containing FOV area files
#' @param bandDir     directory containing all infiltration band assignments and areas 
#' @param studyAnn    flattened table of all study annotation
#' @param id          CellDive_ID of sample of interest; default = NULL, all samples
#'                    loaded
#'
#' @return
setGlobalAnnCells <- function(cellDir, fovAreaDir, bandDir, studyAnn, id = NULL){
    log_info("Setting global variable: annCells")
    annCells <<- tryCatch({
                     loadAnnotatedCells(studyAnn, cellDir,  
                                        bandDir = bandDir,  
                                        fovAreaDir = fovAreaDir,  
                                        cellDiveID = id) 
                   }, error = function(e){
                         log_error(e)
                         log_warn("Global variable 'annCells' not set.")
                 })
}

#' Set variable 'annCells'
#' 
#' Load all annotated cells tables and save in tibble 'annCells'. 
#' In case of error, catch and print error message and warn that variable is
#' not set.
#' 
#' @param cellDir     directory containing all annotated cells tables
#' @param fovAreaDir  directory containing FOV area files
#' @param bandDir     directory containing all infiltration band assignments and areas 
#' @param studyAnn    flattened table of all study annotation
#' @param id          CellDive_ID of sample of interest; default = NULL, all samples
#'                    loaded
#'
#' @return
setAnnCells <- function(cellDir, fovAreaDir, bandDir, studyAnn, id = NULL){
    log_info("Setting study annotated cells table...")
    tryCatch({
        loadAnnotatedCells(studyAnn, cellDir,
                           bandDir = bandDir,
                           fovAreaDir = fovAreaDir,
                           cellDiveID = id)
      }, error = function(e){
        log_error(e)
        log_warn("Study table of annotated cells not set.")
    })
}


#' Set global variable 'nbhdCounts'
#'
#' Read all macrophage neighborhood files and save counts in global table
#' 'nbhdCounts'. In case of error, catch and print error message and warn
#' that variable is not set.
#'
#' @param nbhdDir        directory containing subfolders "C2" and "C3", each
#'                       containing files of cell to cell distances for 
#'                       cell pairs within 30 microns of each other
#' @param nbhdCountsDir  directory that already contains or will contain
#'                       neighborhood cell type counts for each center cell
#'                       of a specific type (e.g., counts of T cells within
#'                       30 microns of MHCIIpos_macro cells, per FOV)
#' @param cells          annotated cells table
#' @param id             CellDive_ID of sample to keep
#' @param analyses       list of all analyses
#' @param studyAnn       flattened table of all study meta data
#' 
#' @return nothing
setGlobalNbhdCounts <- function(nbhdDir = NULL, nbhdCountsDir = NULL, cells = NULL,
                                id = NULL, analyses = NULL, studyAnn = NULL){

    log_info("Setting global variable: nbhdCounts")
    nbhdCounts <<- tryCatch({
                       nbhdDirs <- file.path(nbhdDir, c("C2", "C3"))
                       ## make sure counts files exist
                       ncts <- loadMacroNeighborhoodData(nbhdDirs, analyses, nbhdCountsDir, 
                                                         studyAnn, cellDiveID = id)
                       formatNeighborhoodCounts(nbhdCountsDir, cells, cellDiveID = id)
                     }, error = function(e){
                         log_error(e)
                         log_warn("Global variable 'nbhdCounts' not set.")
                   })
}


#' Set variable 'nbhdCounts'
#'
#' Read all macrophage neighborhood files and save counts in table
#' 'nbhdCounts'. In case of error, catch and print error message and warn
#' that variable is not set.
#'
#' @param nbhdDir        directory containing subfolders "C2" and "C3", each
#'                       containing files of cell to cell distances for 
#'                       cell pairs within 30 microns of each other
#' @param nbhdCountsDir  directory that already contains or will contain
#'                       neighborhood cell type counts for each center cell
#'                       of a specific type (e.g., counts of T cells within
#'                       30 microns of MHCIIpos_macro cells, per FOV)
#' @param cells          annotated cells table
#' @param id             CellDive_ID of sample to keep
#' @param analyses       list of all analyses
#' @param studyAnn       flattened table of all study meta data
#' 
#' @return nothing
setNbhdCounts <- function(nbhdDir = NULL, nbhdCountsDir = NULL, cells = NULL,
                                id = NULL, analyses = NULL, studyAnn = NULL){

    log_info("Setting table of macrophage neighborhood counts...")
    tryCatch({
        nbhdDirs <- file.path(nbhdDir, c("C2", "C3"))
        ncts <- loadMacroNeighborhoodData(nbhdDirs, analyses, nbhdCountsDir,
                                          studyAnn, cellDiveID = id)
        formatNeighborhoodCounts(nbhdCountsDir, cells, cellDiveID = id)
      }, error = function(e){
        log_error(e)
        log_warn("Table of neighborhood counts NOT set.")
    })
}

#' Set global variables 'tumorNbhd' and 'tumorNbhdCells'
#'
#' Read tumor neighborhood file and save all data in global table 'tumorNbhd'. 
#' For convenience, save unique list of neighborhood cell UUIDs in global 
#' vector 'tumorNbhdCells'. In case of error, catch and print error message and warn
#' that variable is not set.
#' 
#' @param nbhdDir    directory containing subfolders "C2" and "C3", each
#'                   containing files of cell to cell distances for 
#'                   cell pairs within 30 microns of each other
#' @param studyAnn   flattened table of all study meta data
#' @param id         CellDive_ID of sample to keep
#' @param questions  global questions list, used to determine whether it is 
#'                   necessary to load tumor neighborhood data
#'
#' @return
setGlobalTumorNeighborhood <- function(nbhdDir, studyAnn, id = NULL, questions = NULL){

    if(!is.null(questions)){
        if(!any(tolower(unlist(questions)) == "neighborhood")){
            log_warn("No questions are restricted to tumor neighborhood cell regions. Global variables 'tumorNbhd' and 'tumorNbhdCells' not set.")
            return()                
        }
    }

    log_info("Setting global variables: tumorNbhd and tumorNbhdCells")

    tryCatch({

        nDir  <- file.path(nbhdDir, "C2")
        tFile <- file.path(nDir, dir(nDir)[grepl("Tumor", dir(nDir))])

        tumorNbhd <<- readRDS(tFile)

        if(!is.null(id)){
            fovs <- studyAnn %>% pull(FOV_ID) %>% unique()
            tumorNbhd <<- tumorNbhd %>% filter(FOV %in% fovs)
        }

        if(nrow(tumorNbhd) > 0){
            tumorNbhdCells <<- tumorNbhd %>% pull(N.UUID) %>% unique()
        }

      }, error = function(e){
        log_error(e)
        log_warn("Global variables 'tumorNbhd' and 'tumorNbhdCells' not set.")
    })
}


#' Set study variables 'tumorNbhd' and 'tumorNbhdCells'
#'
#' Read tumor neighborhood file and save all data in table 'tumorNbhd'. 
#' For convenience, save unique list of neighborhood cell UUIDs in  
#' vector 'tumorNbhdCells'. In case of error, catch and print error 
#' message and warn that variable is not set.
#' 
#' @param nbhdDir    directory containing subfolders "C2" and "C3", each
#'                   containing files of cell to cell distances for 
#'                   cell pairs within 30 microns of each other
#' @param studyAnn   flattened table of all study meta data
#' @param id         CellDive_ID of sample to keep
#' @param questions  study questions list, used to determine whether it is 
#'                   necessary to load tumor neighborhood data
#'
#' @return
setTumorNeighborhood <- function(nbhdDir, studyAnn, id = NULL, questions = NULL){

    if(!is.null(questions)){
        if(!any(tolower(unlist(questions)) == "neighborhood")){
            log_warn("No questions are restricted to tumor neighborhood cell regions. Study variables 'tumorNbhd' and 'tumorNbhdCells' not set.")
            return()
        }
    }

    log_info("Setting study variables: tumorNbhd and tumorNbhdCells")
    tumorNbhd <- NULL
    tumorNbhdCells <- NULL

    tryCatch({

        nDir  <- file.path(nbhdDir, "C2")
        tFile <- file.path(nDir, dir(nDir)[grepl("Tumor", dir(nDir))])

        log_info(paste0("Loading cells in tumor neighborhood from file", tFile))
        tumorNbhd <- readRDS(tFile)

        if(!is.null(id)){
            fovs <- studyAnn %>% pull(FOV_ID) %>% unique()
            tumorNbhd <- tumorNbhd %>% filter(FOV %in% fovs)
        }

        if(nrow(tumorNbhd) > 0){
            tumorNbhdCells <- tumorNbhd %>% pull(N.UUID) %>% unique()
        }

        list('tumorNbhd' = tumorNbhd, 'tumorNbhdCells' = tumorNbhdCells)

      }, error = function(e){
        log_error(e)
        log_warn("Tumor neighborhood table and vector of cells in tumor neighborhood not set.")
    })
}


#' Add tumor microenvironment assignments to global annotated cells table
#' 
#' Read each TME file, and join to global variable 'annCells'. Depenging on
#' assignment level specified, prepend level to column name (e.g., cell_MHCI). 
#' In case of error, catch and print error message and warn that variable 
#' is not set.
#'
#' @param tmeDir           directory containing tumor microenvironment assignment
#'                         files, each containing one column named 'microEnv.*". 
#'                         Column in final table will be named according to assignment
#'                         level and microenvironment (e.g., cell_MHCI)
#' @param assignmentLevel  'sample' or 'cell', the level at which TME assignments were made
#' @param questions        global list of questions
#'
#' @return nothing
globalAddTME <- function(tmeDir, assignmentLevel = "sample", questions = NULL){

    if(!is.null(questions)){
        colPtrn <- paste0(assignmentLevel, "_")
        if(!any(grepl(colPtrn, names(unlist(questions))))){
            log_warn("No questions involve TME on ", assignmentLevel, " level. No TME columns added to 'annCells'.")
            return()
        }
    }

    tryCatch({

        log_info(paste0("Adding ", assignmentLevel, " level TME assignments to global variable 'annCells'"))

        files <- file.path(tmeDir, dir(tmeDir))

        for(fl in files){
            tme <- readRDS(fl)
            mrkrCol <- names(tme)[grepl("microEnv", names(tme))]
            mrkr <- gsub("microEnv\\.","",mrkrCol)
            annCells <<- annCells %>%
                         left_join(tme %>%
                                   select_at(c("UUID", mrkrCol)) %>%
                                   rename_at(vars(mrkrCol),
                                             list(~(. = paste0(assignmentLevel, "_",mrkr)))),
                                   by = "UUID")
        }
      }, error = function(e){
        log_error(e)
        log_warn(paste0("TME assignments on ", assignmentLevel, " level not added to 'annCells'"))
    })

}


#' Add tumor microenvironment assignments to annotated cells table
#' 
#' Read each TME file, and join to study variable 'annCells'. Depenging on
#' assignment level specified, prepend level to column name (e.g., cell_MHCI). 
#' In case of error, catch and print error message and warn that variable 
#' is not set.
#'
#' @param tmeDir           directory containing tumor microenvironment assignment
#'                         files, each containing one column named 'microEnv.*". 
#'                         Column in final table will be named according to assignment
#'                         level and microenvironment (e.g., cell_MHCI)
#' @param annCells         tibble of annotated cells to which TME column should be added
#' @param assignmentLevel  'sample' or 'cell', the level at which TME assignments were made
#' @param questions        study list of questions
#'
#' @return nothing
addTME <- function(annCells, tmeDir, assignmentLevel = "sample", questions = NULL){

    if(!is.null(questions)){
        colPtrn <- paste0(assignmentLevel, "_")
        if(!any(grepl(colPtrn, names(unlist(questions))))){
            log_warn("No questions involve TME on ", assignmentLevel, " level. No TME columns added to 'annCells'.")
            return(annCells)
        }
    }

    tryCatch({

        log_info(paste0("Adding ", assignmentLevel, " level TME assignments to study variable 'annCells'"))

        files <- file.path(tmeDir, dir(tmeDir))

        for(fl in files){
            tme <- readRDS(fl)
            mrkrCol <- names(tme)[grepl("microEnv", names(tme))]
            mrkr <- gsub("microEnv\\.","",mrkrCol)
            annCells <- annCells %>%
                        left_join(tme %>%
                        select_at(c("UUID", mrkrCol)) %>%
                        rename_at(vars(mrkrCol),
                        list(~(. = paste0(assignmentLevel, "_",mrkr)))),
                      by = "UUID")
        }

        annCells

      }, error = function(e){
        log_error(e)
        log_warn(paste0("TME assignments on ", assignmentLevel, " level not added to 'annCells'"))
    })

}


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
#' @param cellsInTumorNeighborhood logical; when TRUE, read tumor neighborhood file; save table in 'tumorNbhd'
#                                  and save UUIDs of all neighborhood cells in global vector 'tumorNbhdCells'
#' @param tmeSampleStatus         logical; when TRUE, add to 'annCells' table columns for all categorized tumor
#'                                microenvironments (e.g., MHCI+/- TME) available in directory referenced by
#'                                'tme_by_sample_dir' key in config
#' @param tmeCellStatus           logical; when TRUE, add to 'annCells' table columns for all categorized tumor
#'                                microenvironments (e.g., MHCI+/- TME) available in directory referenced by
#'                                'tme_by_cell_dir' key in config
#'
#' @return nothing
loadGlobalStudyData <- function(config,
                                all = FALSE,
                                cellTypeDefinitions = TRUE,
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
        cellTypeDefinitions <- markerList <- analyses <- conditions <- sampleAnnotation <- TRUE
        annotatedCells <- questions <- neighborhoodCounts <- cellsInTumorNeighborhood <- TRUE
        #tmeSampleStatus <- TRUE 
        tmeCellStatus <- TRUE
    }

    cellTypes       <<- NULL
    markers         <<- NULL
    allQuestions    <<- NULL
    conds           <<- NULL
    analysisList    <<- NULL
    sampAnn         <<- NULL
    annCells        <<- NULL
    nbhdCounts      <<- NULL
    tumorNbhd       <<- NULL
    tumorNbhdCells  <<- NULL
    

    if(is.null(config$cell_dive_id)){
        config$cell_dive_id <- "All"
    }

    if(!is.null(config$debug) && (config$debug == "yes" || config$debug)){
        log_threshold(DEBUG)
    }

    if(cellTypeDefinitions){
        setGlobalCellTypes(getFiles(path = config$meta_dir, 
                                    files = config$meta_files, 
                                    pattern = "CellTypes.xlsx"))
    }

    if(markerList){
        setGlobalMarkers(getFiles(path = config$meta_dir, 
                                  files = config$meta_files, 
                                  pattern = "Markers.xlsx"))
    }

    if(questions){
        setGlobalQuestions(config$statistics_questions_file, question = config$question)
    }

    if(conditions || analyses){
        setGlobalConds(config$statistics_conditions_file, 
                       config$statistics_conditions_index,
                       config$arrange_annotation)
    }

    if(analyses){
        setGlobalAnalyses(config$statistics_conditions_file, conds)
    }

    if(sampleAnnotation || annotatedCells){
        metaFiles <- getFiles(path = config$meta_dir, pattern = ".xlsx")
        setGlobalSampleAnnotation(metaFiles)
    }

    if(annotatedCells || neighborhoodCounts || tmeSampleStatus || tmeCellStatus){
        setGlobalAnnCells(config$cell_data_dir, 
                          config$fov_area_dir, 
                          config$band_dir, 
                          sampAnn, 
                          id = config$cell_dive_id)
    }

    if(neighborhoodCounts){
        setGlobalNbhdCounts(nbhdDir = config$neighborhood_dir, 
                            nbhdCountsDir = config$neighborhood_counts_dir, 
                            cells = annCells, 
                            id = config$cell_dive_id, 
                            analyses = analysisList, 
                            studyAnn = sampAnn)
    }

    if(cellsInTumorNeighborhood | tmeSampleStatus | tmeCellStatus){
        setGlobalTumorNeighborhood(config$neighborhood_dir, sampAnn, 
                                   id = config$cell_state_id,
                                   questions = allQuestions)
    }

    if(tmeSampleStatus){
        globalAddTME(config$tme_by_sample_dir, 
                     assignmentLevel = "sample", 
                     questions = allQuestions)
    }

    if(tmeCellStatus){
        globalAddTME(config$tme_by_sample_dir, 
                     assignmentLevel = "cell", 
                     questions = allQuestions)
    }
}




#' Load all study data 
#'
#' Given just the study configuration in list form, load into a list all 
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
#' @param cellsInTumorNeighborhood logical; when TRUE, read tumor neighborhood file; save table in 'tumorNbhd'
#                                  and save UUIDs of all neighborhood cells in global vector 'tumorNbhdCells'
#' @param tmeSampleStatus         logical; when TRUE, add to 'annCells' table columns for all categorized tumor
#'                                microenvironments (e.g., MHCI+/- TME) available in directory referenced by
#'                                'tme_by_sample_dir' key in config
#' @param tmeCellStatus           logical; when TRUE, add to 'annCells' table columns for all categorized tumor
#'                                microenvironments (e.g., MHCI+/- TME) available in directory referenced by
#'                                'tme_by_cell_dir' key in config
#'
#' @return nothing
loadStudyData <- function(config,
                          all = FALSE,
                          cellTypeDefinitions = TRUE,
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
        cellTypeDefinitions <- markerList <- analyses <- conditions <- sampleAnnotation <- TRUE
        annotatedCells <- questions <- neighborhoodCounts <- cellsInTumorNeighborhood <- TRUE
        #tmeSampleStatus <- TRUE 
        tmeCellStatus <- TRUE
    }

    stDat <- list()

    if(is.null(config$cell_dive_id)){
        config$cell_dive_id <- "All"
    }

    if(!is.null(config$debug) && (config$debug == "yes" || config$debug)){
        log_threshold(DEBUG)
    }

    if(cellTypeDefinitions){
        stDat$cellTypes <- setCellTypes(getFiles(path = config$meta_dir,
                                        files = config$meta_files,
                                        pattern = "CellTypes.xlsx"))
    }

    if(markerList){
        stDat$markers <- setMarkers(getFiles(path = config$meta_dir,
                                    files = config$meta_files,
                                    pattern = "Markers.xlsx"))
    }

    if(questions){
        stDat$allQuestions <- setQuestions(config$statistics_questions_file, 
                                           question = config$question)
    }

    if(conditions || analyses){
        stDat$conds <- setConditions(config$statistics_conditions_file,
                                     config$statistics_conditions_index,
                                     config$arrange_annotation)
    }

    if(analyses){
        stDat$analysisList <- setAnalyses(config$statistics_conditions_file, stDat$conds)
    }

    if(sampleAnnotation || annotatedCells || neighborhoodCounts || tmeCellStatus){
        metaFiles <- getFiles(path = config$meta_dir, pattern = ".xlsx")
        stDat$sampAnn <- setSampleAnnotation(metaFiles)
    }

    if(annotatedCells || neighborhoodCounts || tmeSampleStatus || tmeCellStatus){
        stDat$annCells <- setAnnCells(config$cell_data_dir,
                                      config$fov_area_dir,
                                      config$band_dir,
                                      stDat$sampAnn,
                                      id = config$cell_dive_id)
    }

    if(neighborhoodCounts){
        stDat$nbhdCounts <- setNbhdCounts(nbhdDir = config$neighborhood_dir,
                                          nbhdCountsDir = config$neighborhood_counts_dir,
                                          cells = stDat$annCells,
                                          id = config$cell_dive_id,
                                          analyses = stDat$analysisList,
                                          studyAnn = stDat$sampAnn)
    }

    if(cellsInTumorNeighborhood || tmeSampleStatus || tmeCellStatus){
        stDat <- c(stDat, 
                   setTumorNeighborhood(config$neighborhood_dir, 
                                        stDat$sampAnn,
                                        id = config$cell_state_id,
                                        questions = stDat$allQuestions))
    }

    if(tmeSampleStatus){
        stDat$annCells <- stDat$annCells %>%
                          addTME(config$tme_by_sample_dir,
                                 assignmentLevel = "sample",
                                 questions = stDat$allQuestions)
    }

    if(tmeCellStatus){
        stDat$annCells <- stDat$annCells %>%
                          addTME(config$tme_by_cell_dir,
                                 assignmentLevel = "cell",
                                 questions = stDat$allQuestions)
    }

    stDat
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
    checkFilesFound(annFiles, annDir, "annotated cells")

    if(!is.null(bandDir)){
        baFiles  <- file.path(bandDir, dir(bandDir))
        checkFilesFound(baFiles, bandDir, "infiltration band info")
    }
    if(!is.null(fovAreaDir)){
        faFiles  <- file.path(fovAreaDir, dir(fovAreaDir))
        checkFilesFound(faFiles, fovAreaDir, "FOV area")
    }

    if(!is.null(cellDiveID) && tolower(cellDiveID[1]) != "all"){
        cdidPat  <- paste0("^", cellDiveID, "_")
        annFiles <- annFiles[grepl(cdidPat, basename(annFiles))]
        checkFilesFound(annFiles, annDir, "annotated cells")
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
            ptrn <- paste0("^", id, "_")
            if(!is.null(faFiles)){ faFile <- faFiles[grep(ptrn, basename(faFiles))] }
            if(!is.null(faFile) && length(faFile) == 1 && file.exists(faFile)){
                log_debug(paste0("joining FOV areas from file: ", faFile))
                fa <- readRDS(faFile)
                sampCells <- sampCells %>%
                             left_join(fa, by = intersect(names(sampCells), names(fa)))
            }

            baFile <- NULL
            if(!is.null(baFiles)){ baFile <- baFiles[grep(ptrn, basename(baFiles))] }
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
              }) %>%
              bind_rows()

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
    annot <- list(Patients = pts, 
                  Samples = samples, 
                  FOVs = fovs, 
                  flat = flat, 
                  IDs = idMap)

    return(annot)
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

    fracConds  <- xlsx::read.xlsx(condFile, sheetName="Fraction", check.names=F) %>%
                  as_tibble() %>%
                  left_join(condIdx %>%
                            filterConditionsByAnalysisType("general"),
                            by=intersect(names(.), names(condIdx)))

    densConds  <- xlsx::read.xlsx(condFile, sheetName="Density", check.names=F) %>%
                  as_tibble() %>%
                  left_join(condIdx %>%
                            filterConditionsByAnalysisType("general"),
                            by=intersect(names(.), names(condIdx)))

    navgcountsConds <- xlsx::read.xlsx(condFile, sheetName="Neighborhood_averagecounts", check.names=F) %>%
                       as_tibble() %>%
                       mutate(`Center Population A` = Center,
                              `Neighborhood Population A` = Neighborhood) %>%
                       unite("Condition", 
                             c("Center Population A", "Neighborhood Population A"), 
                             sep="____", 
                             remove = FALSE) %>%
                       left_join(condIdx %>%
                                 filterConditionsByAnalysisType("spatial"),
                                 by = intersect(names(.), names(condIdx)))

    nfracsConds     <- xlsx::read.xlsx(condFile, sheetName="Neighborhood_fraction", check.names=F) %>%
                       as_tibble() %>%
                       mutate(Condition = paste0(`Center Population A`, "____", `Neighborhood Population A`),
                              Population = paste0(`Center Population B`, "____", `Neighborhood Population B`)) %>%
                       left_join(condIdx %>%
                                 filterConditionsByAnalysisType("spatial"),
                                 by = intersect(names(.), names(condIdx)))

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
#' @param nbhdCountsDir  directory of pre-compiled neighborhood data files (one per samples)
#' @param fovs           character vector of FOV IDs to be included in returned result; default: "all"
#'
#' @return complete table of all neighborhood data associated with analyses to be run
loadMacroNeighborhoodData <- function(nbhdDirs, nbhdAnalyses, nbhdCountsDir, sampAnn, cellDiveID = "All"){

    nameMap <- c("MHCIIpos_macro" = "M1",
                 "MHCIIneg_macro" = "M2")

    fmtFiles <- file.path(nbhdCountsDir, dir(nbhdCountsDir))
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


