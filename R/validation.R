## required columns (constants)
req_cols <- list(patients = c("Patient_ID", "Patient_response"),
                 samples  = c("CellDive_ID", "Patient_ID", "Sample_number"),
                 fovs     = c("CellDive_ID", "FOV_number", "FOV_exclusion", "FOV_exclusion_stage",
                              "Marker_exclusion", "Num_manual_exclusions", "Num_epidermis_exclusions",
                              "Num_interface_areas"),
                 cell_types = c("Category", "Cell_type", "Subtype",
                                "Pos_markers", "Pos_required", "Neg_markers", "Tag"),
                 markers = c("Marker_name", "Description", "Cell_type"),
                 annotated_cells = c("Patient_ID", "Sample_ID", "FOV_ID"),
                 sample_annotation = c("Sample", "FOV"),
                 questions = c("Question", "QuestionNumber", "ComparisonVariable", "Group"))

#' Check that a table includes all columns required for that dataset
#' 
#' @param cols   vector of cols in the table to be checked
#' @param dataset  name of dataset being checked; values allowed:
#'                      ['patients'|'samples'|'fovs'|'cell_types'|'markers'|
#'                        'annotated_cells'|'sample_annotation']
#' @return vector of error messages (NULL if none found)
checkRequiredCols <- function(cols, dataset){
    errs <- c()
    req <- req_cols[[dataset]]
    missing <- setdiff(req, cols) 
    if(length(missing) > 0){
        for(m in missing){ errs <- c(errs, paste0("Column missing from ",dataset," table: ",m)) }
    }
    return(errs)
}

#' Check for existence of all required meta files
#' 
#' Given entire study configuration, ensure that either {meta_files} or {meta_dir} is defined.
#' If {meta_dir} given and not {meta_files}, populate {meta_files} with all xlsx files in 
#' directory defined by {meta_dir}. Minimum required files are: 
#'    *_CellTypes.xlsx
#'    *_Markers.xlsx
#'    *_Patients.xlsx
#'    *_Samples.xlsx
#'    *_FOVs.xlsx
#'
#' @param sCfg    study config in list format
#' @return nested list of errors, warnings and messages
validateMetaFileList <- function(sCfg){

    log_debug("  validating meta file list..\n")
    errors <- warnings <- messages <- list()

    ## meta_dir/files
    metaFiles <- c()
    if(is.null(sCfg$meta_files)){
        if(is.null(sCfg$meta_dir)){
            error[[length(errors) + 1]] <- "Must specify meta_dir OR/AND meta_files"
        } else if(!dir.exists(sCfg$meta_dir) || length(dir(sCfg$meta_dir)) == 0) {
            error[[length(errors) + 1]] <- "meta_dir must exists and must not be empty"
        } else {
            metaFiles <- file.path(sCfg$meta_dir, dir(sCfg$meta_dir))
        }
    } else {
        metaFiles <- sCfg$meta_files
    }

    req <- c("_CellTypes\\.xlsx",
             "_Markers\\.xlsx",
             "_Patients\\.xlsx",
             "_Samples\\.xlsx",
             "_FOVs\\.xlsx")

    for(fs in req){
        if(!any(grepl(fs, metaFiles))){
            errors[[length(errors) + 1]] <- paste0("Can not find required file *",fs," in list of meta files.")
        }
    }

    return(list(errs=errors, wrns=warnings, msgs=messages))
}

#' Check for access to halo boundary information
#' 
#' Given study config if list format, check for either a pre-processed
#' *.rda file or a directory of XML files directly from Halo; if a directory is specified,
#' it must contain one or more subdirectories, each containing XML files for one type of boundary.
#' 
#' @param sCfg    study config in list format
#' @return nested list of errors, warnings and messages
validateHaloBoundaryFiles <- function(sCfg){
    
    log_debug("  validating halo boundary files..\n")
    errors <- warnings <- messages <- list()

    hbf <- sCfg$halo_boundaries_file
    hbd <- sCfg$halo_boundaries_xml_dir

    if(all(is.null(c(hbf,hbd))) | 
        (is.null(hbd) & (!file.exists(hbf) | file.size(hbf) == 0)) | 
           (!is.null(hbd) & (!dir.exists(hbd) | file.size(hbd) == 0))){
        errors[[length(errors) + 1]] <- "Must provide either valid halo_boundaries_file or valid halo_boundaries_xml_dir"
    }

    return(list(errs = errors, wrns = warnings, msgs = messages))

}


#' List files that already exist and will be used as input and files that are to be 
#' written by pipeline
#'
#' Given entire study config, for each *_file paramater print info on whether the file
#' exists and will be used as input or does not yet exist (or is empty) and will be
#' generated during analysis.
#' 
#' @param   sCfg    study config in list format
#' @return nested list of errors, warnings and messages
validateProvidedFiles <- function(sCfg){

    log_debug("  validating provided files...\n")
    warnings <- existing <- non <- list()
    fileFields <- names(sCfg)[grepl("_file$",names(sCfg))]
    for(ff in fileFields){
        if(fileDone(sCfg[[ff]])){
            existing[[ff]] <- sCfg[[ff]]
        } else {
            non[[ff]] <- sCfg[[ff]]
        }
    }
    if(length(existing) > 0){
        warnings <- as.list("The following files already exist and will be used (NOT OVERWRITTEN) in associated analyses: ")
        warnings <- c(warnings, paste0("\t", names(existing), ": ", existing))
        warnings[[length(warnings) + 1]] <- "\n"
    }
    if(length(non) > 0){
        warnings[[length(warnings) + 1]] <- "The following files do NOT exist or are empty and contents will be generated from scratch:"
        warnings <- c(warnings, paste0("\t", names(non), ": ",non))
        warnings[[length(warnings) + 1]] <- "\n"
    }

    return(list(errs=list(), wrns=warnings, msgs=list()))
}

#' Indicate duplicate values within a vector, usually a column from a dataframe
#' 
#' Given a vector of values, print info on which values are found multiple times.
#' 
#' @param vals  vector of values
#' @param label  label (usually column name) to include in message
#' @param rtrn   ['message'|'values'] what should be returned
#' @return  character string with message of duplicates or vector of
#'          duplicated values (NULL if none found)
duplicatesWithinCol <- function(vals, label, rtrn='message'){
    e <- NULL
    if(length(vals[duplicated(vals)]) > 0){
        e <- switch(rtrn,
                    'message' = paste0("Duplicated ",label,"(s): ",
                        paste(vals[duplicated(vals)], collapse=", ")),
                    'values' = vals[duplicated(vals)],
                    NULL)
    }
    e
}


#' Indicate duplicates across two columns in a data frame
#' 
#' Given a data frame of two columns, print info on values that occur in both
#' columns.
#' 
#' @param df  data frame of two columns
#' @return  character string with message of values in both columns (NULL if none found)
duplicatesAcrossCols <- function(df){
    e <- NULL
    if(any(df[,1] %in% df[,2])){
        e <- paste0(paste(df[,1][df[,1] %in% df[,2]], collapse=", "),
                    " in both ",names(df)[1]," and ",names(df)[2]," columns.")
    }
    e
}

#' Check for conflicts between marker combinations of duplicate classes
#' 
#' @param df  data frame of cell types with a single identifier and at minimum, 
#'            columns Pos_markers and Neg_markers
#' @return message indicating that there is a marker combo conflict between 
#'         cell types with the same identifier
duplicateClassConflicts <- function(df){
    pos <- unique(unlist(strsplit(df %>% pull(Pos_markers), ",")))
    neg <- unique(unlist(strsplit(df %>% pull(Neg_markers), ",")))
    if(length(intersect(pos,neg)) > 0){
        return(paste0("Duplicated Subtypes with conflicting marker combinations: ",
                      paste(stDups, collapse=", ")))
    }
    return(NULL)
}

#' Perform various checks to ensure cell types spreadsheet has been filled out correctly.
#' 
#' Criteria to be checked include:
#'    * Values in Subtype and Tag columns must be unique 
#'    * No value may occur in more than one column
#'    * Classification_type must be one of: 'state' or 'type'
#'    * All markers listed in columns Pos_markers and Neg_markers must be in markers XLSX file
#'    * All Pos_markers+Neg_markers combinations must be unique to avoid ambiguous cell type definitions
validateCellTypes <- function(cellTypesXLSX, markerDescFile){

    errs <- wrns <- list()

    cts   <- read.xlsx(cellTypesXLSX, 1, check.names = F) %>% as_tibble()
    mDesc <- read.xlsx(markerDescFile, 1, check.names = F) %>% as_tibble()
    errs <- checkRequiredCols(names(cts), 'cell_types')
    errs <- c(errs, checkRequiredCols(names(mDesc), 'markers'))

    ## 
    ## unique subtypes and tags
    ##
    for(col in c("Subtype","Tag")){
        colDups <- duplicatesWithinCol(cts %>% pull(col), col, rtrn='values')
        ## if there are subtype duplicates, check that their marker combinations do not conflict
        if(!is.null(colDups)){
            colErrs <- c()
            for(dup in colDups){
                cnfl <- duplicateClassConflicts(cts %>%
                                                filter_at(vars(col), any_vars(. == dup)) %>%
                                                filter(Pos_required == 'all'))
                colErrs <- c(colErrs, cnfl)
            }
            if(length(colErrs) > 0){
                errs <- c(errs, colErrs)
            } else {
                wrns <- c(wrns, paste0("Duplicate but apparently valid ",col,"(s) found: ",
                                       paste(colDups, collapse=", ")))
            }
        }
    }

    ##
    ## no duplicates between Category, Cell_type, Subtype
    ##
    errs[[length(errs) + 1]] <- duplicatesAcrossCols(cts %>% select(Category, Cell_type))
    errs[[length(errs) + 1]] <- duplicatesAcrossCols(cts %>% select(Category, Subtype))
    errs[[length(errs) + 1]] <- duplicatesAcrossCols(cts %>% select(Cell_type, Subtype))

    ##
    ## supported classification types are 'state' and 'type'
    ##
    if(!all(unique(tolower(cts$Classification_type)) %in% c('state','type'))){
        e <- paste0("Unsupported Classification_type(s): ",
                    paste(cts$Classification_type[!cts$Classification_type %in% c('state','type')], collapse=", "))
        errs[[length(errs) + 1]] <- e
    }

    ##
    ## all markers are in MarkerDesc
    ##
    allM <- unique(unlist(strsplit(c(cts$Pos_markers, cts$Neg_markers),",")))
    allM <- allM[!is.na(allM)]
    if(!all(allM %in% mDesc$Marker_name)){
        e <- paste0("Unrecognized markers (markers not in Markers file): ",
                    paste(allM[!allM %in% mDesc$Marker_name], collapse=", "))
        errs[[length(errs) + 1]] <- e
    }

    ##
    ## unique marker combinations (one set of classifications per combination)
    ##
    fullCombos <- cts %>%
                  mutate(Negs = sapply(Neg_markers, function(x){ paste(sort(paste0(unlist(strsplit(x,",")),"-")), collapse=",") }),
                         Pos = sapply(Pos_markers, function(x){ paste(sort(unlist(strsplit(x,","))),collapse=",") }))

    if(any(duplicated(fullCombos))){
        e <- paste0("Duplicate marker combinations: \n",
                    paste0(fullCombos %>%
                           unite("Marker combo", Pos, Negs, sep=",") %>%
                           pull(`Marker combo`) %>%
                           duplicated, collapse="\n"))
        errs[[length(errs) + 1]] <- e
    }

    if(length(wrns) > 0){
        for(w in unlist(wrns)){ log_warn(w) }
    }
    if(length(errs) > 0){
        for(e in unlist(errs)){ log_error(e) }
        return(FALSE)
    }
    return(NULL)
}

#' Check for duplicate Patient_ID
#' 
#' @param patientAnnotationXLSX    XLSX file containing columns Patient_ID and Patient_response
#'                                 at least
#' @return list of errors, if any
validatePatientAnnotation <- function(patientAnnotationXLSX){
    errors <- list()
    pa <- read.xlsx(patientAnnotationXLSX, 1, check.names=F) %>% as_tibble()
    errors <- checkRequiredCols(names(pa), 'patients')
    errors[[length(errors) + 1]] <- duplicatesWithinCol(pa %>% pull(Patient_ID), "Patient_ID")
    return(list(errs=errors))
}

#' Check for errors in sample annotation
#' 
#' @param sampleAnnotationXLSX    XLSX file containing all sample annotation
#' @param patientAnnotationXLSX   XLSX file containing all patient annotation
#' @return  list of errors and warnings 
validateSampleAnnotation <- function(sampleAnnotationXLSX, patientAnnotationXLSX){

    errors <- warnings <- list()

    pa <- read.xlsx(patientAnnotationXLSX, 1, check.names=F) %>% as_tibble()
    sa <- read.xlsx(sampleAnnotationXLSX, 1, check.names=F) %>% as_tibble() %>%
          mutate(Sample_ID = paste0(Patient_ID, "_", Sample_number))
    errors <- checkRequiredCols(names(sa), 'samples')

    ##
    ## unique CellDive_ID and Patient_ID + Sample_number (Sample_ID) combinations
    ##
    errors[[length(errors) + 1]] <- duplicatesWithinCol(sa %>% pull(CellDive_ID), "CellDive_ID")
    errors[[length(errors) + 1]] <- duplicatesWithinCol(sa %>% pull(Sample_ID), "Sample_ID")

    ##
    ## warn of discrepancies between Patient_IDs in sample and patient sheets
    ##
    if(!all(sort(unique(pa$Patient_ID)) == sort(unique(sa$Patient_ID)))){
        inSAnotPA <- sa$Patient_ID[!sa$Patient_ID %in% pa$Patient_ID]
        inPAnotSA <- pa$Patient_ID[!pa$Patient_ID %in% sa$Patient_ID]
        if(length(inSAnotPA) > 0){
            warnings[[length(warnings) + 1]] <-
                     paste0("Patient_ID(s) found in *Samples.xlsx but not *Patients.xlsx: ",
                            paste(inSAnotPA, collapse=", "))
        }
        if(length(inPAnotSA) > 0){
            warnings[[length(warnings) + 1]] <-
                     paste0("Patient_ID(s) found in *Patients.xlsx but not *Samples.xlsx: ",
                            paste(inPAnotSA, collapse=", "))
        }
    }

    return(list(errs=errors, wrns=warnings))
}

#' Check for errors in FOV annotation
#' 
#' @param fovAnnotationXLSX    XLSX file containing all FOV-level annotation
#' @param sampleAnnotationXLSX XLSX file containing all Sample-level annotation
#' @return  list of errors and warnings
validateFOVannotation <- function(fovAnnotationXLSX, sampleAnnotationXLSX){

    errors <- warnings <- messages <- list()

    sa <- as_tibble(read.xlsx(sampleAnnotationXLSX, 1, check.names=F))
    fa <- as_tibble(read.xlsx(fovAnnotationXLSX, 1, check.names=F)) %>%
          mutate(CDID_FOV_combo = paste0(CellDive_ID,"_", FOV_number))
    errors <- checkRequiredCols(names(fa), 'fovs')

    ##
    ## 1:1 CellDiveID to FOV_number (unique fovs within each sample)
    ##
    errors[[length(errors) + 1]] <- duplicatesWithinCol(fa %>% pull(CDID_FOV_combo), "CellDive_ID + FOV combo")

    ##
    ## warn of any differences between sample annotation and fov annotation sample lists
    ##
    if(!all(sort(unique(fa$CellDive_ID)) == sort(unique(sa$CellDive_ID)))){
        inSAnotFA <- sa$CellDive_ID[!sa$CellDive_ID %in% fa$CellDive_ID]
        inFAnotSA <- fa$CellDive_ID[!fa$CellDive_ID %in% sa$CellDive_ID]
        if(length(inSAnotFA) > 0){
            warnings[[length(warnings) + 1]] <-
                     paste0("CellDive_ID(s) found in *SampleAnnotations.xlsx but not *FOVannotations.xlsx: ",
                            paste(inSAnotFA, collapse=", "))
        }
        if(length(inFAnotSA) > 0){
            warnings[[length(warnings) + 1]] <-
                     paste0("CellDive_ID(s) found in *FOVannotations.xlsx but not *SampleAnnotations.xlsx: ",
                            paste(inFAnotSA, collapse=", "))
        }
    }

    return(list(errs=errors, wrns=warnings))
}


#' Check for errors in all study annotation files
#' 
#' @param fovAnnotationXLSX 
#' @param sampleAnnotationXLSX
#' @param patientAnnotationXLSX
#' @param FALSE if any errors found, TRUE if not
validateStudyAnnotation <- function(fovAnnotationXLSX, sampleAnnotationXLSX, patientAnnotationXLSX){
    log_debug("Validating patient annotation...\n")
    resP  <- validatePatientAnnotation(patientAnnotationXLSX)
    log_debug("Validating sample annotation...\n")
    resS <- validateSampleAnnotation(sampleAnnotationXLSX, patientAnnotationXLSX)
    log_debug("Validating FOV annotation...\n")
    resF <- validateFOVannotation(fovAnnotationXLSX, sampleAnnotationXLSX)

    errors   <- c(unlist(resP$errs), unlist(resS$errs), unlist(resF$errs))
    warnings <- c(unlist(resP$wrns), unlist(resS$wrns), unlist(resF$wrns))
    messages <- c(unlist(resP$msgs), unlist(resS$msgs), unlist(resF$msgs))

    cat("\n\n")
    if(length(warnings) > 0){
        for(w in unlist(warnings)){ log_warn(w) }
    }
    if(length(messages) > 0){
        for(m in unlist(messages)){ log_info(m) }
    }
    if(length(errors) > 0){
        for(e in unlist(errors)){ log_error(e) }
        return(FALSE)
    }
    return(TRUE)
}

#' Validate study configuration by checking for all required meta files,
#' boundary files, etc. [TO DO: ADD VALIDATION OF DATA FILES BY CHECKING
#' FOR SAMPLE NAMES IN FILE NAMES
#'
#' Validate study configuration. Log all warnings and errors for user review
#' 
#' @param sCfg  study configuration in form of a named list
#' 
#' @return logical indicating whether configuration is valid and safe for pipeline
validateStudyConfig <- function(sCfg){

    errors <- warnings <- messages <- c()

    allVals <- c(validateMetaFileList,
                 validateHaloBoundaryFiles,
                 validateProvidedFiles)

    for(x in seq(allVals)){
        validation <- allVals[x][[1]]
        res <- validation(sCfg)
        errors <- c(errors, unlist(res$errs))
        warnings <- c(warnings, unlist(res$wrns))
        messages <- c(messages, unlist(res$msgs))
    }

    cat("\n")
    if(length(warnings) > 0){
        for(w in warnings){ log_warn(w) }
    }
    if(length(messages) > 0){
        for(m in messages){ log_info(m) }
    }
    if(length(errors) > 0){
        for(e in errors){ log_error(e) }
        return(FALSE)
    }
    return(TRUE)
}

#' Make sure all exclusions were made correctly and that no
#' excluded FOVs or markers still exist in data
#'
#' Use study annotations containing FOV and marker exclusions to 
#' ensure that all exclusions have been filtered out of data prior
#' to this point.
#' 
#' @param dat       tibble containing pre-filtered data
#' @param studyAnn  tibble of study annotation including columns FOV_exclusion and 
#'                  Marker_exclusion
#' 
#' @param list of two items: one vector of error messages and one of warning messages
validateExclusions <- function(dat, studyAnn){

    errors <- warnings <- messages <- list()

    ## full FOV exclusions (FOVs should NOT appear in data AT ALL)
    exFOVs <- studyAnn %>%
              select(FOV, dplyr::matches("FOV_exclusion")) %>%
              filter_all(any_vars(str_detect(., pattern = "X"))) %>%
              pull(FOV)
    errors <- warnings <- messages <- list()
    if(any(exFOVs %in% dat$FOV)){
        errors[[length(errors) + 1]] <- paste0("These FOVs are still in the data: ", paste(exFOVs[exFOVs %in% dat$FOV], collapse=", "))
    }

    ## marker exclusions
    mEx <- studyAnn %>%
           select(FOV, Marker_exclusion) %>%
           filter(!is.na(Marker_exclusion))

    for(r in 1:nrow(mEx)){
        ex <- mEx[r,]
        for(m in unlist(strsplit(ex$Marker_exclusion))){
            mDat <- dat %>%
                    filter(FOV == ex$FOV, 
                           grepl(getClassifierPattern(m, delim=","), FullPosMarkerStr))
            if(!is.null(mDat) && nrow(mDat) > 0){
                errors[[length(errors) + 1]] <- paste0("Data still contains marker ",m," in FOV ",ex$FOV)
            }
        }
    }

    return(list(errs=errors, wrns=warnings))

}

#' Get vector of samples that exist in a directory by removing
#' a certain provided extension from all files in the directory
#'
#' List all files in a directory assumed to contain one file per sample
#' by removing an extension from all files
#'
#' @param datDir   directory to check
#' @param ext      extension to remove in order to get sample names
#' 
#' @return vector of file names with extensions removed
samplesInDir <- function(datDir, ext){
    gsub(ext,"",dir(datDir))
}

#' Cross check all input file lists against list of expected samples
#' 
#' All input data should be divided into one file per sample. Here we check
#' the expected list of samples against the the actual files existing in 
#' directories for Halo images and CSV files. Errors/warnings will be logged
#' to flag any inconsistencies
#'
#' @param samplesXLSX   XLSX file containing meta data pertaining to each sample
#' @param imageDir      directory containing images from Halo, where each file name
#'                      contains the name of a single sample (CellDive_ID) in sampleXLSX
#' @param csvDir        directory containing CSV files from Halo, where each file name
#'                      contains the name of a single sample (CellDive_ID) in sampleXLSX
#' 
#' @return expected list of CellDive_IDs
getStudySampleList <- function(samplesXLSX, imageDir, csvDir){
    meta <- read.xlsx(samplesXLSX, 1, check.names=F) %>% as_tibble() %>% pull(CellDive_ID)
    imgs <- samplesInDir(imageDir, "_scanplan_AllFOVs.tif")
    csvs <- samplesInDir(csvDir, "_v1.csv.gz")

    errs <- wrns <- c()

    if(!all(meta %in% imgs)){
        for(m in setdiff(meta,imgs)){ errs <- c(errs, paste0("Sample in meta data but NOT in image dir: ",m)) }
    }
    if(!all(meta %in% csvs)){
        for(m in setdiff(meta,csvs)){ errs <- c(errs, paste0("Sample in meta data but NOT in CSV dir: ",m)) }
    }
    if(!all(imgs %in% meta)){
        for(i in setdiff(imgs,meta)){ wrns <- c(wrns, paste0("Sample in image dir but NOT in meta data: ",i)) }
    }
    if(!all(csvs %in% meta)){
        for(c in setdiff(csvs,meta)){ wrns <- c(wrns, paste0("Sample in CSV dir but NOT in meta data: ",c)) }
    }

    if(length(wrns) > 0){
        for(w in wrns){ log_warn(w) }
    }
    if(length(errs) > 0){
        for(e in errs){ log_error(e) }
    }
    return(meta)
}

#' Validate a single question/comparison
#' 
#' Validate question design by:
#'   (1) checking that all filter keys (column names) match names of columns in meta data
#'   (2) values given for filters match values in meta data
#'   (3) cell region(s) provided are supported (NA|'fov'|'all'|'interface'|'interface inside'|
#'        'interface outside'|'neighborhood')
#'   (4) there is no overlap between the two groups (i.e., the filtering criteria outlined
#'        actually result in each cell belonging to one of the two groups (or neither). 
#' 
#' @param question    table containing two rows, each row describing how to filter the data to
#'                    get one group of cells (see docs for details)
#' @param validMeta   flattened, previously validated table of all study meta data
#' 
#' @return vector of error messages; if empty, question is valid
validateQuestionDesign <- function(question, validMeta){

    errs <- c()

    nonAnnoCompVars <- list("Cell Region" = c(NA, "All", "interface", "interface inside",
                                              "interface outside", "neighborhood"),
                            "cell_MHCI" = c("pos.env","neg.env","mixed.env"),
                            "cell_B7H3" = c("pos.env","neg.env","mixed.env"))

    ## what is to be compared?
    compVar  <- unique(question$ComparisonVariable)
    grpNames <- question[[compVar]]
    grpVals  <- grpNames %>% strsplit(";;") %>% unlist()

    ## check for invalid cell regions
    if(compVar %in% names(nonAnnoCompVars)){
        if(!all(tolower(grpNames) %in% nonAnnoCompVars[[compVar]])){
            err  <- paste0("The following group names are not valid for comparison variable ",
                           compVar, ": ",
                           paste(grpVals[ ! grpVals %in% nonAnnoCompVars], collapse=", "))
            errs <- c(errs, err)
        }
    } else {
        ## check for invalid annotation values
        invalidGrps <- grpVals[ !grpVals %in% c(NA, validMeta[[compVar]]) ]

        if(length(invalidGrps) > 0){
            err  <- paste0("The following group names are not valid for comparison variable ",
                           compVar, ": ",
                           paste(invalidGrps, collapse=", "))
            errs <- c(errs, err)
        }
    }

    ## check that groups are actually different
    if(grpNames[1] == grpNames[2]){
        grpStr <- paste(paste0(c("Group 1: [", "Group 2: ["),
                        paste0(grpNames, "]")), collapse = ", ")
        errs <- paste0("ComparisonVariable seems to be the same in both groups: ",
                        grpStr, ". Is ComparisonVariable wrong?")
    }
    return(errs)

}

#' Validate questions XLSX file
#' 
#' Check that questions are designed and outlined correctly in questions
#' input (XLSX) file. Errors and warnings are logged for user review.
#'
#' @param questionsXLSX  XLSX file describing all questions/comparisons to be
#'                       analyzed
#' @param validMeta      flattened, previously validated meta data
#' @param markers        required only when one or more question includes
#'                       cell-level filter(s) with column named 'cell_[Marker]'
#' 
#' @return logical; when TRUE, questions are valid
validateQuestions <- function(questionsXLSX, validMeta, markers = NULL){

    errs <- list()

    qs <- read.xlsx(questionsXLSX, 1, check.names=F) %>% as_tibble()
    errs <- checkRequiredCols(names(qs), 'questions')

    ## make sure all column labels match those in meta data
    filtCols  <- setdiff(names(qs), req_cols$questions)
    validCols <- c(names(validMeta), "Cell Region")
    if(!is.null(markers)){ validCols <- c(validCols, paste0("cell_", markers)) }

    if(!all(filtCols %in% validCols)){
        invalid <- setdiff(filtCols, validCols)
        msg <- paste0("Filter column(s) in questions file not found in meta data: ",
                      paste(invalid, collapse = ", "))
        errs[[length(errs) + 1]] <- msg
    }

    ## make sure all values in filter columns match those in meta data
    cols2check <- intersect(filtCols, names(validMeta))
    for(fc in cols2check){
        obsVals <- unique(unlist(strsplit(as.character(qs[[fc]]), ";;")))
        expVals <- c(unique(validMeta[[fc]]), NA)
        if(!all(obsVals %in% expVals)){
            errs[[length(errs) + 1]] <- paste0("Invalid value(s) in column ",fc," of questions file: ",
                                               paste(obsVals[!obsVals %in% expVals], collapse=", "))
        }
    }

    ## check each question
    for(qu in unique(qs$QuestionNumber)){
        log_debug(paste("Validating question:", qu))
        quest <- qs %>% filter(QuestionNumber == qu)
        qErrors <- validateQuestionDesign(quest, validMeta)
        if(!is.null(qErrors)){
            errs[[length(errs) + 1]] <- paste0(qu, ": ",qErrors)
        }
    }

    if(length(errs) > 0){
        for(e in unlist(errs)){ log_error(e) }
        return(FALSE)
    }
    return(TRUE)

}

#' Check for errors in all input data including questions, study annotation and
#' cell type definitions
#' 
#' Validate questions file, study annotation and cell type definitions. Any warnings
#' and errors will be logged. If errors are found, exit with status 1.
#'
#' @param study_config    study configuration in list format, with at minimum 
#'                        'meta_dir' defined
#'
#' @return nothing; quit with exit status 1 if any errors found
validateInput <- function(study_config){

    log_info("Validating meta files...\n")
    mFiles <- getFiles(path = study_config$meta_dir, pattern="*.xlsx")
    saValid <- validateStudyAnnotation(mFiles[grepl("FOVs.xlsx",mFiles)],
                                       mFiles[grepl("Samples.xlsx",mFiles)],
                                       mFiles[grepl("Patients.xlsx",mFiles)])
    qValid <- FALSE
    if(saValid){
        log_info("Validating questions file...\n")
        meta <- loadStudyAnnotations(metaFiles = mFiles)$flat
        markers <- read.xlsx(mFiles[grepl("Markers.xlsx",mFiles)], 1, check.names = F) %>%
                   as_tibble() %>%
                   getAllMarkers()
        qValid <- validateQuestions(study_config$statistics_questions_file, meta, markers = markers)
    }

    log_info("Validating markers and cell types...\n")
    ctValid <- validateCellTypes(mFiles[grepl("CellTypes",mFiles)], mFiles[grepl("Markers",mFiles)])

    if(!all(saValid, qValid, ctValid)){
        log_fatal("Validation FAILED.")
        cat("\n")
        q(save="no", status=1)
    }

}

