#' Convert string representing a marker combination into a 
#' string of each marker's corresponding 'rank' in the list of all markers
#'
#' Each marker is assigned a double digit number and the digits are collapsed
#' into a non-delimited string (e.g., 010509) to allow for correct sorting 
#'
#' @param combo   comma-delimited string of marker combination
#' @param mrkrs   vector of all individual markers ordered by priority
#'
#' @return  character string of formatted integers to allow for proper sorting
customComboOrder <- function(combo, mrkrs){
   if(is.na(combo)){ return(NA) }
   y <- unlist(strsplit(combo,","))
   priorities <- which(mrkrs %in% y)
   ordr <- formatC(priorities, flag = "0", width = 2)
   paste(ordr, collapse="")
}

#' Arrange table of conditions based on specific 'instructions' in list form
#' 
#' Arrange table of conditions according to the format of list 'howToArrange'
#'
#' @param conds          unsorted table of conditions
#' @param howToArrange   list describing how to arrange table:
#'                         - list is named by columns on which table should be sorted
#'                         - the order of list names is the order in which data
#'                           will be arranged
#'                         - if list element is not null, it should contain a list or vector
#'                           of column values in the order in which they should occur in 
#'                           the sorted data table
#'                         - if the list element is null/empty, values in that column will 
#'                           be naturally sorted (numerically for numbers and alphanumerically 
#'                           for character strings
#'
#' @return sorted conditions table
arrangeConditions <- function(conds, howToArrange){
    arrange_by <- c()
    col_order <- names(howToArrange)

    for(col in col_order[col_order %in% names(conds)]){
        if(!is.null(howToArrange[[col]])){
            cOrder <- howToArrange[[col]] %>% unlist(use.names=F)
            conds[[paste0("order_",col)]] <- customOrder(conds[[col]], cOrder)
            arrange_by <- c(arrange_by, paste0("order_",col))
        } else {
            arrange_by <- c(arrange_by, col)
        }
    }
    conds <- conds %>%
             arrange_at(arrange_by) %>%
             select(-dplyr::matches("order_")) %>%
             unique()
    conds
}

#' Parse condition string into its individual parts
#'
#' Parse condition string to determine the Category, Cell_type, Subtype,
#' Proliferation (KI67), Functional marker(s) and combination(s)
#' 
#' @param condStr       comma-delimited summary string of the cell type condition of interest
#' @param cellTypes     expanded table describing all cell types of interest (result from getCellTypes())
#' @param funcMarkers   vector of all functional markers
#' @param funcCombos    list of character vectors, each of which contains functional markers that
#'                      make certain combinations of interest; 
#'                         * NOTE: these vectors should include negatives 
#'                                 (e.g., c(PD1, PD1-, TIM3, TIM3-, LAG3, LAG3-))
#'                         * NOTE 2: markers will be sorted in the order they appear in the vector
#'
#' @return  single-row tibble with columns mentioned above filled in
normalizeCondition <- function(condStr, cellTypes, funcMarkers, funcCombos){
    tibble(CellDesc      = condStr,
           Category      = getCellTypeClass(condStr, cellTypes, "Category"),
           Cell_type     = getCellTypeClass(condStr, cellTypes, "Cell_type"),
           Subtype       = getCellTypeClass(condStr, cellTypes, "Subtype"),
           Functional    = extractMarkerCombos(condStr,
                                            list(unique(c(funcMarkers, paste0(funcMarkers,"-")))),
                                            singlesOnly = TRUE),
           Proliferation = ifelse(grepl(getClassifierPattern("KI67",delim=","), condStr), "KI67", NA),
           `Functional Combination` = extractMarkerCombos(condStr, funcCombos, excludeSingles=TRUE))
}

#' Parse, sort and assign unique IDs to non-spatial conditions
#' 
#' Parse, sort and assign unique IDs to non-spatial conditions
#' 
#' @param condFile           XLSX file containing a separate tab for each calculation type,
#'                           where each table contains a full list of conditions to be
#'                           analyzed by that calculation bype
#' @param sheets             names of sheets considered 'general' analyses
#' @param arrangeAnnotation  list describing how conditions should be arranged (see arrangeConditions())
#' @param firstID            the first ID to use when assigning cell state IDs (default = 1)
#'
#' @return table of parsed & sorted conditions with unique IDs assigned
indexGeneralConditions <- function(condFile, sheets, arrangeAnnotation = NULL, firstID = 1){

    allConds <- tibble()
    for(sheet in sheets){
        conds <- tryCatch({
                     xlsx::read.xlsx(condFile, sheetName = sheet, check.names=F) %>%
                     as_tibble %>%
                     mutate(CalcType = sheet, AnalysisType = "general")
                 }, error = function(e){
                     log_warn(paste0("No sheet [", sheet, "] found. skipping."))
                 })
        if(nrow(allConds) == 0){
            allConds <- conds
        } else {
            joinCols <- intersect(names(allConds), names(conds))
            joinCols <- joinCols[!joinCols == "CalcType"]
            allConds <- allConds %>%
                        full_join(conds, by = joinCols)
        }
    }

    ## sort conditions
    if(!is.null(arrangeAnnotation)){
        allConds <- allConds %>% arrangeConditions(arrangeAnnotation)
    }

    ## QC: report how many conditions exist for all calc types and how many
    ## for each type exclusively 
    allConds <- allConds %>%
                unite("CalcType", dplyr::matches("CalcType"), sep = ",", na.rm = T)

    nums <- allConds %>% group_by(CalcType) %>% summarize(Count = n())
    for(nm in 1:nrow(nums)){
        log_debug(paste0("CalcType(s): ",nums$CalcType[nm],"\t", nums$Count[nm], " conditions"))
    }

    allConds$`Cell State ID` <- as.character(as.numeric(rownames(allConds)) + (firstID - 1))
    allConds
}

#' Parse, sort and assign unique IDs to spatial conditions
#' 
#' Parse, sort and assign unique IDs to spatial conditions
#' 
#' @param condFile           XLSX file containing a separate tab for each calculation type,
#'                           where each table contains a full list of conditions to be
#'                           analyzed by that calculation bype
#' @param sheets             names of sheets considered 'spatial' analyses
#' @param arrangeAnnotation  list describing how conditions should be arranged (see arrangeConditions())
#' @param firstID            the first ID to use when assigning cell state IDs (default = 1)
#'
#' @return table of parsed & sorted conditions with unique IDs assigned
indexSpatialConditions <- function(condFile, sheets, arrangeAnnotation = NULL, firstID = 1){

    allConds <- tibble()
    for(sheet in sheets){
        conds <- tryCatch({
                     xlsx::read.xlsx(condFile, sheetName = sheet, check.names=F) %>%
                     as_tibble %>%
                     mutate(CalcType = sheet, AnalysisType = "spatial")
                 }, error = function(e){
                     log_warn(paste0("No sheet [", sheet, "] found. skipping."))
                 })
        if(is.null(conds)){ next }
        if("Center" %in% names(conds)){
            conds <- conds %>% rename(`Center Population A` = Center)
        }
        if("Neighborhood" %in% names(conds)){
            conds <- conds %>% rename(`Neighborhood Population A` = Neighborhood)
        }
        if(nrow(allConds) == 0){
            allConds <- conds
        } else {
            joinCols <- intersect(names(allConds), names(conds))
            joinCols <- joinCols[!joinCols == "CalcType"]
            allConds <- allConds %>%
                        full_join(conds, by = joinCols)
        }
    }

    ## sort conditions
    if(!is.null(arrangeAnnotation)){
        allConds <- allConds %>% arrangeConditions(arrangeAnnotation)
    }

    ## QC: report how many conditions exist for all calc types and how many
    ## for each type exclusively
    allConds <- allConds %>% unite("CalcType", dplyr::matches("CalcType"), sep = ",", na.rm = T)

    nums <- allConds %>% group_by(CalcType) %>% summarize(Count = n())
    for(nm in 1:nrow(nums)){
        log_debug(paste0("CalcType(s): ",nums$CalcType[nm],"\t", nums$Count[nm], " conditions"))
    }

    allConds$`Cell State ID` <- as.character(as.numeric(rownames(allConds)) + (firstID - 1))
    allConds
}



#' Compile, expand, and assign unique IDs to all conditions
#' 
#' Create table of fully parsed conditions from all lists provided and
#' assign unique IDs to each condition
#' 
#' @param condFile           XLSX file containing a tab for each calculation type, where
#'                           each tab contains all conditions for which that calculation should be made
#' @param arrangeAnnotation  list describing how conditions should be arranged
#' 
#' @return table of parsed & sorted conditions with unique IDs assigned
indexConditions <- function(condFile, arrangeAnnotation){

    generalSheets <- c("Fraction","Density")
    spatialSheets <- c("Neighborhood_averagecounts", "Neighborhood_fraction", "Pos_Neg_Microenvironments")
    general <- indexGeneralConditions(condFile, generalSheets, arrangeAnnotation)

    lastGenID <- max(as.numeric(general$`Cell State ID`))
    spatial <- indexSpatialConditions(condFile, spatialSheets, arrangeAnnotation, firstID = lastGenID + 1) %>%
               mutate(Condition = paste(`Center Population A`, `Neighborhood Population A`, sep="____"),
                      Population = paste(`Center Population B`, `Neighborhood Population B`, sep="____"))

    col_order <- names(arrangeAnnotation)
    allConds <- bind_rows(general, spatial) %>%
                select_at(c("Cell State ID", "Condition", "Population",
                          col_order[col_order %in% c(names(general), names(spatial))], "AnalysisType", "CalcType"))

    allConds
}


#' Convert a string containing underscores to separate words, all capitalized
#' 
#' Convert a string containing underscores to separate words, all capitalized
#' 
#' @param str    character string containing underscores (e.g., 'Column_title')
#'
#' @return capitalized string (e.g., 'Column Title')
underscoreToUpper <- function(str){
    mtch <- str_match_all(str, "_.{1}")[[1]][,1]
    for(m in mtch){
        str <- gsub(m, toupper(gsub("_"," ",m)), str)
    }
    str
}

#' Map all sample and FOV IDs from annotation to halo data
#'
#' Join all IDs from annotation data to halo data so all can be used
#' interchangeably downstream
#'
#' @param dat    table of halo data, including at minimum columns Sample and SPOT,
#'               where 'Sample' corresponds to 'CellDive_ID' in annotation and 
#'               'SPOT' corresponds to 'FOV_number' in annotation
#' @param idMap  table of all IDs including at minimum CellDive_ID and FOV_number
#'
#' @return data table with any missing IDs added
joinIDs <- function(dat, idMap){
    if(all(c("Sample", "SPOT") %in% names(dat))){
        dat <- dat %>%
               select(CellDive_ID = Sample,
                      FOV_number = SPOT,
                      everything())
    } else if(all(c("CellDive_ID", "FOV_number") %in% names(dat))){
        dat <- dat %>% select(CellDive_ID, FOV_number, everything())
    } else {
        msg <- "Could not join IDs."
        log_error(msg)
        stop(msg)
    }
    dat %>%
    left_join(idMap, by=intersect(names(dat), names(idMap)))
}

#' Determine whether all values in a vector are NA
#' 
#' Determine whether there are any actual values (non-NA) in a vector
#'
#' @param x   vector of any value type
#' 
#' @return  logical; TRUE if at least one value in x is NOT NA
not_all_na <- function(x) {!all(is.na(x))}

#' Extract cell classifier from a comma-delimited character
#' string representing a specific cell state
#' 
#' Pull out just the cell type/subtype/tag from a full comma-delimited
#' character string representing a specific state. Specifically, this function
#' returns only the first element in a comma-delimited string (everything
#' to the left of the first comma). If the string contains zero commas,
#' the entire string is returned.
#'
#' @param x   tab-delimited cell state character string
#'
#' @return the first element of x (everything to the left of the first comma) 
extractCT <- function(x) { gsub(",.*", "", x) }


#' Get classifier(s) of a certain type based on classifier(s) of a different type 
#' 
#' Given a vector of classifiers, get the corresponding classifiers of a single type.
#' For example, for a vector with classifiers (Tconv4, MHCIIpos_macro, T cell), return
#' the corresponding Subtype for each of them ("CD4+ T cell", "MHCII+ macrophage", "All").
#' Note, because there are multiple subtypes that fall under classifier "T cell", the string
#' "All" is returned.
#'
#' @param ct         vector of classifiers, of any type
#' @param cellTypes  table map of all classifiers including columns "Category", 
#'                   "Cell_type", "Subtype", "Tag"
#' @param facetCol   Column in cell types containing the classifier type to be returned
#'
#' @return vector of classifiers of type [facetCol] corresponding to each one in [ct]
getClassifier <- function(ct, cellTypes, facetCol){
    ctLvls <- c("Category", "Cell_type", "Subtype", "Tag")
    celltypes <- cellTypes %>% select_at(ctLvls)
    sapply(ct, function(x){
        catType <- celltypes %>%
                   filter_at(vars(ctLvls), any_vars(grepl(paste0("^",x,"$"), .))) %>%
                   pull(facetCol) %>% unique()
        ifelse(is.null(catType) || is.na(catType) || length(catType) > 1, 
               ifelse(!any(catType == x), "All", x), catType) 
    }) %>%
    unlist()
} 


#' Add '+' to positive markers in a comma-delimited cell state string
#' 
#' Given a comma-delimited character string representing a cell state, add
#' a '+' to elements determined to be positive markers. Any element in 
#' the split string that is NOT the first and does NOT contain a '-' is
#' assumed to be a positive marker.
#'
#' @param condStr   comma-delimited string representing a cell state; the
#'                  first element is a cell classifier and the remaining elements
#'                  are any combination of positive and/or negative markers
#' 
#' @return cell state string with positive markers explicitly labeled with a '+'
formatMarkers <- function(condStr){
    sapply(condStr, function(x){
        if(!grepl(",", x)){ return(x) }
        spl <- unlist(strsplit(x, ","))
        posMrkrs <- which(!grepl("\\-$", spl))[-1]
        spl[posMrkrs] <- paste0(spl[posMrkrs], "+")
        paste0(spl, collapse = ",")
    }) %>% unlist()
}



getAbbrev <- function(longAbbrev, subscript){
    sapply(1:length(longAbbrev), function(x){
        gsub(subscript[x], "", longAbbrev[x])
    }) %>% unlist()
}


fillInMissingClassifiers <- function(dat){

    ordr <- c("Category", "Cell_type", "Subtype", "Abbrev_label_for_figures")
    ordr <- ordr[ordr %in% names(dat)]
    
    for(x in 1:(length(ordr) - 1)){
        col <- ordr[x]
        colWblanks <- ordr[x + 1]
        idxs <- which(is.na(dat[[colWblanks]]))
        if(!is.na(idxs) && length(idxs) > 0){
            dat[[colWblanks]][idxs] <- dat[[col]][idxs]
        }
    }
    dat
}

removeRedundantClass <- function(state, pop){
    if(length(pop) != length(state)){
        log_warn("lengths differ between vectors or cell states and populations. there may be errors in the removal of population strings from state strings.")
    }
    sapply(1:length(state), function(x){
        gsub(paste0(pop[x],","), "", state[x])
    }) %>%
    unlist()
}



getAllStatsHeatmapData <- function(questionsFile, resDir, calcType,
                              qNums = "all", resSheet="fractions_complete", resFileSfx="_all.xlsx",
                              filtSheet = NULL, maxFDR = 0.05, cellStates = NULL){

    if(length(qNums) == 1 && tolower(qNums) == "all"){
        qNums <- openxlsx::read.xlsx(questionsFile) %>% as_tibble() %>% pull(QuestionNumber)
    }

    fin_tbl <- tibble()

    for(qNum in qNums){
        if(length(trimws(qNum)) == 0){ next }
        statsFile <- file.path(resDir, paste0(qNum,resFileSfx))
        cat(paste0("Reading stats from file ",statsFile,"\n"))
        if(!file.exists(statsFile)){ print(paste0("WARNING: File does not exist for question ",qNum,".")); next }

        ## get filtered data, if any
        qFilt <- NULL
        if(!is.null(filtSheet)){
            qFilt <- xlsx::read.xlsx(statsFile, sheetName = filtSheet, check.names=F) %>% as_tibble()
        }

        ## get complete list of data
        qDat <- NULL
        qDat <- xlsx::read.xlsx(statsFile, sheetName = resSheet, check.names=F) %>% as_tibble()
        if(is.null(qDat)){ stop(paste0("Could not read stats from file ",statsFile)) }

        if(calcType == "fractions"){
            qDat <- qDat %>%
                    mutate(Question  = qNum,
                           AnalysisCondition = paste0(`Cell State`, ";;", Population),
                           signif    = ifelse(`Fraction adjusted p.value` < maxFDR, "*", NA)) %>%
                    select(`Cell State ID`, `Cell State`, AnalysisCondition, dplyr::matches("median fraction"),
                           `log(Odds Ratio)`, `Fraction adjusted p.value`, signif)
            if(!is.null(qFilt)){
                qFilt <- qFilt %>% mutate(AnalysisCondition = paste0(`Cell State`, ";;", Population))
            }
        } else {
            qDat <- qDat %>%
                    mutate(Question  = qNum,
                           AnalysisCondition = `Cell State`,
                           signif    = ifelse(`Density adjusted p.value` < maxFDR, "*", NA)) %>%
                    select(`Cell State ID`, `Cell State`, AnalysisCondition, dplyr::matches("median density"),
                           `log(Fold Change)`, `Density adjusted p.value`, signif)
            if(!is.null(qFilt)){
                qFilt <- qFilt %>% mutate(AnalysisCondition = `Cell State`)
            }
        }
        if(!is.null(qFilt)){
            qDat <- qDat %>% mutate(PASS = ifelse(AnalysisCondition %in% qFilt$AnalysisCondition, 1, 0))
        }

        qDat$QNum <- qNum
        qTitle <- gsub("\\.", " ",
                      paste(gsub("median fraction", "", names(qDat)[grep("median fraction",names(qDat))]),
                            collapse="\nvs\n"))
        qDat[grepl("median", names(qDat)) ] <- NULL


        ##### special, temp mods to question title
        qTitle <- gsub("Non-Responder__Partial Responder", "Non/Partial\nResponder", qTitle)
        qTitle <- gsub("Complete Responder","Pretreatment\nComplete\nResponder", qTitle)
        qTitle <- gsub("^\\+|\\+$","IL2 Treated", trimws(qTitle))
        qTitle <- gsub("^\\-|\\-$", "Untreated", trimws(qTitle))
        qDat$QTitle <- qTitle
        ######

        if(!is.null(cellStates)){
            qDat <- qDat %>% left_join(cellStates, by=intersect(names(qDat), names(cellStates)))
        }

        if(nrow(fin_tbl) == 0){
            fin_tbl <- qDat
        } else {
            fin_tbl <- fin_tbl %>% bind_rows(qDat)
        }
    }

    if(!is.null(filtSheet)){
        fin_tbl <- fin_tbl %>% group_by(`Cell State ID`, `Cell State`) %>% mutate(`Number Significant` = sum(PASS))
    }

    return(fin_tbl)

}

getCellTypeClass <- function(popStr, cellTypes, ctCol){
   classOrder <- c("Category", "Cell_type", "Subtype", "Tag")
   ct <- unlist(strsplit(popStr,","))[1]
   rows <- cellTypes %>% filter_at(vars(classOrder), any_vars(. == ct))
   if(nrow(rows) == 1 || length(unique(rows[[ctCol]])) == 1){
       rows %>% pull(ctCol) %>% unique()
   } else {
       NA
   }
}

numPositives <- function(markerCombo){
    if(is.na(markerCombo)){ return(NA) }
    t <- unlist(strsplit(markerCombo,","))
    length(t[!grepl("\\-",t)])
}

toMatrix <- function(df, rowOrder=NULL, colName="Condition"){
    rNames <- df[,colName]
    df <- as.matrix(df[,-which(names(df) == colName)])
    rownames(df) <- rNames
    if(!is.null(rowOrder)){ return(as.matrix(df[rowOrder,])) }
    as.matrix(df)
}

extractMarkerCombos <- function(popStr, comboList, excludeSingles=FALSE, singlesOnly=FALSE,
                                sortCombo=FALSE, dropNegatives=FALSE){
    pop <- unlist(strsplit(popStr,","))
    combo <- NA
    for(x in seq(comboList)){
        xMarkers <- comboList[[x]]
        found <- pop[pop %in% xMarkers]
        combo <- NA
        num <- length(found)

        if(singlesOnly && excludeSingles){
            stop("excludeSingles and singlesOnly can not both be TRUE")
        }
        if((singlesOnly && num == 1) || (excludeSingles && num > 1) || (!excludeSingles & !singlesOnly)){
            combo <- found
            combo <- sapply(combo, function(x){ which(xMarkers == x) }) %>%
                     unlist() %>%
                     sort() %>%
                     names()
        }
        if(!is.na(combo)){ break } ## assume there is only one possible option
    }
    if(!is.na(combo)){
        if(dropNegatives){
            combo <- combo[!grepl("\\-",combo)]
            if(combo == ""){ combo <- "Negative" }
        }
        if(sortCombo){ combo <- sort(combo) }
        combo <- paste(combo, collapse = ",")
    }
    combo

}



spreadIndivDat <- function(indivDat, comparison, grp1, grp2, calc, effect, idMap, currentIDcol, newIDcol){
    tbl <- indivDat %>% 
           select(`Cell State ID`, `Cell State`, Population,
                  !!as.name(paste(grp1, "Overall Median", calc)) := OverallGroup1Median,
                  !!as.name(paste(grp2, "Overall Median", calc)) := OverallGroup2Median,
                  !!as.name(paste(comparison, "Overall", effect)) := OverallEffect,
                  !!as.name(paste(comparison, "Overall CI.low")) := !!as.name(paste(calc, "CI.low")),
                  !!as.name(paste(comparison, "Overall CI.high")) := !!as.name(paste(calc, "CI.high")),
                  !!as.name(paste(comparison, "Overall p.value")) := OverallPval,
                  !!as.name(paste(comparison, "Overall FDR")) := OverallFDR,
                  !!as.name(paste(comparison, "HMP")) := HMP,
                  !!as.name(paste(comparison, "adjusted HMP")) := FDR) %>%
            unique()

    for(pl in unique(indivDat[[currentIDcol]])){
        grp1 <- idMap %>% filter(!!as.name(currentIDcol) == pl, Group == "Group 1") %>% pull(newIDcol)
        grp2  <- idMap %>% filter(!!as.name(currentIDcol) == pl, Group == "Group 2") %>% pull(newIDcol)
        iComp <- paste0(grp1, " vs ", grp2)
        tmp <- indivDat %>%
               filter(!!as.name(currentIDcol) == pl) %>%
               select(`Cell State ID`, `Cell State`, Population,
                      !!as.name(paste(iComp, effect)) := IndivEffect,
                      !!as.name(paste(iComp, "p.value")) := IndivPval,
                      !!as.name(paste(iComp, "p.value one side")) := `p.val one side`,
                      !!as.name(paste(iComp, "FDR")) := IndivFDR,
                      !!as.name(paste(iComp, "vs Overall")) := Status)
        tbl <- tbl %>% left_join(tmp, by = intersect(names(tmp), names(.)))
    }
    tbl
}

getIDmap <- function(allQuestions, qPre = NULL, qNames = NULL){
    if(is.null(qNames)){
        if(is.null(qPre)){ 
            msg <- "Please provide either 'qPre' or 'qNames' to get ID map"
            log_error(msg)
            stop(msg)
        }
        qNames     <- names(allQuestions)[grepl(paste(qPre, collapse = "|"), names(allQuestions))]
    }

    idMap <- lapply(allQuestions[qNames], function(x){
                    qu <- x$qNum
                    qTbl <- tibble()
                    for(grp in c("Group 1", "Group 2")){
                        qTbl <- qTbl %>%
                        bind_rows(lapply(x[[grp]], function(y) paste(y, collapse = ",") ) %>%
                                         as_tibble() %>%
                                         mutate(Question = qu, Group = grp))
                    }        
                    qTbl
                }) %>%
             bind_rows() %>%
             select(Question, Group, everything())

    idMap
}


#' Format neighborhood counts table for statistical analyses
#' 
#' Filter table for cells in region of interest and join center cell annotation
#' including Sample_ID, Band and PositiveMarkers. Fill in table zeros where 
#' applicable, as raw counts tables include only center cell/neighborhood cell
#' type pairs that exist (i.e., N.Count >= 1)
#' 
#' @param nbhdCountsDir    directory of files each containing table of neighborhood counts 
#'                         including columns 'C.UUID' and 'FOV'
#' @param annCells         annotated cell data for ONLY cells inside region of interest
#' @param cellDiveID       CellDive_ID for which data should be formatted; default: All
#' 
#' @return  filtered & annotated neighborhood counts table
formatNeighborhoodCounts <- function(nbhdCountsDir, annCells, cellDiveID = "All"){

    ncFiles <- file.path(nbhdCountsDir, dir(nbhdCountsDir))
    if(tolower(cellDiveID) != "all"){ 
        ncFiles <- ncFiles[grepl(paste0("^", cellDiveID, "_"), basename(ncFiles))] 
        checkFilesFound(ncFiles, nbhdCountsDir, "neighborhood counts")
    }

    if("Band" %in% names(annCells)){
        dat <- annCells %>% select(Band, C.UUID = UUID, C.PositiveMarkers = PositiveMarkers)
    } else {
        dat <- annCells %>% select(C.UUID = UUID, C.PositiveMarkers = PositiveMarkers)
    }

    allNbhdCounts <- tibble()

    for(ncFile in ncFiles){
        log_debug(paste0("Formatting nbhd counts in file: ",ncFile))
        nCts <- readRDS(ncFile)
        if(!"FOV_ID" %in% names(nCts)){
            nCts <- nCts %>% rename(FOV_ID = FOV)
        }

        FOV_IDs <- unique(nCts$FOV_ID)

        ## fill in missing neighborhood cell types
        nCts <- nCts %>%
                spread(NeighborhoodCellType, N.Count, fill=0) %>%
                gather(unique(nCts$NeighborhoodCellType), 
                       key="NeighborhoodCellType", 
                       value="N.Count")

        tmp <- nCts %>%
               select(FOV_ID, CenterCellType, NeighborhoodCellType) %>%
               unique() %>%
               mutate(N.Count = 1) %>%
               spread(FOV_ID, N.Count, fill=0)

        ## make sure all center cell type/neighborhood cell type pairs
        ## are accounted for in all FOVs
        if(!all(FOV_IDs %in% names(tmp))){
            for(f in FOV_IDs[!FOV_IDs %in% names(tmp)]){
                tmp[[f]] <- 0
            }
        }

        ## remove redundant rows (those with N.Count != 0)
        tmp <- tmp %>%
               gather(FOV_IDs, key="FOV_ID", value="N.Count") %>%
               filter(N.Count == 0)

        allNbhdCounts <- allNbhdCounts %>%
                         bind_rows(nCts %>%
                                   bind_rows(tmp) %>%
                                   left_join(dat, by=c("C.UUID")))

    }

    allNbhdCounts
}

#' Get vector of integers representing the order in which a character vector 
#' should be sorted 
#'
#' Given an unordered character vector and a vector of the same strings ordered
#' in a specific (non alphanumeric) way, create a vector of integers representing the 
#' index at which each item should be when sorted correctly. 
#'  
#' @param vals    vector of unordered character strings
#' @param cOrder  vector of the same character strings, sorted in the custom order
#' @return vector of integers reprsenting the order in which vals should be sorted
#'         to get the custom order
customOrder <- function(vals, cOrder){
    ordr <- sapply(sort(cOrder), function(x){ which(cOrder == x) })
    sapply(vals, function(x){ ordr[x] }) %>% unlist(use.names=F)
}

#' Push values 'All' and NA to end of a vector
#'
#' Push values 'All' and NA to end of a vector; this is generally used
#' when ordering lists of cell types for figures
#' 
#' @param vals   vector of character strings
all_and_na_last <- function(vals){
    c(vals[!vals %in% c("All", NA)], "All", NA)
}

