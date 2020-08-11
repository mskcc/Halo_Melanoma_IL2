#' Convert string representing a marker combination into a 
#' string of each marker's corresponding 'rank' in the list of all markers
#'
#' Each marker is assigned a double digit number and the digits are collapsed
#' into a non-delimited string (e.g., 010509) to allow for correct sorting 
#'
#' @param combo   comma-delimited string of marker combination
#' @param mrkrs   vector of all individual markers ordered by priority
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

not_all_na <- function(x) {!all(is.na(x))}

extractCT <- function(x) { gsub(",.*", "", x) }

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

formatMarkers <- function(condStr){
    sapply(condStr, function(x){
        if(!grepl(",", x)){ return(x) }
        spl <- unlist(strsplit(x, ","))
        posMrkrs <- which(!grepl("\\-$", spl))[-1]
        spl[posMrkrs] <- paste0(spl[posMrkrs], "+")
        paste0(spl, collapse = ",")
    }) %>% unlist()
}

getFractionConditionLabel <- function(numerator, denominator, delim = "|"){
    if(!length(numerator) == length(denominator)){
        log_warn("length differ between input vectors. results may be incorrect.")
    }
    sapply(1:length(numerator), function(x){
        #cnd <- gsub(paste0(denominator[x], ","), "", numerator[x]) 
        cnd <- numerator[x]
        pop <- denominator[x] 
        paste(cnd, delim, pop)
    }) %>% unlist()
}

getDensityConditionLabel <- function(value){
    sapply(1:length(value), function(x){
        ifelse(grepl(",", value[x]), 
            paste0(unlist(strsplit(value[x], ","))[-1], collapse = ","),
            "All")
    }) %>% unlist()
}

getCountConditionLabel <- function(value){
    value
}

assignClasses <- function(dat, facets, cellTypes, facetOrder){ 
    grps <- dat %>% 
            select_at(unique(c(facets, "Condition", "Population"))) %>% 
            mutate(ct = extractCT(Condition))
    for(classType in facets[facets %in% names(cellTypes)]){
        grps <- grps %>%
                mutate(!!as.name(classType) := getClassifier(ct, cellTypes, classType))
        grps[[classType]][is.na(grps[[classType]])] <- ""
    }
    grps <- grps %>% select(-ct) %>% unique()

    for(classType in facets){
        grps[[paste0(classType, "_order")]] <- customOrder(grps[[classType]], facetOrder[[classType]])
        grps[[classType]] <- factor(grps[[classType]], levels = facetOrder[[classType]])
    }
    dat %>%  
    select_at(names(dat)[!names(dat) %in% facets]) %>%
    left_join(grps, by=intersect(names(.), names(grps))) %>%
    arrange_at(paste0(facets, "_order")) 
}

getCondOrderByPopulation <- function(dat, popOrder = NULL){
    if(!is.null(popOrder)){
        vals <- extractCT(dat$Population)
        dat$popOrder <- customOrder(vals, popOrder)
    } else {
        dat$popOrder <- as.numeric(dat$`Cell State ID`)
    }
    dat
}


getAbbrev <- function(longAbbrev, subscript){
    sapply(1:length(longAbbrev), function(x){
        gsub(subscript[x], "", longAbbrev[x])
    }) %>% unlist()
}

setConditionOrder <- function(conds, ids, cellTypes, statsFile = NULL, sheet = NULL, 
                              orderBy = "Cell State ID", calcType = "fractions",
                              facetY = NULL, facetOrder = NULL, popOrder = NULL, idOrder = NULL){

    ids <- as.numeric(ids)

    ## eliminate any irrelevant columns
    cnds <- conds %>%
            mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
            filter(`Cell State ID` %in% ids) %>%
            select(`Cell State ID`, Condition, Population, Cell_type, Subtype) %>%
            mutate(Condition = formatMarkers(Condition),
                   Population = formatMarkers(Population)) %>%
            left_join(cellTypes %>% select(Cell_type, Subtype, Tag, Abbrev, Subscript),
                      by = intersect(names(.), c("Cell_type", "Subtype", "Tag")))

    ## data already contains Cell_type column associated with NUMERATOR condition, 
    ## but we need to change that in this case to match the DENOMINATOR condition
    if(!is.null(facetY)){
        cnds <- assignClasses(cnds, facetY, cellTypes, facetOrder = facetOrder)
    } else {
        cnds$facet <- NA
    }

    if(calcType == "fractions"){
        cnds <- cnds %>%
                mutate(Cond = getFractionConditionLabel(Condition, Population))
    } else if(calcType == "densities"){
        cnds <- cnds %>%
                mutate(Cond = getDensityConditionLabel(Condition))
    }

    if(!is.null(idOrder)){
        cnds <- cnds %>%     
                mutate(ID_order = customOrder(as.character(`Cell State ID`), as.character(idOrder))) %>%
                arrange(desc(ID_order)) %>%
                mutate(labelY = row_number()) %>%
                group_by_at(facetY) %>%
                mutate(y = row_number()) %>%
                ungroup() %>%
                select_at(c("Cell State ID", "Condition", "Population", "Cond", facetY, "y",
                         "Tag", "Abbrev","Subscript", "labelY"))
        cnds$Tag[cnds$Tag == ""] <- "All"
        cnds$y <- factor(cnds$y, levels = unique(cnds$y))
        return(cnds)
    }


    ## if rows are to be ordered based on denominator population, we need a new column
    ## on which we can sort
    cnds <- cnds %>% getCondOrderByPopulation(popOrder = popOrder)

    if(!orderBy %in% names(cnds)){
        if(any(is.null(c(statsFile, sheet)))){ 
            stop(paste0("Can not order by ", orderBy, ". Did you forget to provide statsFile or sheet?")) 
        }
        dat <- xlsx::read.xlsx(statsFile, sheet, check.names = F) %>% as_tibble() 
        cnds <- cnds %>% 
                left_join(dat %>% 
                          mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>% 
                          select_at(c("Cell State ID", orderBy)),
                          by = intersect(names(.), c("Cell State ID", orderBy)))
    }

    ##### TO DO: FIX THIS TO NOT IMPROPERLY ASSUME THAT 1) FACETY IS NOT NULL AND 2) LENGTH FACETY IS 2
    cnds <- cnds %>%
            mutate(!!as.name(paste0(orderBy,"_order")) := ifelse(!!as.name(orderBy) <= 1, 1, 0)) %>%
            arrange(desc(!!as.name(paste0(facetY[1],"_order"))), desc(!!as.name(paste0(facetY[2],"_order"))),
                    desc(!!as.name(paste0(orderBy, "_order"))), abs(log(!!as.name(orderBy)) )) %>%
            mutate(labelY = row_number()) %>%
            group_by_at(facetY) %>%
            mutate(y = row_number()) %>%
            ungroup() %>%
            select_at(c("Cell State ID", "Condition", "Population", "Cond", facetY, "y",
                        "Tag", "Abbrev","Subscript", "labelY")) 

    cnds$Tag[cnds$Tag == ""] <- "All"
    cnds$y <- factor(cnds$y, levels = unique(cnds$y))
    cnds

}

formatConditionsForPlotting <- function(conds, includeIDs, calcType, cellTypes, orderBy = NULL,
                                        statsFile = NULL, sheet = NULL, idOrder = NULL, 
                                        facetY = NULL, facetOrder = NULL, popOrder = NULL){

    cnds <- setConditionOrder(conds, includeIDs, cellTypes, orderBy = orderBy,
                              statsFile = statsFile, sheet = sheet, idOrder = idOrder,
                              facetY = facetY, facetOrder = facetOrder)

    cnds <- cnds %>%
            gather(c("Cond", "Cell State ID"), key = "Column", value = "Value") %>%
            mutate(hjust = ifelse(Column == "Cond", 0.5, 1))

    if(ct == "densities"){
        cnds$Value = gsub(" \\|.*", "", cnds$Value)
    }

    cnds$Column <- factor(cnds$Column, levels = c("Cond", "Cell State ID"))
    cnds
}


formatForOverallvsIndiv <- function(overallFile, indivFiles, idsToPlot, max.fdr = 0.05, 
                                   dataCol = "Fraction", effectCol = "Odds Ratio", 
                                   statUnit = "Patient"){

    sheet <- "all_fractions_and_densities"
    signifCol <- paste(dataCol, "adjusted p.value")
    pvalCol   <- paste(dataCol, "p.value")
    overall <- read.xlsx(overallFile, sheet, check.names = F) %>% as_tibble() %>%
               mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
               filter(`Cell State ID` %in% idsToPlot) 

    countCols <- names(overall)[grepl("total cell count", names(overall))]
    medCols   <- names(overall)[grepl(paste("median",tolower(dataCol)), names(overall))]

    overall <- overall %>%
               select(`Cell State ID`, `Cell State`, Population, 
                      OverallSubpopulationCount = !!as.name(countCols[1]), 
                      OverallPopulationCount = !!as.name(countCols[2]),
                      OverallGroup1Median = !!as.name(medCols[1]), 
                      OverallGroup2Median = !!as.name(medCols[2]),
                      OverallEffect = !!as.name(effectCol), 
                      OverallPval = !!as.name(pvalCol), 
                      OverallFDR = !!as.name(signifCol))

    
    match_status <- function(indivEffect, allEffect, indivFDR, max.fdr = 0.05){
        lapply(1:length(indivEffect), function(x){
            if(is.na(indivEffect[x])){ return("stats unavailable") }
            if(is.na(allEffect[x])){ return(NA) }

            if(all(c(indivEffect[x], allEffect[x]) > 1) || all(c(indivEffect[x], allEffect[x]) < 1)){
                ifelse(indivFDR[x] < max.fdr, 
                       "LOR same direction, FDR < 0.05",
                       "LOR same direction, FDR \u2265 0.05")
            } else {
                "LOR opposite direction"
            }

         }) %>% unlist()
    }

    get_value <- function(indivEffect, allEffect){
        lapply(1:length(indivEffect), function(x){

            if(is.na(indivEffect[x])){ return(NA) }
            if(is.na(allEffect[x])){ return(NA) }

             if(all(c(indivEffect[x], allEffect[x]) > 1) || all(c(indivEffect[x], allEffect[x]) < 1)){
                 indivEffect[x] 
             } else {
                 -1
             }

        }) %>% unlist()                 
    } 

    dat <- tibble()
    for(ind in indivFiles){
        #print(ind)
        id <- gsub(".xlsx", "", basename(ind))
        prfx <- "Indiv"
        indivEffect <- "IndivEffect" 
        indivSignif <- "IndivFDR" 
        indivPval   <- "IndivPval"
        overallEffect <- "OverallEffect"

        stats <- read.xlsx(ind, sheet, check.names = F) %>%
                 as_tibble() %>%
                 mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
                 filter(`Cell State ID` %in% idsToPlot) %>%
                 rename(!!as.name(indivEffect) := effectCol,
                        !!as.name(indivPval) := pvalCol,
                        !!as.name(indivSignif) := signifCol) %>%
                 left_join(overall, by = intersect(names(overall), names(.))) %>%
                 mutate(!!as.name(statUnit) := id,
                        Status = match_status(!!as.name(indivEffect),
                                              !!as.name(overallEffect),
                                              !!as.name(indivSignif),
                                              max.fdr = max.fdr),
                        Value = get_value(!!as.name(indivEffect),
                                          !!as.name(overallEffect)),
                        signif = ifelse(!!as.name(indivSignif) < max.fdr, "*", NA)
                 ) %>%
                 select_at(c("Cell State ID", "Cell State", "Population", 
                             names(overall)[grepl("Overall", names(overall))], 
                             indivEffect, indivPval, indivSignif, 
                             statUnit, "Status", "Value", "signif")) 
                         
        dat <- dat %>% bind_rows(stats)
    }
        
    dat$Status <- factor(dat$Status, 
                         levels = c("LOR same direction, FDR < 0.05", 
                                    "LOR same direction, FDR \u2265 0.05",
                                    "LOR opposite direction",
                                    "stats unavailable"))
    dat

}

formatStatsForPlottingEffects <- function(statsFile, sheet, calcType, dataCol, allConds, ids, 
                                          facetY = NULL, facetOrder = NULL, cellTypes = NULL,
                                          popOrder = NULL, orderBy = NULL, idOrder = NULL){

    cnds <- setConditionOrder(conds, ids, cellTypes, orderBy = orderBy, 
                              statsFile = statsFile, sheet = sheet, idOrder = idOrder,
                              facetY = facetY, facetOrder = facetOrder)

    effect <- ifelse(dataCol == "Fraction", "Odds Ratio", "Fold Change")
    stats <- xlsx::read.xlsx(statsFile, sheetName = sheet, check.names = F) %>%
             as_tibble() %>%
             select(`Cell State ID`, `Cell State`, Population, dplyr::matches(dataCol), 
                    !!as.name(effect), paste0(dataCol," adjusted p.value"), 
                    paste0(dataCol, " CI.low"), paste0(dataCol," CI.high")) %>%
             mutate(Condition = `Cell State`, `Cell State ID` = as.numeric(`Cell State ID`)) %>%
             filter(`Cell State ID` %in% ids)

    if(is.null(stats)){ stop("Couldn't find appropriate sheet in stats file") }

    names(stats) <- gsub(paste0("^",dataCol," "),"",names(stats))

    if(nrow(stats) <= 1){  
        log_warn("    WARNING: no statistics available for this question.")
        return(NULL) 
    }

    ## format y labels
    stats <- stats %>%
             unique() %>%
             ungroup() %>%
             select(-Condition, -Population) %>%   ## do this because they may have changed for facetting purposes
             mutate(signif = signifStars(`adjusted p.value`)) 

    stats %>% left_join(cnds, by = intersect(names(.), names(cnds)))
}

formatFOVdataForPlottingDetail <- function(question, sGrps, dat, condsToPlot, allConds, analyses, 
                                           calcType, dataCol, markers, metricsDir, orderBy = NULL,
                                           statsFile = NULL, sheet = NULL,
                                           nbhdCounts = NULL, tumorNbhdCells = NULL, calcUnit = "FOV_ID",
                                           facetY = NULL, facetOrder = NULL, cellTypes = NULL,
                                           popOrder = NULL, idOrder = NULL){

    condsToPlot <- as.numeric(condsToPlot)
    cnds <- setConditionOrder(allConds, condsToPlot, cellTypes, orderBy = orderBy,
                              statsFile = statsFile, sheet = sheet, idOrder = idOrder, 
                              facetY = facetY, facetOrder = facetOrder)

    qDat <- suppressMessages(getQuestionData(question, sGrps, dat, analyses, calcType, markers, metricsDir,
                                             nbhdCounts = nbhdCounts, calcUnit = calcUnit,
                                             tumorNbhdCells = tumorNbhdCells)) 
    if(is.null(qDat)){
        log_warn("   WARNING: One or more groups in this question contains no data. Skipping.\n")
        return(NULL)
    }

    qcQuestionGroups(qDat, names(analyses), question$groupVar)

    log_debug("Calculating Q1,Q3,Median & IQR")

    maxx <- ifelse(calcType == "fractions", 1, 10^4) 
    minx <- ifelse(calcType == "densities", 10^-1, 0)

    selCols <- c("Cell State ID", facetY[facetY %in% names(.)], "Condition", "Population",
                 question$groupVar, calcUnit, dataCol)

    fovDat <- qDat[[calcType]] %>% 
              ungroup() %>%
              mutate(`Cell State ID` = as.numeric(`Cell State ID`)) %>%
              filter(`Cell State ID` %in% condsToPlot) %>%
              select_at(selCols[selCols %in% names(.)]) %>%
              group_by_at(c("Cell State ID", question$groupVar)) %>%
              mutate(Median = median(!!as.name(dataCol), na.rm=T),
                     Q1 = summary(!!as.name(dataCol), na.rm=T)[2],
                     Q3 = summary(!!as.name(dataCol), na.rm=T)[4],
                     IQR = Q3 - Q1,
                     min = max(Q1 - (1.5*IQR), minx, na.rm=T),
                     max = min(Q3 + (1.5*IQR), maxx, na.rm=T))

    if(!"Condition" %in% names(fovDat)){
        fovDat <- fovDat %>% mutate(Condition = Population)
    }

    ## adjust group labels
    grpLbls <- fovDat %>% 
               ungroup() %>%
               select_at(c(calcUnit, question$groupVar)) %>%
               unique() %>%
               group_by_at(quest$groupVar) %>%
               summarize(Count = n()) %>%
               mutate(GroupLabel = paste0(!!as.name(quest$groupVar), " (FOVs=", Count, ")")) %>%
               select(-Count)
    grpLbls$GroupLabel <- factor(grpLbls$GroupLabel, levels = unique(grpLbls$GroupLabel))

    fovDat <- fovDat %>% 
              left_join(grpLbls, by = c(question$grpVar)) %>% 
              ungroup() %>%
              select(-Condition, -Population) %>%   ## do this because they may have changed for facetting purposes
              left_join(cnds, by = intersect(names(.), names(cnds)))

    fovDat

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

#' Format cell type data to be used as condition annotation
#' 
#' Format cell type data to be used as condition annotation
#' 
#' @param allData               table of all data to be plotted, including columns Cell_type, Condition, Category,
#'                              and any columns included in annotation_order
#' @param duplicate_conds_file  file to which all duplicate conditions will be written (as a warning)
#' @param annotation_order      order in which annotation column ordering should be done 
#' @param functional_order      order in which Functional values should appear in data
#' @param cell_types_order      order in which Cell_type values should appear in data
#' @param subtypes_order        order in which Subtype values should appear in data
#' @param combine_all_functional collapse 'Functional' and 'Functional Combinations' into one column, since no
#'                               condition should have both values
#' @return annotation data sorted as specified
formatAnnotationData <- function(allData, duplicate_conds_file, annotation_order = NULL, functional_order = NULL,
                                 cell_types_order = NULL, subtypes_order = NULL, combine_all_functional = FALSE){

    ## pull out annotation data
    #annotCols <- c("Cell State ID", "Cell_type", "Subtype", "Abbrev_label_for_figures", "Number Significant", cfg$annotion_order)

    annot <- allData %>%
             filter(!is.na(Cell_type)) %>%
             select(unique(c("Condition", "Category", annotation_order))) %>%
             unique() %>%
             ungroup()

    ## warn of duplicates 
    duplicateConditions(annot, annotation_order, outFile = duplicate_conds_file)

    annotCols <- c()

    if(!is.null(annotation_order)){
        for(ann in annotation_order){
            if(ann == "Cell State ID"){
                annotCols <- c(annotCols, "as.numeric(`Cell State ID`)")
            } else if(ann == "Functional"){
                funcOrder <- sapply(sort(functional_order), function(x){ which(functional_order == x) })
                annot$FuncOrder <- sapply(annot$Functional, function(x){ funcOrder[x] })
                annotCols <- c(annotCols, "FuncOrder")
            } else if(ann == "Functional Combination"){
                annot$NumPos <- sapply(annot$`Functional Combination`, function(x){ numPositives(x) })
                annotCols <- c(annotCols, "desc\\(NumPos\\)")
            } else if(ann == "Cell_type"){
                ctOrder <- sapply(sort(cell_types_order), function(x){ which(cell_types_order == x) })
                annot$CellTypesOrder <- sapply(annot$Cell_type, function(x){ ctOrder[x] })
                annotCols <- c(annotCols, "CellTypesOrder")
            } else if(ann == "Subtype"){
                stOrder <- sapply(sort(cell_types_order), function(x){ which(subtypes_order == x) })
                annot$SubtypesOrder <- sapply(annot$Subtype, function(x){ stOrder[x] })
                annotCols <- c(annotCols, "SubtypesOrder")
            } else {
                annotCols <- c(annotCols, ann)
            }
        }
    }

    annotCols <- sapply(annotCols, function(x){ ifelse(grepl(" ",x) && !grepl("`", x), add_ticks(x), x) }) %>%
                        as.vector()

    ## SORT
    annot <- annot %>% arrange_(annotCols)

    ## if functional cols to be combined, do this AFTER sorting
    if(combine_all_functional){
        ## error if any overlap btwn the two cols
        overlap <- intersect(which(!is.na(annot$FuncOrder)), which(!is.na(annot$NumPos)))
        if(length(overlap) > 0){ stop("Can't combine single functional markers and marker combinations.") }
        idxs <- which(is.na(annot$Functional) & !is.na(annot$`Functional Combination`))
        annot$Functional[idxs] <- annot$`Functional Combination`[idxs]
        annotation_order <- annotation_order[!annotation_order == "Functional Combination"]
    }

    ## warn of any conditions dropped after sorting
    missingConditions(annot, annot %>% pull(Condition))

    ## convert to matrix for Heatmap()
    annDat <- annot %>%
              select_at(c("Condition", annotation_order[annotation_order != "Abbrev_label_for_figures"])) %>%
              data.frame() %>%
              toMatrix(annot %>% pull(Condition))

    if("Cell_type" %in% names(annDat)){
        annDat[,"Cell_type"] = unlist(sapply(annDat[,"Cell_type"], function(x){
                                  paste(strwrap(x, width=15), collapse="\n")
                               }))
    }
    annOrder <- gsub(" ",".",gsub("`","",annotation_order))

    figAnn <- annDat[,annOrder]

    figAnn

}


sortRows <- function(dat, ctList, arrangeBy){
    rws <- list(order = c(), ctBlocks = c())

    dat$Abbrev_label_for_figures[which(dat$Cell_type == dat$Subtype | is.na(dat$Subtype))] <- "All"

    for(ct in names(ctList)){
        subtypes <- ctList[[ct]]
        for(st in subtypes){

            if(is.na(st)){
                stConds <- dat %>% filter(Cell_type == ct, is.na(Subtype))
            } else {
                stConds <- dat %>% filter(Cell_type == ct, Subtype == st)
            }
            stConds <- stConds %>%
                       arrange_(.dots = arrangeBy) %>%
                       select(Condition, Abbrev_label_for_figures)
            stConds$Abbrev_label_for_figures[is.na(stConds$Abbrev_label_for_figures)] <- "All"
            #### TEMPORARY::
            ct2 <- ifelse(ct == "Natural killer cell overall", "Natural killer cell", ct)
            rws$ctBlocks <- c(rws$ctBlocks, rep(paste(strwrap(ct2,15), collapse="\n"), nrow(stConds)))
            rws$order <- c(rws$order, stConds$Condition)
        }
    }
    rws$ctBlocks <- factor(rws$ctBlocks, levels = unique(rws$ctBlocks))
    rws
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


formatMedians <- function(dat, c2p, groups, calcCol, layout){
    meds <- dat %>%
            filter(`Cell State ID` %in% c2p) %>%
            select(`Cell State ID`, dplyr::matches("Overall Median")) %>%
            rename_at(vars(names(.)[grepl("Overall Median", names(.))]), 
                      list(~ gsub(" Overall Median.*", "", .))) %>%
            gather(2:ncol(.), key = "Group", value = "raw median") %>%
            left_join(layout, by = "Cell State ID") %>%
            group_by(`Cell State ID`) %>%
            mutate(maxMed = max(`raw median`, na.rm=T),
                   median = `raw median`/maxMed,
                   y2 = as.numeric(y),
                   bottom = y2 - 0.5,
                   top = (y2 - 0.5) + median*0.85)

    medIdxs <- grep("Overall Median", names(meds))
    names(meds)[medIdxs] <- gsub(" Overall Median.*", "", names(meds)[medIdxs]) 


    meds$x     <- customOrder(meds$Group, grps)
    meds$Group <- factor(meds$Group, levels = grps)
    meds$y     <- factor(as.numeric(meds$y), levels = seq(0.5, max(as.numeric(meds$y), na.rm=T) + 0.5, 0.5))

    meds
}

formatIndivDatForPlotting <- function(dat, idMap, layout){
    iDat <- dat %>%
            left_join(idMap, by = intersect(names(.), names(idMap))) %>%
            left_join(layout, by = intersect(names(.), names(layout))) %>%
            arrange(desc(labelY)) %>%
            mutate(status = ifelse(Value == -1, "opposite",
                                   ifelse(Value > 1, "up",
                                          ifelse(Value < 1, "down",
                                                 ifelse(Value == 1, "no change", NA)))))
    iDat$x <- customOrder(iDat$`Lesion ID`, idMap$`Lesion ID`[idMap$`Lesion ID` %in% iDat$`Lesion ID`])
    iDat$x <- factor(iDat$x, levels = sort(unique(iDat$x)))

    iDat
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


#' Format raw neighborhood counts table for statistical analyses
#' 
#' Filter table for cells in region of interest and join center cell annotation
#' including Sample_ID, Band and PositiveMarkers
#' 
#' @param nbhdCountsDir    directory of files each containing table of neighborhood counts 
#'                         including columns 'C.UUID' and 'FOV'
#' @param dlDat            annotated cell data for ONLY cells inside region of interest
#' @param cellDiveID       CellDive_ID for which data should be formatted
#' 
#' @return  filtered & annotated neighborhood counts table
formatNeighborhoodCounts <- function(nbhdCountsDir, dlDat, cellDiveID = "All"){

    ncFiles <- file.path(nbhdCountsDir, dir(nbhdCountsDir))
    if(tolower(cellDiveID) != "all"){ 
        ncFiles <- ncFiles[grepl(paste0("^", cellDiveID, "_"), basename(ncFiles))] 
    }

    if("Band" %in% names(dlDat)){
        dlDat <- dlDat %>% select(Band, C.UUID = UUID, C.PositiveMarkers = PositiveMarkers)
    } else {
        dlDat <- dlDat %>% select(C.UUID = UUID, C.PositiveMarkers = PositiveMarkers)
    }

    allNbhdCounts <- tibble()

    for(ncFile in ncFiles){
        log_debug(paste0("Formatting nbhd counts in file: ",ncFile))
        nCts <- readRDS(ncFile)
        if(!"FOV_ID" %in% names(nCts)){
            nCts <- nCts %>% rename(FOV_ID = FOV)
        }

        FOV_IDs <- unique(nCts$FOV_ID)

        ## fill in missing neighborhood cell types for 
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

        if(!all(FOV_IDs %in% names(tmp))){
            for(f in FOV_IDs[!FOV_IDs %in% names(tmp)]){
                tmp[[f]] <- 0
            }
        }

        tmp <- tmp %>%
               gather(FOV_IDs, key="FOV_ID", value="N.Count") %>%
               filter(N.Count == 0)

        allNbhdCounts <- allNbhdCounts %>%
                         bind_rows(nCts %>%
                                   bind_rows(tmp) %>%
                                   left_join(dlDat, by=c("C.UUID")))

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

