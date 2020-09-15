#' Run statistical test on one or more columns of data in a table
#'
#' Given a table containing a column or columns of groups (e.g., 'Lesion Response') and 
#' one or more columns of data, run specified statistical test on one or more
#' of the data columns
#' 
#' @param fun        statistical function to run (default=t.test)
#' @param df         data table
#' @param vars       names of data column(s) on which to run the statistical test
#' @param groupVar   name of column containing the groups each value in data column belongs
#'                   to; column must contain exactly two unique values and will be converted
#'                   to a factor whose levels are sorted alphanumerically
#'
#' @return a list of stats results where each element contains results for one item in {vars}
#' @export
multi.tests <- function(fun = t.test, df, vars, groupVar, ...) {
    if(!is.factor(df[[groupVar]])){
        df[[groupVar]] <- factor(df[[groupVar]], levels=unique(sort(df[[groupVar]])))
        log_warn(paste0("Factor levels and order not set for ", groupVar,". Setting order to : ", levels(df[[groupVar]])))
    }
    groupVar <- paste0("`",groupVar,"`")

    lst <- list()
    fails <<- 0
    for(x in 1:length(vars)){
        var <- vars[x]
        newName <- paste0("`",var,"`")

        formula <- as.formula(paste(newName, "~", groupVar))
        lst[[var]] <- tryCatch({
                         fun(data = df, formula, conf.int=TRUE)
                      }, error = function(e){
                         fails <<- fails + 1
                         log_warn(paste0("[",fails,"] Could not calculate stats for ",var))
                         #print(e)
                         NULL
                      })
    }
    return(lst)
}


#' Format numbers for XLSX output
#' 
#' Format numbers for XLSX output
#' 
#' @param df   data frame to be formatted
#'
#' @return data frame with all numeric values formatted
formatFinalNumbers <- function(df){
    types <- sapply(df, class)

    ints <- unique(c(names(df)[grepl("count",names(df),ignore.case=TRUE)],
                     names(types)[which(types == "integer")]))
    nums <- unique(c(names(df)[grepl("median|density$|area|estimate|low|high|p\\.value|ratio|mean", 
                               names(df), 
                               ignore.case=TRUE)],
                     names(types)[which(types == "numeric")]))

    nums <- nums[!nums %in% ints]
    for(n in nums){
        df[[n]] <- as.numeric(df[[n]])
        df[[n]] <- as.numeric(formatC(df[[n]], digits=3))
    }

    return(df)
}

#' Test that sample manifests from different stats calculations are identical
#' 
#' When running stats on different types of calculations, each set of results
#' includes a sample manifest indicating exactly which samples were used in the
#' calculations. Manifests for all types of calculations run for a single 'question'
#' should be identical.
#' 
#' @param allManifests  list of manifests from different calculation types
#'
#' @return logical indicating whether manifests are identical
identicalManifests <- function(allManifests){
    combined <- unique(bind_rows(allManifests))
    for(calc in names(allManifests)){
        if(!identical(allManifests[[calc]], combined)){ return(FALSE) }
    }
    return(TRUE)
}

#' Write stats results for a single question and a single cell region to XLSX file
#' 
#' Write stats results for a single question and a single cell region to XLSX file
#' 
#' @param tblList                list of stats results where each item is a table of results 
#'                               for one calculation type (e.g., densities or fractions), and/or 
#'                               sample information; each item will be written to its own tab
#' @param fileName               XLSX file to which results will be saved
#' @param filterName             name of filter used (simply for naming tabs of filtered data)
#'
#' @return nothing
writeStatsQuestionXLSX <- function(tblList, fileName){
    countStyle <- openxlsx::createStyle(numFmt = 'COMMA')
    wb = openxlsx::createWorkbook()
    for(i in seq(tblList)){
        tblName <- names(tblList)[i]
        tbl <- tblList[[tblName]]
        if(nrow(tbl) > 0){
            tbl <- formatFinalNumbers(tbl)
        }
        openxlsx::addWorksheet(wb, tblName)
        openxlsx::writeData(wb, i, tbl)
        openxlsx::freezePane(wb, tblName, firstRow = TRUE)
        countColIdxs <- grep("count", names(tbl), ignore.case=T)
        if(length(countColIdxs) > 0){
            for(cc in countColIdxs){
                openxlsx::addStyle(wb, i, style=countStyle, cols=cc, rows=1:nrow(tbl)+1)
            }
        }
    }
    openxlsx::saveWorkbook(wb, fileName, overwrite=TRUE)
}


#' Adjust 1s and 0s to force them to be included by wilcox.test()
#' 
#' Add a very small number to 0s or subtract it from 1s so that they
#' can be converted to log scale and properly handled by the wilcox test
#' 
#' @param dat  data with a column for the groups to be compared (e.g., 'Treated' or
#'             'Patient Response') and one for Fraction
#' @param groupVar  character string; column name of groups to be compared
#' @param valueCol  column in dat to be scaled
#'
#' @return dat with column 'psValue' with scaled data values
pseudoScale <- function(dat, groupVar, valueCol){
    tmp <- dat
    ## assume fractions at first
    tmp$OM <- 1 - tmp[[valueCol]]
    tmp <- tmp %>%
           ungroup() %>%
           group_by_at(groupVar) %>%
           filter_at(all_of(valueCol), any_vars(. > 0 & . < 1))
    if(nrow(tmp) > 0){
        ps <- tmp %>%
              summarize_at(vars(valueCol,"OM"), min) %>%
              select_at(all_of(c(valueCol,"OM"))) %>%
              min()/2

        scaled <- dat[[valueCol]]
        scaled[scaled == 0] <- ps
        scaled[scaled == 1] <- 1 - ps
    } else {  ## nothing in between 0 and 1, assume densities
        ps <- (dat %>%
               ungroup() %>%
               filter_at(all_of(valueCol), any_vars(. > 0)) %>%
               summarize_at(vars(valueCol), min) %>%
               pull(valueCol) - 1)/2
        scaled <- dat[[valueCol]]
        scaled[scaled == 0] <- ps
    }
    scaled
}


#' Create table of neighborhood average counts statistics
#'
#' Given a table (or tables) of average counts of neighborhood types for each center cell in 
#' two sample groups, run wilcox.test to get odds ratios and conficence intervals for each neighborhood
#' type
#'
#' @param avgDat      table of data with following columns at minimum: 
#'                     Population, Samp1, Samp2, ..., SampN
#'                        OR
#'                     Population, Sample, AvgCount, [Group]
#' @param allSamples  a vector containing the union of group1 samples and group2 samples ONLY
#' @param groupVar    string representing name of variable to be compared; will be either a 
#'                    column name from sampAnnot or arbitrary string when comparing non-clinical 
#'                    variable (e.g., inside/outside, pos/neg microenv)
#' @param countDat    table of neighborhood counts with at least columns for calcUnit, 'CenterCellType',
#'                    'NeighborhoodCellType', and groupVar
#' @param calcUnit    data column name containing values for which a single area value
#'                    should be calculated (default: Sample)
reportAverageCountsStats <- function(avgDat, allSamples, groupVar, countDat, calcUnit = "FOV_ID"){

    dat <- avgDat %>%
           mutate(CondTitle = nbhdAvgCond(CenterCellType, NeighborhoodCellType)) %>%
           removeDuplicateRows(c("CondTitle", calcUnit, groupVar))

    dat$psCount <- pseudoScale(dat, groupVar, "AvgNbhdCellTypeCount")

    logCounts <- dat %>%
                 mutate(logAvgCount = log(psCount)) %>%
                 select_at(c(groupVar, calcUnit, "CondTitle", "logAvgCount")) %>%
                 spread(CondTitle, logAvgCount)

    conds <- dat %>% 
             select(CondTitle, CenterCellType, NeighborhoodCellType, `Cell State ID`) %>% 
             unique()
    stats <- logCounts %>% getStatsFC(groupVar, unique(dat$CondTitle))
    counts <- countDat %>% getNeighborhoodCountsSummary()
    meds   <- dat %>% getGroupMedianNeighborhoodAvgs(groupVar)

    conds %>%
    full_join(counts, by = c("CenterCellType", "NeighborhoodCellType")) %>%
    full_join(meds, by = "CondTitle") %>%
    full_join(stats, by = "CondTitle") %>%
    fixPvals(., paste(unique(dat[[groupVar]]), "median")) %>%
    ungroup() %>%
    select(`Cell State ID`, 
           `Center Population` = CenterCellType, 
           `Neighborhood Population` = NeighborhoodCellType,
           `Total Center Cell Count` = Total.C.Count,  
           `Total Neighborhood Cell Count` = Total.N.Count, 
           `Overall mean count`,
           contains(" median"), 
           `log(Fold Change)`, `Fold Change`, `Inverse Fold Change`,
           CI.high, CI.low, p.value, `adjusted p.value`) %>%
    arrange(as.numeric(`Cell State ID`))

}


getStatsLO <- function(loDat, groupVar, varCols, func = wilcox.test){
    multi.tests(fun = func, df = loDat, vars = varCols, groupVar = groupVar) %>%
    wilcoxResToTable(., log2exp = TRUE) %>%
    dplyr::rename(`Odds Ratio` = estimate) %>%
    mutate(`log(Odds Ratio)` = log(as.numeric(`Odds Ratio`)),
           `Inverse Odds Ratio` = 1/as.numeric(`Odds Ratio`),
           `adjusted p.value` = p.adjust(`p.value`, method = "bonferroni")) 
}

getStatsFC <- function(logDat, groupVar, varCols, func = wilcox.test){
    multi.tests(fun = func, df = logDat, vars = varCols, groupVar = groupVar) %>% 
    wilcoxResToTable(., log2exp = TRUE) %>%
    dplyr::rename(`Fold Change` = estimate) %>%
    mutate(`log(Fold Change)` = log(as.numeric(`Fold Change`)),
           `Inverse Fold Change` = 1/as.numeric(`Fold Change`),
           `adjusted p.value` = p.adjust(`p.value`, method = "bonferroni"))
}

getTotalCountsForLOReport <- function(allCountDat){
    allCountDat %>%
    group_by(`Cell State ID`, Condition, Population, CondTitle) %>%
    summarize(`Cell State total cell count` = sum(CondCount, na.rm=T),
              `Population total cell count` = sum(PopCount, na.rm=T)) %>%
    unique()
}

getNeighborhoodCountsSummary <- function(allCountDat){
    centers <- allCountDat %>%
               select(C.UUID, CenterCellType) %>%
               unique() %>%
               group_by(CenterCellType) %>%
               summarize(Total.C.Count = n())
    allCountDat %>%
    left_join(centers, by = "CenterCellType") %>%
    group_by(CenterCellType, NeighborhoodCellType, Total.C.Count) %>%
    summarize(Total.N.Count = sum(N.Count, na.rm = T),
              `Overall mean count` = mean(N.Count, na.rm = T)) %>%
    ungroup()
}

getTotalCountsForNeighborhoodLOReport <- function(allCountDat, conds){
    counts <- allCountDat %>% getNeighborhoodCountsSummary()

    conds %>%
    left_join(counts %>%
              select(`Center Population A` = CenterCellType,
                     `Neighborhood Population A` = NeighborhoodCellType,
                     Total.C.Count, Total.N.Count),
              by = c("Center Population A", "Neighborhood Population A")) %>%
    rename(`Neighborhood Population A total cell count` = Total.N.Count) %>%
    left_join(counts %>%
              select(`Center Population B` = CenterCellType,
                     `Neighborhood Population B` = NeighborhoodCellType,
                     Total.C.Count, Total.N.Count),
              by = c("Center Population B", "Neighborhood Population B", "Total.C.Count")) %>%
    rename(`Center cell type total count` = Total.C.Count,
           `Neighborhood Population B total cell count` = Total.N.Count)    
}

getTotalCountsForDensityFCReport <- function(dat){
    dat %>%
    select(`Cell State ID`, CondTitle, Count, TotalArea, Density) %>%
    group_by(`Cell State ID`, CondTitle) %>%
    summarize(`Cell State total cell count` = sum(Count, na.rm = T),
              `Total Area` = sum(TotalArea, na.rm = T),
              `Overall Median Density` = median(Density, na.rm = T))
}

getGroupMedianFractions <- function(allCountDat, groupVar, fracCol){
    groupNames <- unique(as.vector(allCountDat[[groupVar]]))

    allCountDat %>%
    group_by_at(c("CondTitle", groupVar)) %>%
    summarize(Median = median(!!as.name(fracCol), na.rm=T)) %>%
    spread_(groupVar, "Median", fill = 0) %>%
    dplyr::rename_at(all_of(groupNames), list(~ paste(., "median fraction")))
}

getGroupMedianDensities <- function(dat, groupVar){
    groupNames <- unique(as.vector(dat[[groupVar]]))

    dat %>%
    group_by_at(c("Cell State ID", "CondTitle", groupVar)) %>%
    summarize(Median = median(Density, na.rm=T)) %>%
    spread(groupVar, Median, fill=0) %>%
    dplyr::rename_at(all_of(groupNames), list(~ paste(., "median density")))
}

getGroupMedianNeighborhoodAvgs <- function(dat, groupVar){
    groupNames <- unique(as.vector(dat[[groupVar]]))

    dat %>%
    group_by_at(c(groupVar, "CondTitle")) %>%
    summarize(GroupMedian = median(AvgNbhdCellTypeCount, na.rm = T)) %>%
    spread(groupVar, GroupMedian, fill=0) %>%
    dplyr::rename_at(all_of(groupNames), list(~ paste(., "median")))
}

fractionsToLO <- function(fracDat, fracCol, groupVar){
    dat <- fracDat
    dat$psFraction <- pseudoScale(dat, groupVar, fracCol)
    dat %>% mutate(LO = logOdds(psFraction)) 
}

fractionCond <- function(conds, pops){
    lapply(1:length(conds), function(x){
        gsub("/| ","_",paste(conds[x], pops[x], sep="___"))
    }) %>% 
    unlist()
}

nbhdAvgCond <- function(center, nbhd){
    paste0(center, "___", nbhd)
}

nbhdFractionCond <- function(centerA, nbhdA, centerB, nbhdB){
    paste0(centerA, "_", nbhdA, "__vs__", centerB, "_", nbhdB)
}

#' Run wilcox test and format log odds statistics 
#' 
#' Run wilcox test and format log odds statistics
#'
#' @param fracDat      pre-filtered table of fraction data with, at minimum, columns: Condition, Population,
#'                     calculation unit ("FOV_ID" by default), Fraction and groupVariable (e.g., Patient_response) 
#' @param groupVar     string representing name of variable to be compared; will be either a column name from sampAnnot
#'                     or arbitrary string when comparing non-clinical variable (e.g., inside/outside, pos/neg microenv)
#' @param calcUnit     data column name containing values for which a single area value
#'                     should be calculated (default: Sample)
#' @param max.fdr      maximum adjusted p.value to consider a condition 'Passed'
#' @param min.or       minimum odds ratio to consider a condition 'Passed'
#' @param  full fraction statistics table
reportLogOdds <- function(fracDat, groupVar, calcUnit = "FOV_ID",
                          max.fdr = 0.05, min.or = 2, fracCol = "Fraction"){

    dat <- fracDat %>% mutate(CondTitle = fractionCond(Condition, Population)) 

    ## transform fractions to log odds for stats
    loTbl <- dat %>% 
             fractionsToLO(fracCol, groupVar) %>%
             select_at(c("CondTitle", calcUnit, groupVar, "LO")) %>%
             spread(CondTitle, LO)
    
    stats  <- loTbl %>% getStatsLO(groupVar, unique(dat$CondTitle))
    counts <- dat %>% getTotalCountsForLOReport()
    meds   <- dat %>% getGroupMedianFractions(groupVar, fracCol)

    ## compile report
    counts %>%
    left_join(meds, by = intersect(names(.), names(meds))) %>%
    left_join(stats, by = intersect(names(.), names(stats))) %>%
    fixPvals(., paste(unique(dat[[groupVar]]), "median fraction")) %>%
    mutate(`Fraction Passed` = ifelse(`adjusted p.value` < max.fdr &
                                      abs(`log(Odds Ratio)`) > log(min.or),
                                          "X", "")) %>%
    select(`Fraction Passed`, `Cell State ID`, `Cell State` = Condition, Population,
           `Cell State total cell count`, `Population total cell count`,
           contains("median fraction"),
           `log(Odds Ratio)`, `Odds Ratio`, `Inverse Odds Ratio`,
           everything()) %>%
           select(-CondTitle) %>%
           arrange(as.numeric(`Cell State ID`)) %>%
           ungroup()

}

#' Run wilcox test and format log odds statistics 
#' 
#' Run wilcox test and format log odds statistics
#'
#' @param fracDat      pre-filtered table of data with following columns at minimum: 
#'                     Condition, CondTitle, MainCount, SubCount, Samp1, Samp2, ..., SampN
#'                     Alternate format: Group, ContTitle, Subpopulation, Population, SubCount, MainCount, Sample, Fraction
#' @param countDat     tibble of counts for each population within each calcUnit
#' @param allSamples   a vector containing the union of group1 samples and group2 samples ONLY
#' @param groupVar     string representing name of variable to be compared; will be either a column name from sampAnnot
#'                     or arbitrary string when comparing non-clinical variable (e.g., inside/outside, pos/neg microenv)
#' @param calcUnit     data column name containing values for which a single log odds value
#'                      should be calculated (default: FOV_ID)
#' 
#' @return tibble containing population count medians in each group, log odds, confidence intervals, 
#'         p.value and FDR 
reportNeighborhoodLogOdds <- function(fracDat, countDat, allSamples, groupVar, calcUnit = "FOV_ID",
                                      fracCol = "Fraction"){ 

    dat <- fracDat %>%
           mutate(CondTitle = nbhdFractionCond(`Center Population A`, `Neighborhood Population A`,
                                               `Center Population B`, `Neighborhood Population B`))

    loTbl <- dat %>%
             fractionsToLO(fracCol, groupVar) %>%
             select_at(c("CondTitle", calcUnit, groupVar, "LO")) %>%
             spread(CondTitle, LO)

    ## pull out IDs and condition parts
    conds <- dat %>% select(CondTitle, `Cell State ID`,
                            `Center Population A`, `Neighborhood Population A`,
                            `Center Population B`, `Neighborhood Population B`) %>%
                            unique()

    stats  <- loTbl %>% getStatsLO(groupVar, unique(dat$CondTitle))
    counts <- countDat %>% getTotalCountsForNeighborhoodLOReport(conds)
    meds   <- dat %>% getGroupMedianFractions(groupVar, fracCol)

    counts %>% 
    left_join(meds, by = "CondTitle") %>%
    left_join(stats, by = intersect(names(.), names(stats))) %>%
    fixPvals(., paste(unique(dat[[groupVar]]), "median fraction")) %>%
    select(`Cell State ID`,
           `Center Population A`, `Neighborhood Population A`,
           `Center Population B`, `Neighborhood Population B`,
           `Center cell type total count`,
           `Neighborhood Population A total cell count`,
           `Neighborhood Population B total cell count`,
           contains("median fraction"),
           `log(Odds Ratio)`, `Odds Ratio`, `Inverse Odds Ratio`,
           everything()) %>%
    select(-CondTitle) %>%
    filter(!is.na(`Center Population B`)) %>%
    arrange(as.numeric(`Cell State ID`))

}

#' Remove rows of duplicate conditions from conditions index
#'
#' Conditions index contains some duplicates in the case of densities
#' and neighborhood averages as a result of joining with multiple
#' fractions conditions; Here we find and remove those duplicates
#'
#' @param dat         tibble that possibly contains duplicate conditions
#' @param uniqueCols  vector of column names required to distinguish unique conditions
#' 
#' @return tibble with duplicate conditions removed
removeDuplicateRows <- function(dat, uniqueCols){
    tmp <- dat %>% select_at(uniqueCols)
    if(any(duplicated(tmp))){
        log_warn(paste0("Removing ",length(which(duplicated(tmp))), 
                        " rows with duplicate conditions before analysis."))
        dupIDs <- sort(as.numeric(unique(dat$`Cell State ID`[which(duplicated(tmp))])))
        log_warn(paste("   IDs removed: ", paste(dupIDs, collapse=",")))
        dat <- dat[-which(duplicated(tmp)),]
    }
    dat
}

#' Run wilcox test and format log fold change statistics 
#' 
#' Run wilcox test and format log fold change statistics
#'
#' @param denDat       table of data with following columns at minimum: 
#'                     Population, Samp1, Samp2, ..., SampN
#'                        OR
#'                     Population, Sample, Count, Area, Density, [Group] 
#' @param allSamples   a vector containing the union of group1 samples and group2 samples ONLY
#' @param groupVar     string representing name of variable to be compared; will be either a column name from sampAnnot
#'                     or arbitrary string when comparing non-clinical variable (e.g., inside/outside, pos/neg microenv)
#' @param calcUnit     data column name containing values for which a single area value
#'                      should be calculated (default: Sample)
#' @param max.fdr      maximum adjusted p.value for 'passing' results
#' @param min.fc       minimum absolute fold change for 'passing' results
#' 
#' @return tibble of density stats including mean densities per group, fold changes, confidence intervals,
#'         p values and FDRs for all conditions
reportDensityStats <- function(denDat, allSamples, groupVar, calcUnit = "FOV_ID", 
                               max.fdr = 0.05, min.fc = 2){

    dat <- denDat %>% 
           rename(CondTitle = Population) %>%
           removeDuplicateRows(c("CondTitle", calcUnit, groupVar)) 

    dat$psDensity <- pseudoScale(dat, groupVar, "Density")
    dat <- dat %>% mutate(logDen = log(psDensity))

    logDens <- dat %>%
               select_at(c(groupVar, "CondTitle", calcUnit, "logDen")) %>%
               spread(CondTitle, logDen)

    counts <- dat %>% getTotalCountsForDensityFCReport()
    stats  <- logDens %>% getStatsFC(groupVar, unique(dat$CondTitle))
    meds   <- dat %>% getGroupMedianDensities(groupVar)

    counts %>%
    left_join(meds, by = intersect(names(.), names(meds))) %>%
    left_join(stats, by = intersect(names(.), names(stats))) %>%
    fixPvals(., paste(unique(dat[[groupVar]]), "median density")) %>%
    mutate(`Density Passed` = ifelse(`adjusted p.value` < max.fdr &
                                     abs(`log(Fold Change)`) > log(min.fc),
                                         "X", "")) %>%
    select(`Density Passed`, 
           `Cell State ID`, 
           `Cell State` = CondTitle, 
           `Cell State total cell count`, 
           `Total Area`, 
           `Overall Median Density`, 
           contains("median density"), 
           `log(Fold Change)`, 
           `Fold Change`, 
           `Inverse Fold Change`,
           everything()) %>%
    ungroup() %>%
    arrange(as.numeric(`Cell State ID`)) 

}

#' Assign sample group name based on all values of the comparison variable 
#' included in the group
#' 
#' Given filtering critera for a study question in list form and the group number, 
#' if there are multiple values to be included in a single group, name the group
#' by all of the values separated by a "/". For example, if a question is comparing
#' Patient_response and Group 1 includes both "Non-Responder" and "Partial Responder", 
#' the name will be "Non-Responder/Partial Responder"
#' 
#' @param quest     question description in list form (results of parseStatsQuestions()
#' @param groupNum  the group number to be assign a name (1 or 2)
#'
#' @return character string of group name
getGroupName <- function(quest, groupNum){
    paste(quest[[paste0("Group ",groupNum)]][[quest$groupVar]], collapse="/")
}


#' Divide sample IDs into two groups based on filtering criteria in list form
#'
#' Given a list \code{question} of filtering criteria, where list names match a column
#' in \code{sampAnnot} and values match one or more value(s) in that column
#' 
#' @param question   a single nested list with two inner lists: `Group 1` and `Group 2`; each inner list
#'                   contains filtering criteria for one group, based on information in \code{sampAnnot};
#'                   any list names NOT in \code{sampAnnot} will be assumed to be in data, which will be
#'                   filtered after sample filtering
#' @param sampAnnot  table of sample annotation where each row is a sample and each column contains values
#'                   corresponding to a single clinical variable
#' @param tumorNbhdCells vector of UUIDs for cells that are within 30um of one or more tumor cells
#' 
#' @return list of two items, each containing a vector of calc units belonging to one of the two
#'         groups of data to be compared
getQuestionGroups <- function(question, sampAnnot, calcUnit = "FOV_ID", tumorNbhdCells = NULL){
    grpVar <- question$groupVar
    smps <- list()

    for(g in c(1,2)){ 
        grpAnn <- sampAnnot
        grp <- question[[paste0("Group ",g)]]
        grpName <- getGroupName(question, g)
        for(filt in names(grp)){
            if(filt %in% names(grpAnn)){
                grpAnn <- grpAnn %>% 
                          filter_at(all_of(filt), any_vars(. %in% grp[[filt]]))
            } else if(filt == "Cell Region") {
                if(tolower(grp[[filt]]) %in% c("interface", "interface inside", "interface outside")){
                    grpAnn <- grpAnn %>% filter(FOV_type == "Interface")
                }
            }
        }
        smps[[grpName]] <- unique(grpAnn[[calcUnit]])
    }
    return(smps)
}



#' Check if all calculation units in a data set apply to 
#' a single combination of variables
#'
#' Determine whether a calculation unit (default FOV_ID) contains
#' a single value for a combination of variables
#' 
#' @param dat       tibble of data including at least columns for calculation
#'                  unit and the variable columns to check
#' @param varCols   vector of column names for which every calculation unit 
#'                  is expected to map to a single calculation unit
#' @param calcUnit  column name for which a single value was calculated
#'
#' @return logical; when TRUE, a single unique row was found in data for each
#'         calculation unit
singleValPerCalcUnit <- function(dat, varCols, calcUnit = 'FOV_ID'){
    dat %>% 
    select_at(c(calcUnit, varCols)) %>%
    unique() %>%
    group_by_at(calcUnit) %>%
    summarize(Count = n()) %>%
    filter(Count > 1) %>%
    nrow() == 0
}

#' Determine whether it is safe to use pre-calculated data for a 
#' particular calculation unit
#' 
#' Given a question description in list form and the raw cell data,
#' determine whether it is safe to use pre calculated data by finding
#' out if the combinations of filter parameters map to a single
#' calculation unit
#'
#' @param questionList list describing filtering criteria for a single question
#' @param dat          all cell level data to be used in calculations
#' @param calcUnit     column name containing values for which a single value
#'                     should be calculated
#' 
#' @return logical; when true, it is safe to use pre-calculated data for this calcUnit
useStoredData <- function(questionList, dat, calcUnit = 'FOV_ID'){
    varCols <- lapply(questionList, function(x) names(x)) %>%
               unlist(use.names=F) %>%
               unique()
    varCols <- intersect(varCols, names(dat))
    singleValPerCalcUnit(dat, varCols, calcUnit = calcUnit)
}


#' Filter out population data relative to a single question
#' 
#' Using group descriptions in the form of a list, for each group filter for data
#' that applies to only that group and add group assignments based either on data itself
#' or on sample annotation
#' 
#' @param question      parsed question data in the form of a list that describes how to filter data
#' @param smps          vector of all samples (both groups) to be included
#' @param analysisList  tibble containing all analyses to be run (may include fractions, densities, etc.)
#' @param markers       vector of all markers in analysis
#' @param metricsDir    output/rda directory where either previously-generated population *.rda files
#'                      were created or newly generated ones should be written
#' @param calcTypes     vector of calculation types to get ("fractions", "densities", etc.)
#' @param nbhdCounts        default = NULL; if running neighborhood analyses, this is the neighborhood data 
#'                      from Nick, modified to be on the correct dataLevel
#' @param calcUnit      data column name containing values for which a single area value
#'                      should be calculated (default: Sample)
#' @return list of filtered and group-assigned population data
getQuestionData <- function(question, smps, dat, analysisList, calcTypes, markers, 
                            metricsDir, calcUnit = "FOV_ID",
                            nbhdCounts=NULL, tumorNbhdCells = NULL){

    popDat <- list()

    dataStored <- useStoredData(question, dat, calcUnit = 'FOV_ID')
    usedUUIDs <- c()

    grpNames <- c()
    for(g in c(1,2)){
        grpName  <- getGroupName(question, g)
        grpNames <- c(grpNames, grpName)
        grpLabel <- paste0("Group ",g)
        log_debug(paste0("Filtering for group: ",grpName))

        cellRegion    <- getCellRegion(question[[grpLabel]]$`Cell Region`)
        metricsSubdir <- getMetricsSubDir(cellRegion, metricsDir) 
        if(!is.null(metricsSubdir)){
            log_debug(paste0("Checking for pre-computed metrics in ",metricsSubdir)) 
        }

        log_debug(paste0("Group dat dimensions before filtering: ",paste(dim(dat), collapse=" x ")))
        grpDat <- dat %>% 
                  filterForComparisonGroup(question[[grpLabel]], 
                                           tumorNbhdCells = tumorNbhdCells)

        ## make sure no uuids overlap with previously filtered group data
        if(any(intersect(grpDat$UUID, usedUUIDs))){ 
            msg <- "Some UUIDs have been included in multiple groups"
            log_error(msg)
            stop(msg)
        }

        usedUUIDs <- c(usedUUIDs, grpDat$UUID)

        areas <- NULL
        if("densities" %in% calcTypes){
            areas <- filterForAreas(grpDat, cellRegion, groupBy = calcUnit)
        }
        if(!all(grpDat[[calcUnit]] %in% smps[[grpName]])){
            log_warn("GROUP DATA INCORRECTLY FILTERED. Skipping.")
            return(NULL)
        } else if(!all(smps[[grpName]] %in% grpDat[[calcUnit]])){
            log_warn("Not all group samples found in data.")
        }
        log_debug(paste0("Group dat dimensions after filtering: ",paste(dim(grpDat), collapse=" x ")))
 
        grpPopDat <- getPopulationMetrics(grpDat, 
                                          analysisList, 
                                          markers, 
                                          areas      = areas, 
                                          outDir     = metricsSubdir,
                                          calcUnit   = calcUnit, 
                                          include    = smps[[grpName]], 
                                          nbhdCounts = nbhdCounts,
                                          numThreads = 6, 
                                          dataStored = dataStored)

        for(calc in names(grpPopDat)){
            grpPopDat[[calc]] <- grpPopDat[[calc]] %>% 
                                 mutate(group = grpName) %>%
                                 rename_at(vars(group), list(~ (. = question$groupVar)))
            if(is.null(popDat[[calc]])){
                popDat[[calc]] <- tibble()
            }
            popDat[[calc]] <- popDat[[calc]] %>% bind_rows(grpPopDat[[calc]])
        }

    }
    
    for(calc in names(popDat)){
        popDat[[calc]][[question$groupVar]] <- factor(popDat[[calc]][[question$groupVar]], levels = grpNames)
    }

    popDat
}


#' Determine if a particular analysis table is included in results list AND
#' that the table is not empty
#'
#' @param res            list containing multiple tables of data named by their analysis
#' @param analysisName   name of analysis for which results should be checked
#' @return logical; TRUE if results table exists and is not empty, FALSE otherwise
resGood <- function(res, analysisName){
    if(!analysisName %in% names(res)){
        warning(paste0("No results returned for analysis type: ",analysisName))
        return(FALSE)
    } else if(nrow(res[[analysisName]]) == 0){
        warning(paste0("Empty table return for analysis type: ",analysisName))
        return(FALSE)
    }
    return(TRUE)
}

#' Filter neighborhood log odds report
#' 
#' Filter neighborhood log odds report
#' 
#' @param res   results table of neighborhood log odds statistics
#' @param filt  filter in list form including values for min.estimate, max.fdr and min.cell.count
#' @return filtered results table
filterNeighborhoodLOreport <- function(res, filt){
    res %>%
    filter(`Neighborhood Population B total cell count` > filt$min.cell.count,
           (`Odds Ratio` > filt$min.estimate | `Inverse Odds Ratio` > filt$min.estimate),
           `adjusted p.value` < filt$max.fdr) %>%
    arrange(desc(abs(`log(Odds Ratio)`)), `adjusted p.value`)
}

#' Filter log odds report
#' 
#' Filter log odds report
#' 
#' @param res   results table of log odds statistics
#' @param filt  filter in list form including values for min.estimate, max.fdr and min.cell.count
#' @return filtered results table
filterLOreport <- function(res, filt){
    res %>%
    filter(`Population total cell count` > filt$min.cell.count,
           (`Odds Ratio` > filt$min.estimate | `Inverse Odds Ratio` > filt$min.estimate),
           `Fraction adjusted p.value` < filt$max.fdr) %>%
    arrange(desc(abs(`log(Odds Ratio)`)), `Fraction adjusted p.value`)
}

#' Filter density statistics
#' 
#' Filter density statistics
#' 
#' @param res   results table of density fold change and statistics
#' @param filt  filter in list form including values for min.estimate, max.fdr and min.cell.count
#' @return filtered results table
filterDensityReport <- function(res, filt){
    res %>%
    filter(`Cell State total cell count` > filt$min.cell.count,
           (`Fold Change` > filt$min.estimate | `Inverse Fold Change` > filt$min.estimate),
           `Density adjusted p.value` < filt$max.fdr,
           !duplicated(`Cell State`)) %>%
    arrange(desc(abs(`log(Fold Change)`)), `Density adjusted p.value`)
}

#' Filter neighborhood average counts statistics
#' 
#' Filter neighborhood average counts statistics
#' 
#' @param res   results table of neighborhood average counts statistics
#' @param filt  filter in list form including values for min.estimate, max.fdr and min.cell.count
#' @return filtered results table
filterNavgCountReport <- function(res, filt){
    res %>%
    filter(`Total Neighborhood Cell Count` > filt$min.cell.count,
           (`Fold Change` > filt$min.estimate | `Inverse Fold Change` > filt$min.estimate),
           `adjusted p.value` < filt$max.fdr) %>%
    arrange(desc(abs(`log(Fold Change)`)), `adjusted p.value`)
}


#' Convert list of results from multiple wilcox tests to a table
#' 
#' Results from multi.tests() are returned in list form, each element containing
#' raw test results for one condition. When test is wilcox.test, convert that
#' list to a table with columns CondTitle, estimate, CI.low, CI.high and p.value.
#' 
#' @param statsRes   results of multiple wilcox.test runs in list form
#' @param log2exp    logical; if tests were run on log values, set this value to TRUE to
#'                   convert estimate and confidence interval values back to normal scale
#' @return table of all wilcox test results, one row per test 
wilcoxResToTable <- function(statsRes, log2exp=FALSE){
    res <- tibble()
    for(cond in names(statsRes)){
        sr <- statsRes[[cond]]
        row <- tibble(CondTitle = cond,
                      estimate = sr$estimate,
                      CI.low = sr$`conf.int`[1],
                      CI.high = sr$`conf.int`[2],
                      p.value = sr$`p.value`)
        if(log2exp){
            row <- row %>%
                   mutate(estimate = exp(estimate),
                          CI.low = exp(CI.low),
                          CI.high = exp(CI.high))
        }
        res <- res %>% bind_rows(row)
    }
    res
}


#' Check for missing or extra samples in question groups
#' 
#' Loop through all tables in question data list and check that samples
#' assigned to each group are correct. If samples are missing, print a warning.
#' If there are extra samples in either group, throw an error.
#' 
#' @param qDat      list of tables, each containing all data for two groups to 
#'                  be compared; minimum columns include the calculation unit and
#'                  the groupVariable
#' @param sGrps     list of two vectors, each named by a group name and containing the
#'                  calculation units (default: FOV_ID) included in the group
#' @param groupVar  name of column containing group assignments (corresponding to names
#'                  of sGrps)
#' @param calcUnit  data column name containing values for which a single calculation should 
#'                    be made (default: FOV_ID)
#'
#' @return nothing 
qcQuestionGroups <- function(qDat, sGrps, groupVar, calcUnit="FOV_ID"){
    for(x in names(qDat)){
        log_debug(paste0("Checking samples in [",x,"]"))
        for(grp in names(sGrps)){
            smps <- qDat[[x]] %>% 
                    filter_at(all_of(groupVar), all_vars(. == grp)) %>%
                    pull(calcUnit) %>% sort() %>% unique()
            if(!identical(sort(smps), sort(sGrps[[grp]]))){
                notInData <- setdiff(sGrps[[grp]], smps)
                notInGroup <- setdiff(smps, sGrps[[grp]])
                log_warn(paste0("Discrepancies in ", x, " table in group ", grp, ": "))
                if(length(notInData) > 0 ){
                    log_warn(paste0("    Samples in group sample list but NOT in filtered data :",
                             paste(notInData, collapse=", ")))
                }
                if(length(notInGroup) > 0){
                    log_error(paste0("    Samples in filtered data but NOT group sample list: ", 
                              paste(notInGroup, collapse=", ")))
                    stop(paste0("    Samples in filtered data but NOT group sample list: ",
                              paste(notInGroup, collapse=", ")))
                }
            }
        }
    }
}

#' Parse table of conditions, separating single comma-delimited condition
#' strings into individual columns 
#' 
#' Parse single comma-delimited condition strings into the appropriate columns
#' such as cell type, subtype, functional, proliferation, etc.
#' 
#' @param condTbl     table containing comma-delimited condition strings to be parsed
#' @param condCol     name of column containing condition strings
#' @param condIdx     table of pre-indexed conditions corresponding to those in condTbl
#' @param cellTypes   table of cell type definitions
#' @param funcMarkers vector of all functional markers
#' @param funcCombos  vector of comma-delimited functional marker combinations
#'
#' @return parsed, normalized table
normalizeConditionTable <- function(condTbl, condCol, condIdx, cellTypes, 
                                    funcMarkers=NULL, funcCombos=NULL){
    tbl <- tibble()
    for(x in 1:nrow(condTbl)){
        row <- condTbl[x,] %>%
               cbind( normalizeCondition(condTbl[x,] %>% pull(condCol),
                                         cellTypes, funcMarkers, funcCombos)
               )
        tbl <- tbl %>% bind_rows(row)
    }
    tbl %>% 
    full_join(condIdx, by = intersect(names(tbl), names(condIdx))) %>% 
    select(ID, dplyr::matches("opulation"), everything())
}



#' Print QC metrics of statistics report
#' 
#' Print results counts for each tab of the report in addition
#' to any missing or extra conditions
#' 
#' @param qRes         list of tables each containing a full statistics report
#'                     for a single calculation type
#' @param allAnalyses  list of tables each containing the conditions for which a report
#'                     should be generated
#' @param samples      vector of samples that should be included in results
#' @param calcUnit     data column name containing values for which a single count should be made
#'                       (default: FOV_ID)
#'
#' @return nothing 
logStatsQC <- function(qRes, allAnalyses, samples, calcUnit = "FOV_ID"){

    resTabs2origTabs <- list(all_fractions_and_densities = c('fractions', 'densities'),
                             fractions_complete = 'fractions',
                             densities_complete = 'densities',
                             nbhd_fracs_complete = 'nfracs',
                             nbhd_averages_complete = 'navgcounts')

    qc <- tibble()
    for(tab in names(qRes)){
        tabQC    <- tibble()
        numRes   <- ifelse(!is.null(qRes[[tab]]), nrow(qRes[[tab]]), 0)
        missing  <- extra <- NA
        nmissing <- nextra <- 0

        if(tab %in% names(resTabs2origTabs)){
            for(origTab in resTabs2origTabs[[tab]]){
                if(!origTab %in% names(allAnalyses)){
                    tabQC <- tabQC %>% bind_rows(tibble(Sheet = origTab))
                    next
                }
                allConds <- allAnalyses[[origTab]]
                ## check that number of rows in results matches up with rows in 
                ## conditions input sheet
                expected <- nrow(allConds)
                if(numRes < expected){
                    missing <- allConds$`Cell State ID`[!allConds$`Cell State ID` %in% qRes[[tab]]$`Cell State ID`]
                    missing <- paste(unique(sort(as.numeric(missing))), collapse=",")
                    nmissing <- expected - numRes 
                } else if(numRes > expected){
                    extra <- qRes[[tab]]$`Cell State ID`[!qRes[[tab]]$`Cell State ID` %in% allConds$`Cell State ID`]
                    extra <- paste0(unique(extra), collapse = ",")
                    nextra <- numRes - expected
                }
                tabQC <- tabQC %>% 
                         bind_rows(tibble(Sheet = origTab,
                                          `Number rows in table` = numRes,
                                          `Number missing` = nmissing,
                                          `Number extra` = nextra,
                                          `Missing` = missing,
                                          `Extra` = extra))
            }
        } else if(tab == "samples"){
            ## check that sample list is as expected 
            numSamples <- length(samples)
            if(numRes < numSamples){
                nmissing <- numSamples - numRes
                missing  <- samples[ !samples %in% qRes[[tab]][[calcUnit]] ]
            } else if(numSamples < numRes){
                nextra <- numRes - numSamples
                extra  <- qRes[[tab]][[calcUnit]][ !qRes[[tab]][[calcUnit]] %in% samples ]
            }
            tabQC <- tibble(Sheet = tab,
                            `Number rows in table` = numRes,
                            `Number missing` = nmissing,
                            `Number extra` = nextra,
                            `Missing` = missing,
                            `Extra` = extra)
        } else {
            tabQC <- tibble(Sheet = tab,
                            `Number rows in table` = numRes)
        }

        qc <- qc %>% bind_rows(tabQC)
    }

    ## print QC info
    log_info("---------------------------------------------------------")
    longestName <- max(nchar(qc$Sheet))
    qc$Sheet <- sapply(qc$Sheet, function(x){ 
                     sp <- paste(rep(" ",longestName - nchar(x)), collapse="")
                     paste0(x, sp) 
                 }) %>%
                unlist()

    sp <- paste(rep(" ", longestName - nchar(names(qc)[1])), collapse="")
    names(qc)[1] <- paste0(names(qc)[1], sp) 
    log_info(paste(names(qc), collapse="\t"))
    for(n in 1:nrow(qc)){
        log_info(paste(unlist(qc[n,]), collapse="\t"))
    }
    log_info("---------------------------------------------------------")

}



getTMEstats <- function(fracDat, groupVar, calcUnit, fracCol = "Fraction"){

    dat <- fracDat %>%
           mutate(CondTitle = nbhdFractionCond(`Center Population A`, `Neighborhood Population A`,
                                               `Center Population B`, `Neighborhood Population B`))

    loTbl <- dat %>%
             fractionsToLO(fracCol, groupVar) %>%
             select_at(c("CondTitle", calcUnit, groupVar, "LO")) %>%
             spread(CondTitle, LO)

    ## pull out IDs and condition parts
    conds <- dat %>% 
             ungroup() %>%
             select(CondTitle, `Cell State ID`,
                    `Center Population A`, `Neighborhood Population A`,
                    `Center Population B`, `Neighborhood Population B`) %>%
             unique()

    stats  <- loTbl %>% 
              getStatsLO(groupVar, unique(dat$CondTitle)) %>%
              select(-`adjusted p.value`)

    counts <- tibble(CondTitle = unique(dat$CondTitle),
                     `Center cell type total count` = nrow(dat),
                     `Neighborhood Population A total cell count` = sum(dat$`N.pop.A_Count`),
                     `Neighborhood Population B total cell count` = sum(dat$`N.pop.B_Count`))

    meds   <- dat %>% getGroupMedianFractions(groupVar, fracCol)

    conds %>% 
    left_join(counts, by = "CondTitle") %>% 
    left_join(meds, by = "CondTitle") %>% 
    left_join(stats, by = "CondTitle")

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

    tmeStats <- lapply(1:nrow(analyses), function(x){

                    pop <- unlist(analyses[x,])
                    log_debug(paste0("[", pop[3], "]    ", pop[1], "/", pop[2]))
                    csid <- pop[3]

                    fracs <- getTMEfractions(tme, csid, "Tumor", pop[1], pop[2])
                    stats <- reportTMELogOdds(fracs, groupVar, calcUnit)

                    stats
                 })

    tmeStats %>% bind_rows()
}


#' Generate full statistics report for all analyses within a single study question
#' 
#' A single study question seeks to find interesting differences between two groups
#' of samples (or FOVs, etc.). For each calculation type falling under both general and 
#' spatial analyses, run statistics and report results. Return all results in list form
#' where each element in the list is results for a single calculation or supporting question
#' information like a list of samples included in each comparison group, filters that
#' were used, etc.
#' 
#' @param question       a list describing how to divide samples for one study question
#'                       (see docs for XLSX file format, to be parsed with parseStudyQuestions())
#' @param dat            table of all annotated cell data
#' @param sampAnn        all sample annotation necessary for assigning samples to a group
#' @param allAnalyses    parsed and formatted list describing all analyses to be run for 
#'                       each calculation type (see docs for XLSX file format)
#' @param markers        vector of all markers used to distinguish markers from other classifiers
#'                       in condition strings
#' @param metricsDir     directory containing all pre-computed metrics to be loaded; this directory 
#'                       is required but may be empty; if that is the case, metrics will be calculated
#'                       for only samples involved in current question
#' @param filtList       list containing values for filtering; must include: max.fdr, min.estimate, 
#'                       and min.cell.count
#' @param calcUnit       data column name containing values for which a single count should be made
#'                       (default: FOV_ID)
#' @param nbhdCounts     table of formatted neighborhood counts to be used for spatial analyses; may be
#'                       NULL if allAnalyses does not contain any spatial analyses
#' @param tumorNbhdCells vector of UUIDs for cells that are inside the neighborhood of at least one 
#'                       tumor cell 
compareSampleGroups <- function(question, dat, sampAnn, allAnalyses, markers, metricsDir, 
                                filtList, filtName = "filtered", calcUnit = "FOV_ID", 
                                nbhdCounts = NULL, tumorNbhdCells = NULL){

    toRename <- c("CI.low", "CI.high", "p.value", "adjusted p.value") ## results columns to be renamed

    ###
    ### Prep question data
    ###
    sGrps    <- getQuestionGroups(question, dat, calcUnit = calcUnit, tumorNbhdCells = tumorNbhdCells)
    allSmps  <- unique(unlist(sGrps, use.names=F))
    groupVar <- question$groupVar

    log_info(question$question)

    qDat <- suppressMessages(getQuestionData(question, sGrps, dat, allAnalyses, names(allAnalyses),
                                             markers, metricsDir, calcUnit = calcUnit,
                                             nbhdCounts = nbhdCounts, tumorNbhdCells = tumorNbhdCells))
    if(is.null(qDat)){
        log_warn("One or more groups in this question contains no data. Skipping.")
        return(NULL)
    }
    qcQuestionGroups(qDat, sGrps, question$groupVar)

    qRes <- list()

    ###
    ### Fractions/log odds
    ###
    if("fractions" %in% names(allAnalyses) && resGood(qDat, "fractions")){
        log_info("Running fractions stats")
        fltr <- filtList$fractions
        lo   <- reportLogOdds(qDat$fractions, groupVar, calcUnit = calcUnit,
                              max.fdr = fltr$max.fdr, min.or = fltr$min.estimate)
        lof  <- lo %>% 
                filter(`Fraction Passed` == "X") %>% 
                select(-`Fraction Passed`) %>%
                arrange(desc(abs(`log(Odds Ratio)`)), `adjusted p.value`)

        qRes$fractions_complete <- lo
        qRes[[paste0("fractions_", filtName)]] <- lof
    } else {
        log_info("No fractions stats to report.")
    } 

    ###
    ### Densities
    ###
    if("densities" %in% names(allAnalyses) && resGood(qDat, "densities")){
        log_info("Running density stats")
        fltr    <- filtList$densities
        den     <- reportDensityStats(qDat$densities, allSmps, groupVar, calcUnit = calcUnit) 
        denf    <- den %>% 
                   filter(`Density Passed` == "X") %>% 
                   select(-`Density Passed`) %>%
                   arrange(desc(abs(`log(Fold Change)`)), `adjusted p.value`)

        qRes$densities_complete <- den
        qRes[[paste0("densities_", filtName)]] <- denf
    } else {
        log_debug("No densities stats to report")
    }

    ###
    ### Macro Neighborhood Fractions
    ###
    if("nfracs" %in% names(allAnalyses) && resGood(qDat, "nfracs")){
        log_info("Running neighborhood fraction stats")
        nlo     <- reportNeighborhoodLogOdds(qDat$nfracs, qDat$ncounts, allSmps, 
                                             groupVar, calcUnit = calcUnit, fracCol = "MeanFraction")
        nlof    <- filterNeighborhoodLOreport(nlo, filtList$nFracs) %>%
                   arrange(desc(abs(`log(Odds Ratio)`)), `adjusted p.value`)

        qRes$nbhd_fracs_complete <- nlo
        qRes[[paste0("nbhd_fracs_", filtName)]] <- nlof
    } else {
        log_debug("No neighborhood fractions stats to report.")
    }

    ###
    ### Macro Neighborhood Average Counts
    ###
    if("navgcounts" %in% names(allAnalyses) && resGood(qDat, "navgcounts")){
        log_info("Running neighborhood average count stats")
        navgs     <- reportAverageCountsStats(qDat$navgcounts %>% 
                                              select_at(c("Cell State ID","CenterCellType", "NeighborhoodCellType",
                                              calcUnit, "AvgNbhdCellTypeCount", groupVar)) %>% 
                                              unique(), 
                                              allSmps, groupVar, qDat$ncounts)
        navgsf    <- filterNavgCountReport(navgs, filtList$nAvgCounts) %>%
                     arrange(desc(abs(`log(Fold Change)`)), `adjusted p.value`)

        qRes$nbhd_averages_complete <- navgs
        qRes[[paste0("nbhd_averages_", filtName)]] <- navgsf
    } else {
        log_debug("No neighborhood average counts stats to report.")
    }

    ###
    ### Combine full fraction and density results into one 'mega sheet'
    ###
    allFracDen <- tibble()
    lo2 <- tibble()
    if(!is.null(qRes$fractions_complete)){ 
        idxs <- which(names(lo) %in% toRename)
        lo2 <- lo
        names(lo2)[idxs] <- paste("Fraction", names(lo2)[idxs])
        allFracDen <- lo2 %>% 
                      select(`Fraction Passed`, 
                             `Cell State ID`, 
                             `Cell State`, 
                             Population, 
                             everything()) %>%
                      arrange(as.numeric(`Cell State ID`))
        qRes$fractions_complete <- NULL
    }
    if(!is.null(qRes$densities_complete)){
        idxs <- which(names(den) %in% toRename)
        den2 <- den
        names(den2)[idxs] <- paste("Density", names(den2)[idxs])

        allFracDen <- allFracDen %>%
                      full_join(den2, by = intersect(names(lo2), names(den2))) %>% 
                      select(any_of(c("Fraction Passed", 
                                      "Density Passed",
                                    "Cell State ID", 
                                    "Cell State", 
                                    "Population",
                                    names(.)))) %>%
                      arrange(as.numeric(`Cell State ID`))
        qRes$densities_complete <- NULL
    }
    qRes$all_fractions_and_densities <- allFracDen

    log_info("compiling final report")

    ###
    ### log filter criteria
    ###
    log_debug("  reporting filters used")
    filtTbl <- tibble()
    for(f in names(filtList)){
        filtTbl <- filtTbl %>%
                   bind_rows(filtList[[f]] %>% as_tibble() %>% mutate(Calculation = f))
    }
    filtTbl <- filtTbl %>% gather(1:3, key='filter', value='value') %>% spread(Calculation, value, fill=0)

    qRes[[filtName]] <- filtTbl

    ###
    ### summarize samples
    ###
    log_debug("  reporting summary of samples included in comparison groups")
    samples        <- sampAnn %>% filter(!!as.name(calcUnit) %in% allSmps) 
    sample_summary <- samples

    if(groupVar %in% names(samples)){
        sample_summary <- sample_summary %>%  
                          group_by_at(question$groupVar) %>% 
                          summarize(NumberSamples = n())
    }

    qRes$samples <- samples
    qRes$sample_summary <- sample_summary

    sortedTabs <- c("all_fractions_and_densities",
                    "fractions_complete", 
                    paste0("fractions_",filtName),
                    "densities_complete",
                    paste0("densities_",filtName),
                    "nbhd_fracs_complete",
                    paste0("nbhd_fracs_",filtName),
                    "nbhd_averages_complete",
                    paste0("nbhd_averages_",filtName),
                    filtName,
                    "samples",
                    "sample_summary")

    log_debug("  QCing results counts")
    logStatsQC(qRes, allAnalyses, allSmps, calcUnit = calcUnit)

    qRes[sortedTabs[sortedTabs %in% names(qRes)]]
}


#' Filter statistics results for biologically significant results
#' 
#' Remove from statistics results any conditions that do not pass
#' criteria including minimum population and subpopulation counts,
#' minimum difference in group medians and FDR.
#'
#' @param dat               tibble of statistics to be filtered
#' @param calc              calulation type ['fractions'|'densities'] (default: fractions)
#' @param min.subpop.count  minimum subpopulation cell count (default: 300)
#' @param min.pop.count     minimum population cell count (default: 1000)
#' @param min.median.diff   minimum difference between group median fraction or 
#'                          density (default = 0.1 or 10%)
#' @param pop.col           column header of population cell count 
#'                          (default = "Population total cell count")
#' @param subpop.col        column header of subpopulation cell count
#'                          (default = "Cell State total cell count")
#' @param fdr.col           column header of FDR column (default = "adjusted p.value")
#' @param median.cols       vector of column headers of median values
#' 
#' @return tibble of conditions that pass all criteria
bioFilter <- function(dat, calc = "fractions", 
                      min.subpop.count = 300, min.pop.count = 1000,
                      min.median.diff = 0.1, max.fdr = 0.05,
                      subpop.col = 'Cell State total cell count',
                      pop.col = 'Population total cell count',
                      fdr.col = 'adjusted p.val', median.cols = NULL){

    dat <- dat %>%
           filter(!!as.name(subpop.col) >= min.subpop.count,
                  !!as.name(pop.col) >= min.pop.count,
                  !!as.name(fdr.col) < max.fdr)
    if(calc == "fractions"){
         return(dat %>%
                filter(abs(!!as.name(median.cols[2]) - !!as.name(median.cols[1])) >= min.median.diff))
    } else if(calc == "densities"){
         return(dat %>%
                filter(abs((!!as.name(median.cols[2]) - !!as.name(median.cols[1]))/!!as.name(median.cols[2])) >= min.median.diff))
    }
    msg <- paste0("Unrecognized calc: [", calc, "]")
    log_error(msg)
    stop(msg)
}

#' Calculate harmonic mean p-value
#' 
#' Given a vector of p values (pre-adjustment), calculate the
#' harmonic mean p-value
#' 
#' @param pvals   numeric vector of p-values
#' @param na.rm   logical indicating whether NA values should be removed
#'                prior to HMP calculation
#' 
#' @return a single numeric value representing harmonic mean p-value
HMP <- function(pvals, na.rm = FALSE){
    if(na.rm){
        pvals <- pvals[!is.na(pvals)]
    }
    1/sum(1/pvals)
}

#' Calculate harmonic mean p-value and FDR for multiple iterations of
#' the same statistical test, each with different groups of samples
#' 
#' Used, in our case, to test consistency of patient level comparisons,
#' use p-values from patient level comparisons of a certain condition to
#' calculate a harmonic mean p-value and FDR
#' 
#' @param  dat           tibble including all patient level results of all conditions 
#'                       in question (conditions distinguished by Cell State ID)
#' @param  indivPvalCol  column name of patient level p-values (pre-adjustment)
#' @param  na.rm         logical indicating whether to remove NA results prior to 
#'                       HMP calculation
#' 
#' @return tibble from dat with columns for HMP and FDR added
getHMP <- function(dat, indivPvalCol, na.rm = FALSE){
    dat %>%
    mutate(!!as.name(indivPvalCol) := ifelse(is.na(!!as.name(indivPvalCol)), 1, 
                                            !!as.name(indivPvalCol)),
           `p.val one side` = ifelse(Value != -1,
                                     !!as.name(indivPvalCol)/2,
                                     1 - !!as.name(indivPvalCol)/2)) %>%
    group_by(`Cell State ID`, `Cell State`, Population) %>%
    mutate(HMP = HMP(`p.val one side`, na.rm = na.rm),
           FDR = p.adjust(HMP, method = "bonferroni")) %>%
    select(`Cell State ID`, `Cell State`, Population, HMP, FDR, everything()) %>%
    unique()
}


#' Get a vector of Cell State IDs for which zero patient-level stats 
#' are determined significant
#'
#' @param  dat         tibble of patient level stats
#' @param  indivIDCol  column containing IDs of patient level comparisons
#' @param  signifCol   column indicating significance; only rows with NA in this column
#'                     are condidered to NOT be significant
#' 
#' @return vector of Cell State IDs for which no patient level stats
#'         were determined to be significant
condsNoSignif <- function(dat, indivIDCol, signifCol){
    dat %>%
    select_at(c("Cell State ID", indivIDCol, signifCol)) %>%
    unique() %>%
    spread(indivIDCol, signifCol) %>%
    filter_at(all_of(unique(dat[[indivIDCol]])), all_vars(is.na(.))) %>%
    pull(`Cell State ID`)
}



