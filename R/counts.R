#' Summarize counts of cell type marker combinations and their classifiers (cell types)
#' 
#' Count cell type marker combinations in a number of ways, by Sample and FOV for example. 
#' Output will include classification based on Classifiers column in dat and cell type class columns 
#' specified to be output explicitly in their own columns (e.g., Cell_type, Tag)
#' 
#' 
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param cellTypes       tibble of cell type assignments
#' @param classCols       column names from cellTypes tibble indicating which categories should be 
#'                        summarized; default=c("Category","Cell_type","Subtype","Tag")
#' @param summarizeBy     vector of ways in which data should be grouped; values must be equal to
#'                        column names in dat; default=c("Sample","FOV")
#' @param xlsxFile        name of XLSX file to which counts summaries should be written; if NULL, no
#'                        XLSX file will be written
#' @return list of summary tables
#' @export
getMarkerComboCounts <- function(dat, cellTypes, classCols=NULL, summarizeBy=NULL, xlsxFile=NULL){

    if(is.null(classCols)){ classCols <- c("Category", "Cell_type", "Subtype", "Tag") }
    if(is.null(summarizeBy)){ summarizeBy <- c("Sample_ID", "FOV_ID") }

    ctClasses <- cellTypes %>%
                 filter(Category != "Functional") %>%
                 select_at(classCols) %>%
                 unlist() %>%
                 unique()
    ctClasses <- ctClasses[!is.na(ctClasses)]

    ## remove any non-cell type classifiers from Classifiers column
    classifiers <- lapply(dat$Classifiers, function(x){ unlist(strsplit(x,";")) })
    dat$ctClassifiers <- unlist(lapply(classifiers, function(x){ paste(x[x %in% ctClasses],collapse=";") }))
    dat$IDMarkers[which(dat$IDMarkers == "")] <- "superNeg"
    dat$ctClassifiers[which(dat$IDMarkers == "superNeg")] <- "superNeg"
    dat$ctClassifiers[which(dat$ctClassifiers == "")] <- "Unreasonable"

    #### get all unique cell type combinations with counts and classifiers
    ctsTotal <- dat %>%
                group_by(IDMarkers,ctClassifiers) %>%
                summarize(Total = n()) %>%
                arrange(desc(Total))

    sheetLst <- list()

    for(sb in summarizeBy){
        ## get counts for each markerPosTag+ctClassifier combo by sample
        sheetLst[[paste0("combos_by_",sb)]] <- ctsTotal %>%
                                               left_join(dat %>%
                                                         group_by_at(c(sb, "IDMarkers", "ctClassifiers")) %>%
                                                         summarize(Count = n()) %>%
                                                         spread(sb, Count),
                                                     by = c("IDMarkers", "ctClassifiers"))
    }
    return(sheetLst)
}

#' Count cells in each class/cell type
#'
#' Summarize cell counts by cell types and data level (e.g., all, Sample, FOV)
#' 
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param cellTypes       table of cell type definitions
#' @param classCols       column names from cellTypes tibble indicating which categories should be 
#'                        summarized; default=c("Category","Cell_type","Subtype","Tag")
#' @param summarizeBy     vector of ways in which data should be grouped; values must be equal to
#'                        column names in dat; default=c("Sample","FOV")
#' @return list of tables, each containing a different summary of counts
#' @export
getClassCountSummaries <- function(dat, cellTypes, classCols=NULL, summarizeBy=NULL){

    if(is.null(classCols)){ classCols <- c("Category", "Cell_type", "Subtype", "Tag") }
    if(is.null(summarizeBy)){ summarizeBy <- c("Sample_ID", "FOV_ID") }

    tbls <- list()
    for(sb in summarizeBy){
        for(cl in classCols){
            ## summarize standard classes (those found in individual class columns in data)
            sheet <- summarizeCellCounts(dat, summarizeBy = c(cl, sb)) %>%
                     spread(sb,CellCount) %>%
                     ungroup() %>%
                     mutate(Total = rowSums(.[,-1], na.rm=T)) %>%
                     select(!!as.name(cl), Total, everything()) %>%
                     filter_at(vars(cl), all_vars(!is.na(.)))

            ## summarize additional misc classes that are only in Classifiers column
            addlVals <- cellTypes %>% 
                        filter_at(vars(cl), all_vars(!. %in% sheet[[cl]])) %>%
                        pull(cl) %>%
                        unique()
            sheet <- sheet %>%
                     bind_rows(lapply(addlVals, function(x){
                                  summarizeClassCounts(dat, x, cl, summarizeBy = sb) 
                              })
                     ) 
            tbls[[paste0(cl,"_by_",sb)]] <- sheet 
        }
    }

    return(tbls)
}

summarizeCellCounts <- function(dat, summarizeBy=NULL){
    
    checkGrouping(names(dat), summarizeBy)

    dat %>%
    group_by_at(summarizeBy) %>%
    summarize(CellCount = n()) 

}

summarizeClassCounts <- function(dat, class, classType, summarizeBy=NULL){

    checkGrouping(names(dat), summarizeBy)

    dat <- dat %>%
           getClassifierCounts(class, summarizeBy = summarizeBy) 

    if(!is.null(summarizeBy)){
        dat <- dat %>% spread(summarizeBy, Count) 
    }

    dat <- dat %>%
           mutate(Total = rowSums(.), Class = class) %>%
           select(Class, Total, everything()) %>%
           rename(!!classType := Class)
    dat
}

#' Get count(s) for a single marker
#' 
#' Get number of cells that are positive for a single marker
#' 
#' @param dat      table of cell level data where a row represents a single unique cell 
#'                 and at minimum, columns include any values in {groupBy} vector and PositiveMarkers
#' @param marker   marker to count
#' @param groupBy  vector of column names to group data by before counting; default=NULL; set to
#'                 NULL to get total count
#' @return tibble of count(s) with one row per group
#' @export
countMarker <- function(dat, marker, groupBy=NULL){

    checkGrouping(names(dat), groupBy)

    dat %>% 
    filter(grepl(getClassifierPattern(marker, delim=","), PositiveMarkers)) %>% 
    group_by_at(groupBy) %>% 
    summarize(Count=n()) %>% 
    dplyr::rename_at(vars(Count), list(~ (. = marker)))

}

#' Get counts for each individual marker
#'
#' Get number of cells that are positive for each individual marker
#' 
#' @param dat      table of cell level data where a row represents a single unique cell 
#'                 and at minimum, columns include any values in {groupBy} vector and PositiveMarkers
#' @param marker   vector of markers to count
#' @param groupBy  column name to group data by before counting; default=NULL; set to
#'                 NULL to get total counts
#' @return table containing a single row for each group and a column for each marker
#' @export
getAllIndividualMarkerCounts <- function(dat, markers, groupBy=NULL){

    allCounts <- lapply(markers, function(x){ countMarker(dat, x, groupBy=groupBy) }) %>%
                 bind_rows() %>%
                 gather(markers, key="Marker", value="Count") %>%
                 filter(!is.na(Count)) %>%
                 spread(Marker, Count) %>%
                 arrange_at(groupBy)
                 
    allCounts
}


#' Count and summarize cell types by looking at positive marker combinations
#' 
#' Determine the markers expressed in each cell and count and summarize the unique combinations of
#' expressed markers, cell types/subtypes, and cells that touch or do not touch macrophages
#'
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param cellTypes       tibble of cell type assignments
#' @param classCols       column names from cellTypes tibble indicating which categories should be 
#'                        summarized; default=c("Category","Cell_type","Subtype","Tag")
#' @param summarizeBy     vector of ways in which data should be grouped; values must be equal to
#'                        column names in dat; default=c("Sample","FOV")
#' @param xlsxFile        name of XLSX file to which counts summaries should be written; if NULL, no
#'                        XLSX file will be written
#' @param markerCombos    logical; when TRUE, all cell type marker combinations will be counted
#' @param classes         logical; when TRUE, counts will be summarized by classes
#' @param touchingMacros  logical; when TRUE, class counts will be repeated on the TM cells only 
#' @return list of tibbles, each a different summary of counts
#' @export 
getCountSummaries <- function(dat, cellTypes, classCols=NULL, summarizeBy=NULL, xlsxFile=NULL,
                                    markerCombos=TRUE, classes=TRUE, touchingMacros=TRUE){

    if(is.null(classCols)){ classCols <- c("Category", "Cell_type", "Subtype", "Tag") }
    if(is.null(summarizeBy)){ summarizeBy <- c("Sample_ID", "FOV_ID") }

    allSheets <- list()
   
    if(markerCombos){
        comboSheets <- getMarkerComboCounts(dat, cellTypes, classCols=classCols, summarizeBy=summarizeBy)
        allSheets <- c(allSheets, comboSheets)
    }
    if(classes){
        classSheets <- getClassCountSummaries(dat, cellTypes, classCols=classCols, summarizeBy=summarizeBy)
        allSheets <- c(allSheets, classSheets)
    
        if(touchingMacros){
            ## break TM category by Cell_type and by Subtype
            TMdat <- dat %>% filter(grepl(getClassifierPattern("TM"), Classifiers))
            TMclassSheets <- getClassCountSummaries(TMdat, cellTypes, classCols=c("Cell_type", "Subtype"), 
                                            summarizeBy=summarizeBy)
            names(TMclassSheets) <- paste0("TM_",names(TMclassSheets))
            allSheets <- c(allSheets, TMclassSheets)
        }
    }
    if(!is.null(xlsxFile)){
        openxlsx::write.xlsx(allSheets, file=xlsxFile)
    }

    return(allSheets)

}

#' Summarize counts for a single cell type included in Classifiers column
#' 
#' Get counts for a single cell type for entire data set, or summarized by one or
#' more columns
#'
#' @param dat         table where column headers include Classifiers, and each row
#'                    represents a single cell
#' @param cellType    character string to pull from Classifiers column
#' @param summarizeBy name of column in dat by which counts should be grouped; default
#'                    is NULL, so counts returned will be for complete data set
#' @return table of counts 
#' @export
getClassifierCounts <- function(dat, cellType, summarizeBy=NULL){

    checkGrouping(names(dat), summarizeBy)

    dat %>% 
    filter(grepl(getClassifierPattern(cellType), Classifiers)) %>%
    group_by_at(summarizeBy) %>%
    summarize(Count=n())
}


