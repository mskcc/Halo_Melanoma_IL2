#' Filter annotated cell data for a single condition
#' 
#' Get all cells that meet all criteria outlined in condition string
#' 
#' @param dat           table of annotated cell data, where one row represents a single cell
#'                      and, at minimum, columns for Classifiers and PositiveMarkers 
#' @param conditionStr  comma-delimited string containing a combination of cell types (classifiers) 
#'                      and/or markers
#' @param markers       vector of individual markers used to distinguish markers from classifiers
#'
#' @return subset of dat including only cells that match ALL criteria included in conditionStr
filterDataForCondition <- function(dat, conditionStr, markers){

    if(conditionStr == "DAPI"){ return(dat) }

    pts <- unlist(strsplit(conditionStr,","))
    if("DAPI" %in% pts){ pts <- pts[!pts == "DAPI"] }

    negMarkers <- gsub("-","",pts[grep("-",pts)])
    posMarkers <- pts[which(pts %in% markers)]

    if(length(c(negMarkers, posMarkers)) > 0){
        cats <- pts[-which(gsub("-","",pts) %in% c(negMarkers, posMarkers))]
    } else {
        cats <- pts
    }

    ## filter for cells that satisfy ALL categories
    if(!is.null(cats) && length(cats) > 0){
        for(cat in cats){
            dat <- dat %>%
                   filter(grepl(getClassifierPattern(cat), Classifiers, perl=T))
        }
    }

    ## filter for cells positive for ALL positive markers
    if(!is.null(posMarkers) && length(posMarkers) > 0){
        if("DAPI" %in% posMarkers){ posMarkers <- posMarkers[!posMarkers == "DAPI"] }
        if(!is.null(posMarkers) && length(posMarkers) > 0){
            pats <- unlist(sapply(posMarkers, function(x){ paste(getClassifierPattern(x,delim=",")) }))
            for(pat in pats){
                dat <- dat %>%
                       filter(grepl(pat, PositiveMarkers, perl=TRUE))
            }
        }
    }
    ## filter for cells negative for ALL negative markers
    if(!is.null(negMarkers) && length(negMarkers) > 0){
        pats <- unlist(sapply(negMarkers, function(x){ getClassifierPattern(x, delim=",") }))
        dat <- dat %>%
               filter(!grepl(paste(pats, collapse="|"), PositiveMarkers, perl=TRUE))
    }

    return(dat)
}

#' Filter data tibble for a selection of calculation units
#'
#' Filter data tibble for a selection of calculation units
#' 
#' @param dat       data tibble to be filtered
#' @param calcUnit  name of column corresponding to the calculation unit of 
#'                  interest (default: FOV_ID)
#' @param include   vector of calculation units to keep
#'
#' @return subset of dat consisting of only calcuation units in 'include' vector
filterForCalculationUnits <- function(dat, calcUnit = "FOV_ID", include = NULL){
    if(!is.null(include)){
       if(!all(include %in% dat[[calcUnit]])){
           wrn <- paste0("Calculation units missing from data: ",
                           paste(include[!include %in% dat[[calcUnit]] ], collapse=", "))
           log_warn(wrn)
           warning(wrn)
       }
       return( dat %>% filter_at(vars(calcUnit), all_vars(. %in% include)) )
    }
    dat
}


#' Filter data for cells that are part of a single comparison group
#'
#' Given a list that defines a single group of cells, filter data for only those cells
#'
#' @param dat             data table to be filtered
#' @param group           list describing data to be included in group; list names are
#'                        either column names from dat or 'Cell Region' and values are
#'                        the values to keep
#' @param tumorNbhdCells  vector of UUIDs for cells that fall within the neighborhood of at least
#'                        one tumor cell
#'
#' @return filtered data table
filterForComparisonGroup <- function(dat, group, tumorNbhdCells = NULL){
    grp <- dat
    for(filt in names(group)){
        if(filt == "Cell Region"){
            grp <- grp %>%
                   filterDataForCellRegion(tolower(group[[filt]]), tumorNbhdCells = tumorNbhdCells)
        } else if(filt %in% names(dat)){
            grp <- grp %>%
                   filter_at(vars(filt), all_vars(. %in% group[[filt]]))
        } else {
            warning(paste0("Can not filter on column that does not exist: ",filt,"...skipping."))
        }
    }
    grp
}

#' Filter for cells that are in a specific cell region
#' 
#' Remove from data any cells that fall outside the region of interest
#' 
#' @param dat             any halo data table; in order to filter for cells within
#'                        the interface, table must have a 'Band' column
#' @param cellRegion      the region to filter for
#'                          - fov: region is total FOV, return entire dat table
#'                          - interface: keep any row that is assigned to a Band
#'                          - interface inside: keep any row that is in a negative distance Band
#'                          - interface outside: keep any row that is in a positive distance Band
#' @param tumorNbhdCells  character vector of cell UUIDs for cells that are in the neighborhood of
#'                        at least one tumor cell; to be used for filtering data for 'neighborhood'
#'                        cell region
#'
#' @return  filtered data table
filterDataForCellRegion <- function(dat, cellRegion, tumorNbhdCells = NULL){
    valid <- c("fov","interface","interface inside", "interface outside", "neighborhood")
    if(!tolower(cellRegion) %in% valid){
        stop(paste0("Unrecognized cell region: ",cellRegion))
    }
    if(tolower(cellRegion) == "neighborhood" && is.null(tumorNbhdCells)){
        stop(paste0("Attempted to filter for neighborhood but neighborhood UUIDs not provided."))
    }
    dat <- switch(tolower(cellRegion),
                  "interface" = dat %>% filter(!is.na(Band)),
                  "interface inside" = dat %>% filter(!is.na(Band), grepl("\\-", Band)),
                  "interface outside" = dat %>% filter(!is.na(Band), !grepl("\\-", Band)),
                  "neighborhood" = dat %>% filter(UUID %in% tumorNbhdCells),
                  dat)
    dat
}


#' Extract area information from annotated cell data for each FOV
#' for a specific cell region
#'
#' Pull out only area information from annotated cell data for each FOV
#' depending on the cell region of interest
#' 
#' @param dat           annotated cell data that has been joined with cell region areas
#'                      for cell region of interest
#' @param cellRegion    region of interest [fov|interface|interface inside|interface outside]
#' @param groupBy       data column name(s) containing values for which a single area value
#'                      should be calculated (default: FOV_ID)
#'
#' @return  a table containing columns [groupBy], Area 
filterForAreas <- function(dat, cellRegion = "fov", groupBy = "FOV_ID"){
    groupCols  <- unique(c("Sample_ID","FOV_ID","Band",groupBy)) ## cols required to identify unique values;
                                                           ## must keep them so that two rows
                                                           ## that happen to be identical without any
                                                           ## one of them are all kept (e.g., a sample 
                                                           ## in which Band1 in FOV1 has the same area
                                                           ## as Band1 in FOV2
    areaCol    <- "FOVBandArea"

    if(tolower(cellRegion) == "fov"){
        groupCols <- groupBy
        areaCol <- "FOVArea"
    }

    dat %>%
    select_at(c(groupCols, areaCol)) %>%
    unique() %>%
    rename_at(vars(areaCol), list(~ (. = "Area"))) %>%
    filterDataForCellRegion(cellRegion) %>%
    group_by_at(groupBy) %>%
    summarize(Area = sum(Area, na.rm = T))
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

#' Filter table of neighborhood data for neighborhood cells
#' that satisfy all filter criteria including classifiers and/or positive markers
#' 
#' @param dat   data tibble of neighborhood data including columns for Center cells
#'              and columns for Neighborhood cells
#' @param nbhd  comma-delimited character string describing neighborhood cells to keep
#'              (e.g., "Tconv8,PD1,LAG3,TIM3-")
#' 
#' @return  filtered data tibble
filterForNeighborhood <- function(dat, nbhd){
    n <- unlist(strsplit(nbhd,","))
    classes <- n[!n %in% c(markers, paste0(markers, "-"))]
    mrkrs <- n[n %in% c(markers, paste0(markers,"-"))]

    if(length(classes) > 0){
        for(cls in classes){
            if(cls == "EXM0"){
                dat <- dat %>%
                       filter(!grepl("EXM",N.Classifiers))
            } else {
                dat <- dat %>%
                       filter(grepl(getClassifierPattern(cls,delim=";"),N.Classifiers))
            }
        }
    }
    if(length(mrkrs) > 0){
        for(mrk in mrkrs[!grepl("\\-", mrkrs)]){ dat <- dat %>%
                                                  filter(grepl(getClassifierPattern(mrk,delim=","),
                                                             N.FullPosMarkerStr)) }
        for(mrk in mrkrs[grepl("\\-", mrkrs)]){ dat <- dat %>%
                                                filter(!grepl(getClassifierPattern(gsub("-","",mrk), delim=","),
                                                              N.FullPosMarkerStr)) }
    }
    return(dat)
}


