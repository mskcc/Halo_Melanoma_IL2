#' Make sure a cell can be classified into only one of a pair or list of groups
#' 
#' Given a tibble with UUIDs and classification for each cell, check that
#' each cell fits into just one of the groups being analyzed
#checkForDoubleCounts <- function(cellClassification, groups, delim=";"){  #### ADD UNIT TEST

#    classes <- strsplit(cellClassification$Classifiers, delim)

#    for(g in groups){
#        cellClassification[[g]] <- 0
#        cellClassification[[g]][g %in% classes]
#    }
    
#}

#' Generate vector of positive marker combinations
#' 
#' Given a vector of markers, expand into vector of all combinations of those
#' markers, with the option of including the empty combination (no positive markers)
#' 
#' @param  markers     vector of individual markers
#' @param  into        character string; type of vector 'into' which markers should be expanded; currently function
#'                     only supports expanding into 'allCombinations'
#' @param  emptyCombo  logical; TRUE indicates the empty combo should be included in final vector; default=FALSE
#' @return  list of expanded cell types
#' @export
expandMarkerCombos <- function(markers, into="allCombinations", emptyCombo=FALSE){
    markers <- unlist(strsplit(markers,","))

    combos <- switch(tolower(into),
                     "allcombinations" = getAllCombos(markers),
                     "singlemarkers" = markers,
                     NULL)

    if(is.null(combos)){
        log_error("Code currently can only expand markers into 'allCombinations' or 'singleMarkers'")
        stop(paste0("Unsupported 'into' argument: ",into))    
    }

    if(emptyCombo){
        combos <- c(combos, "")
    }

    return(combos)
}

#' Assign a label to a group of marker combinations (cell types) 
#' 
#' Assign a label to a group of marker combinations (cell types) either by extracting
#' the positive markers in each combo (e.g., 'CD4,CD163') or by using a tag, plus the 
#' number of positive markers in a combo (e.g., EXM1, EXM2)
#' 
#' @param combos      vector of comma-delimited marker combinations to be labeled
#' @param labelUsing  either 'positive markers' or 'tag+number positives', indicating how
#'                    to label all combinations
#' @param tags        if labeling using 'tag+number positives', a vector of tags, one for each
#'                    combination in {combos} <--- DOUBLE CHECK THAT THAT IS CORRECT 
#' @return vector of labels
#' @export
labelCellTypes <- function(combos, labelUsing="positive markers", tags=NULL){

    validVals <- c("positivemarkers", "tag+numberpositives")
    labelUsing <- gsub(" ","",labelUsing)

    if(!labelUsing %in% validVals){
        log_warn(paste("WARNING: Invalid value for 'LabelCombinationsBy': ",
                    labelUsing,
                    ". Labeling using positive markers instead."))
        labelUsing <- "positivemarkers"
    }

    combos <- sapply(combos, function(x){ tmp <- unlist(strsplit(x,","))
                                          paste(tmp[!grepl("-",tmp)], collapse=",")
                              }) %>%
              unlist()
 
    if(labelUsing %in% c("positivemarkers","")){
        return(combos[!is.na(combos)])
    } else if(labelUsing == "tag+numberpositives"){
        return(paste(tags, 
                     sapply(combos, function(x){ length(unlist(strsplit(x,","))) }),
                     sep=""
                )
         )
    } 
}


#' Parse cell types XLSX file
#' 
#' Parse cell types XLSX file, expanding cell types and assigning labels when necessary
#' 
#' @param cellTypesXLSX   xlsx file in cell types format (see docs)
#' 
#' @return table containing parsed/expanded cell types to be used for cell annotation
#' @export
getCellTypes <- function(cellTypesXLSX){

    if(!file.exists(cellTypesXLSX) || file.size(cellTypesXLSX) == 0){
        stop(paste0("Cell Types file ", cellTypesXLSX, " is missing or empty."))
    }

    cellTypes   <- as_tibble(read.xlsx(cellTypesXLSX,1))
    toExpand    <- cellTypes %>% filter(Pos_required != 'all')
    allExpanded <- cellTypes %>% filter(Pos_required == 'all')

    for(te in 1:nrow(toExpand)){
        curRow <- toExpand[te,]
        markers <- curRow %>% pull(Pos_markers)
        into <- NULL
        into <- ifelse(curRow$Pos_required == "+" & 
                          (is.na(curRow$Combo_label) | grepl("number positive",curRow$Combo_label)), 
                       "allCombinations", 
                       "singleMarkers")
        combos <- expandMarkerCombos(markers, into=into, emptyCombo=FALSE)
        neg <- curRow$Neg_markers[!is.na(curRow$Neg_markers)]
        expanded <- combos %>%
                    lapply(., function(x){
                                   x <- sort(unlist(strsplit(x,",")))
                                   exp <- tibble(Pos_markers = paste(x[!grepl("-",x)],
                                                                     collapse=","),
                                          Neg_markers = paste(c(neg, x[grepl("-",x)]),
                                                            collapse=","))
                                   bind_cols(curRow %>% select(-c(Pos_markers,Neg_markers)), exp)
                              }
                    ) %>%
                    bind_rows()
        empt <- which(is.null(expanded$Neg_markers) | expanded$Neg_markers == "")
        expanded$Neg_markers[!empt] <- paste(expanded$Neg_markers[!empt], curRow$Neg_markers, sep=",")
        expanded$Neg_markers[empt] <- curRow$Neg_markers

        labelUsing <- curRow %>% pull(Combo_label) %>% tolower()
        if(!labelUsing %in% c("tag","tags") && !is.na(labelUsing)){
            expanded$Tag <- labelCellTypes(expanded$Pos_markers, labelUsing=labelUsing, tags=expanded$Tag)
        }

        allExpanded <- bind_rows(allExpanded, expanded)

    }

    allExpanded$Neg_markers <- gsub("-","",allExpanded$Neg_markers)


    #### LATE MODIFICATIONS
    rw <- which(allExpanded$Cell_type == "B cell")
    allExpanded[rw, c("Subtype", "Abbrev_label_for_figures", "Subscript")] <- c("All", "B", "All")
    rw <- which(allExpanded$Cell_type == "Other leukocyte")
    allExpanded[rw, c("Abbrev_label_for_figures", "Subscript")] <- c("Leuk", "Other")
    rw <- which(allExpanded$Cell_type == "Tumor")
    allExpanded[rw, c("Subtype", "Abbrev_label_for_figures", "Subscript")] <- c("All", "Tumor", "All")
    rw <- which(allExpanded$Subtype == "Natural killer cell")
    allExpanded[rw, c("Abbrev_label_for_figures", "Subscript")] <- c("NK", "")
    rw <- which(allExpanded$Subtype == "Natural killer T cell")
    allExpanded[rw, c("Abbrev_label_for_figures", "Subscript")] <- c("NKT", "")      

 
    allCT <- tibble(Category = "Immune",
                    Cell_type = c("Immune", "T cell", "Natural killer cell overall",
                                 "Macrophage/monocyte"),
                   Subtype = "All",
                   Abbrev = c("Immune", "T", "NK", "MΦ"),
                   Subscript = c("All", "All", "All", "All"))

    allExpanded <- allExpanded %>%
                   mutate(Abbrev = getAbbrev(Abbrev_label_for_figures, Subscript)) %>%
                   select(-Abbrev_label_for_figures) %>%
                   bind_rows(allCT) %>%
                   select(Category, Cell_type, Subtype, Tag, Abbrev, Subscript) %>%
                   unique()

    allExpanded$Abbrev[allExpanded$Abbrev == "MΦ+"] <- "MΦ"
    ########################

    return(allExpanded)
}


#' Add class(es) to pre-existing classes in Classifier column of annotated cell data
#' 
#' Given a vector of character strings where each element is a delimited collection of
#' classifiers for a single cell, add to every element a delimited string of new classes
#'
#' @param  existingClasses   a character vector where each element is a delimited (default delim=";")
#'                           string of classes to which a single cell belongs
#' @param  newClasses        a delimited character string of new classes to add to every item in 
#'                           {existingClasses}; NOTE: setting unique at this point will NOT 
#' @param  delim             delimiter separating individual classes in both {existingClasses} and {newClasses}
#' @return  vector of expanded class lists
#' @export
appendClassifiers <- function(existingClasses = NULL, newClasses, unique=FALSE, delim=";"){

    exC <- existingClasses
    if(!is.null(exC) && !is.na(exC)){
        exC <- unlist(strsplit(existingClasses, delim))
    } 
    if(unique){
        newClasses <- unlist(strsplit(newClasses, delim))
        nwC <- paste(unique(c(exC, newClasses)), collapse=delim)
    } else {
        nwC <- paste(c(exC, newClasses), collapse=delim)
    }
    return(nwC)
}

#' Get vector of all markers used in study
#'
#' Get vector of all markers used in study
#'
#' @param markerDesc  read from Markers meta file
#' @return vector of all markers under Marker_name column of markers file
getAllMarkers <- function(markerDesc){
    markerDesc %>% pull(Marker_name)
}

#' Get vector of all 'cell type' markers used in study
#'
#' Get vector of all 'cell type' markers used in study
#'
#' @param markerDesc  read from Markers meta file
#' @return vector of marker names identified as a cell type marker by 
#'         a 'X' under the Cell_type column of markers file 
getCellTypeMarkers <- function(markerDesc){
    markerDesc %>% filter(Cell_type == "X") %>% pull(Marker_name)
}

#' Get vector of all 'identity' markers used in study
#'
#' Get vector of all 'identiry' markers used in study
#'
#' @param markerDesc  read from Markers meta file
#' @return vector of markers classified as 'Identity' markers by the word 'Identity'
#'         in the Description column of markers file 
getIdentityMarkers <- function(markerDesc){
    markerDesc %>% filter(Description == "Identity") %>% pull(Marker_name)
}

#' Add columns to indicate all positive markers and positive cell identity markers
#'
#' Given tibble where each row is a unique cell and columns for each marker have
#' a value of 1 (marker is positive in that cell) or 0 (marker is negative in that cell),
#' add columns 'PositiveMarkers' and 'CTMarkers'. 'PositiveMarkers' is a comma-delimited
#' string of ALL positive markers in a cell and 'CTMarkers' is a comma-delimited string
#' of ONLY cell identity markers
#'
#' @param dat         tibble where each row is a unique cell and columns for each marker have
#'                    a value of 1 (marker is positive in that cell) or 0 (marker is negative 
#'                    in that cell)
#' @param markerDesc  tibble with at minimum columns Marker_name, Cell_type; markers where Cell_type == 'X'
#'                    will be included in values under CTMarkers
#' @return  original data table with PositiveMarkers and CTMarkers columns added
addCellPositiveMarkers <- function(dat, markerDesc){

    markers   <- getAllMarkers(markerDesc)
    ctMarkers <- getCellTypeMarkers(markerDesc) 
    idMarkers <- getIdentityMarkers(markerDesc)

    if(!any(markers %in% names(dat)) && all(c("Marker", "Value") %in% names(dat))){
        dat <- dat %>% spread(Marker, Value)
    }

    if(!"PositiveMarkers" %in% names(dat)){
        dat$PositiveMarkers <- apply(dat[,which(names(dat) %in% markers)], 1,
                                      function(x){
                                         paste(sort(names(x)[x==1]), collapse=",")
                                      })

        dat$CTMarkers <- unlist(lapply(dat$PositiveMarkers,
                                         function(x){ 
                                             x <- unlist(strsplit(x,","))
                                             paste0(x[x %in% ctMarkers], collapse=",")
                                       }))

        dat$IDMarkers <- unlist(lapply(dat$PositiveMarkers,
                                         function(x){
                                             x <- unlist(strsplit(x,","))
                                             paste0(x[x %in% idMarkers], collapse=",")
                                       }))
    }

    dat

}

#' Throw error if attempting to assign one or more classifications to 
#' cells that already contain conflicting classes
#'
#' Given a vector of existing cell class assignments (and/or NAs), if any classes
#' exist and are different from the ones to be added, throw an error 
#'
#' @param dat      vector of cell classifications to be checked for conflicts
#' @param classes  a vector of classes named by their classification type (Category, 
#'                 Cell_type, Subtype, etc.) to be assigned to cells
#' @return nothing
cellClassConflicts <- function(dat, classes){
    for(cl in names(classes)){  ## Category Cell_type Subtype Tag
        if(any(!is.na(dat[[cl]]))){
            vals <- unique(dat[[cl]])
            if(vals[!is.na(vals)] != classes[cl]){
                log_error(paste("Conflict betwen", cl, vals[!is.na(vals)], "and", classes[cl]))
                stop(paste("Conflict betwen", cl, vals[!is.na(vals)], "and", classes[cl]))
            }
        }
    }
}


#' Set all class columns in data table to NA
#' 
#' If class columns exist in data, make sure all data
#' in those columns are removed. If they do not exist, initialize
#' them to NA
#'
#' @param  dat         cell-level data table
#' @param  classTypes  column names of all classification types to be assigned
#'                     (generally Category, Cell_type, Subtype, Tag, Classifiers)
#' @return data table with empty columns for all class types
resetDataClassification <- function(dat, classTypes){
    for(type in classTypes){
        dat[[type]] <- as.character(NA)
    }
    dat
}

#' Get all marker combinations that are not sufficient to identify a cell type
#' 
#' Cells with a postive marker combination that contains no cell type markers and
#' no identity markers are classified as "Unknown super negative" cells. Pull out
#' the unique set of these combinations from halo data
#' 
#' @param dat       cell level halo data including columns CTMarkers and PositiveMarkers
#' @param idMarkers vector of Identity markers pulled from marker descriptions xlsx file
#' @return vector of unique marker combinations existing in data that are not sufficient
#'         to identify the cell as a certain type
getAllUnknownSupernegCombos <- function(dat, idMarkers){
    dat %>%
    filter(IDMarkers == "", !grepl(getClassifierPattern(idMarkers, delim=","), PositiveMarkers)) %>%
    select(PositiveMarkers) %>% 
    unique() %>% 
    pull()
}

#' Add all cell classifiers for each cell to Classifiers column
#'
#' Based on cellTypes table, determine all proper classifiers for each cell. Add columns for
#' Category, Cell_type, Subtype, Tag, and Classifiers (aggregation of all classifiers)
#' 
#' @param  dat        tibble where each row is a unique cell and PositiveMarkers column 
#'                    has already been added
#' @param  cellTypes  pre-processed (expanded) table of cell type and cell state definitions (see docs for details)
#' @param  markerDesc table of marker information including Marker_name, Cell_type and Description (see docs for details)
#' @param  reclassify logical; when TRUE, existing Classifiers column will be set to NULL before classifying cells; 
#'                    otherwise, classifiers will be appended to Classifiers column
#' @return newly-annotated cell data table 
addCellTypes <- function(dat, cellTypes, markerDesc, reclassify=FALSE){

    markers <- getAllMarkers(markerDesc) 
    ctMarkers <- getCellTypeMarkers(markerDesc)
    idMarkers <- getIdentityMarkers(markerDesc)

    if(!"PositiveMarkers" %in% names(dat)){
        log_error("Data must contain 'PositiveMarkers' column before cell types can be added")
        return(dat)
    }

    if(reclassify || !'Classifiers' %in% names(dat)){ 
        dat <- resetDataClassification(dat, c("Category", "Cell_type", "Subtype", "Tag", "Classifiers"))
    }

    for(x in 1:nrow(cellTypes)){
        curRow <- cellTypes[x,]

        pos <- unlist(strsplit(curRow$Pos_markers,","))
        neg <- unlist(strsplit(curRow$Neg_markers,","))

        tmp <- dat
        ## filter for cells that do NOT contain any negative markers
        tmp <- tmp %>%
               filter(!grepl(getClassifierPattern(neg, delim=","), PositiveMarkers))

        ## filter for cells that contain ALL positive markers
        for(p in pos){
            tmp <- tmp %>%
                   filter(grepl(getClassifierPattern(p, delim=","), PositiveMarkers))
        }

        idxs <- which(dat$UUID %in% tmp$UUID) 

        classes <- unlist(curRow %>% select(Category, Cell_type, Subtype, Tag))

        if(curRow$Classification_type == "type"){
    
            cellClassConflicts(dat[idxs,], classes)

            dat$Category[idxs]    <- curRow$Category
            dat$Cell_type[idxs]   <- curRow$Cell_type
            dat$Subtype[idxs]     <- curRow$Subtype
            dat$Tag[idxs]         <- curRow$Tag
        }

        idxs <- which(dat$UUID %in% tmp$UUID & !is.na(dat$Classifiers))
        if(length(idxs) > 0){
            dat$Classifiers[idxs] <- unlist(sapply(dat$Classifiers[idxs], function(x){ 
                                         appendClassifiers(x, classes, unique=FALSE)
                                     }))
        }

        newClassifierIdxs <- which(dat$UUID %in% tmp$UUID & is.na(dat$Classifiers))
        if(length(newClassifierIdxs) > 0){
            dat$Classifiers[newClassifierIdxs] <- paste(classes, collapse=";")
        }
        log_debug(paste0(curRow$Subtype, ": ", length(unique(c(idxs,newClassifierIdxs))), " cells"))
    }

    dat$Classifiers[dat$PositiveMarkers %in% getAllUnknownSupernegCombos(dat, idMarkers)] <- "Unknown;superNeg"
    
    ## finally, mark any cells that do not have a cell type classifier as unreasonable
    ctClasses <- c("Immune","Tumor","Other","Unknown;superNeg")
    dat$Classifiers[!grepl(paste(paste0("^",ctClasses),collapse="|"), dat$Classifiers)] <- "Unreasonable;Unreasonable"

    dat
}


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
#'                        column names in dat; default=c("Sample_ID","FOV_ID")
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
#' Summarize cell counts by cell types and data level (e.g., all, Sample_ID, FOV_ID)
#' 
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param cellTypes       table of cell type definitions
#' @param classCols       column names from cellTypes tibble indicating which categories should be 
#'                        summarized; default=c("Category","Cell_type","Subtype","Tag")
#' @param summarizeBy     vector of ways in which data should be grouped; values must be equal to
#'                        column names in dat; default=c("Sample_ID","FOV_ID")
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

checkGrouping <- function(cols, groupBy){
    if(!all(groupBy %in% cols)){
        stop(paste0("Invalid groups: ",groupBy))
    }
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


#' Get vector of UUIDs for cells that fall into a certain class
#' 
#' Get vector of UUIDs for cells that fall into a certain class; this is
#' especially helpful for ensuring that there is no overlap/double counting when
#' comparing two groups of cells
#'
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param class           class/cell type to find
#' @return vector of UUIDs of cells in class
#' @export
getClassUUIDs <- function(dat, class){
    dat %>% 
    filter(grepl(getClassifierPattern(class), Classifiers)) %>%
    pull(UUID)
}

#' Build regex pattern needed to properly pull cell types out of Classifiers column
#' 
#' Build regex pattern needed to properly pull cell types out of Classifiers column; this
#' is to ensure special characters in the class name are escaped and that the class and 
#' ONLY the class are being counted
#' 
#' @param c   character string of class being searched for
#' @return  a character string to be used when grepping for class in Classifiers column of data
#' @export
getClassifierPattern <- function(c, delim=";"){

    pat <- paste(c(paste0("^",c,delim),paste0(delim,c,delim),paste0(delim,c,"$"),paste0("^",c,"$")), collapse="|")
    pat <- gsub("\\+","\\\\+",pat)
    pat <- gsub("\\(","\\\\(",pat)
    pat <- gsub("\\)","\\\\)",pat)

    return(pat)
}

#' Check that there are no overlaps in UUIDs of cells in multiple groups
#' 
#' In order to ensure there is no double counting of cells, check for overlaps
#' in lists of UUIDs in each group being counted
#' 
#' @param UUIDlist    list containing vectors of UUIDs, each representing a cell group
#'                    that should have no intersections with any other vectors in the list
#' @return list where each element is a pair of vectors that contains one or more overlapping cells
#' @export 
cellGroupOverlaps <- function(UUIDlist){

    allOverlaps <- list()

    for(x in 1:(length(UUIDlist)-1)){ 
        s1 <- UUIDlist[[x]]
        tmp <- UUIDlist[-x]
        for(y in 1:length(tmp)){
            s2 <- tmp[[y]]
            overlap <- intersect(s1,s2)
            if(!is.null(overlap) && length(overlap) > 0){
                if(!is.null(names(UUIDlist))){
                    allOverlaps[[length(allOverlaps) + 1]] <- c(names(UUIDlist)[x], names(tmp)[y])
                } else {
                    idx2 <- which(unlist(lapply(UUIDlist, function(lst){ identical(lst, s2) })))
                    allOverlaps[[length(allOverlaps) + 1]] <- c(x, idx2) 
                }
            }
        }
    }

    allOverlaps
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
#'                        column names in dat; default=c("Sample_ID","FOV_ID")
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

#' Get vector of UUIDs that match one or more class/cell type
#' 
#' Get vector of UUIDs that match one or more class/cell type
#' 
#' @param dat        table of halo data where a row represents a single cell and contains,
#'                   at minimum, columns UUID and Classifiers
#' @param cellTypes  vector of one or more cell types to search for in Classifiers column
#' @return  vector of all UUIDs that match one or more value in {cellTypes}
#' @export
uuidsAnyClass <- function(dat, cellTypes){
    fullPat <- getClassifierPattern(cellTypes[1])
    if(length(cellTypes) > 1){
        for(ct in cellTypes[2:length(cellTypes)]){
            fullPat <- paste(fullPat, getClassifierPattern(ct), sep="|")       
        }
    }    
    dat %>% filter(grepl(fullPat, Classifiers)) %>% pull(UUID)
} 

#' Get vector of UUIDs that match ALL classes/cell types
#' 
#' Get vector of UUIDs that match ALL classes/cell types
#' 
#' @param dat        table of halo data where a row represents a single cell and contains,
#'                   at minimum, columns UUID and Classifiers
#' @param cellTypes  vector of one or more cell types to search for in Classifiers column
#' @return  vector of all UUIDs that match ALL values in {cellTypes}
#' @export
uuidsAllClasses <- function(dat, cellTypes){
    tmp <- dat    
    for(ct in cellTypes){
        tmp <- tmp %>% filter(grepl(getClassifierPattern(ct), Classifiers))
    }
    if(nrow(tmp) > 0){
        return(tmp$UUID)
    } else {
        return(NULL)
    }
}

#' Separate Classifiers column into individual category columns according to *CellTypes.xlsx file
#' 
#' Separate Classifiers column in annotated cells file into columns described
#' in cell types table
#' 
#' @param annDat    table of annotated cell data where each row represents one unique cell and 
#'                  contains a Classifiers column containing all classifiers assigned to a cell, 
#'                  delimited by semi-colons
#' @param cellTypes parsed and expanded table version of *CellTypes.xlsx
#' @param classCols vector of columns from cellTypes table; default=c("Category", "Cell_type", "Subtype", "Tag")
#' @return  annDat table with columns for Category, CellType, Subtype and Tag
#' @export
separateClassifiers <- function(annDat, cellTypes, classCols=NULL){

    classCols <- c("Category","Cell_type","Subtype","Tag")

    for(col in classCols){
        annDat[[col]] <- NA
        cts <- unique(cellTypes[[col]])
        for(ct in cts){
            uuids <- uuidsAllClasses(annDat, ct)
            annDat[[col]][annDat$UUID %in% uuids] <- ct
        }
    }
    annDat

}

#' Get all possible combinations of markers
#'
#' Given a vector of single markers, return a list
#' of all possible combinations of those markers,
#' including both positive and negative variations
#' 
#' @param markers  vector of single markers
#' @return  list of all combinations
#' @export
getAllCombos <- function(markers){
    all.combos = list()
    ## get all combinations of markers
    for(x in 1:length(sort(markers))){
        ## get all combos of x
        combos = combn(markers,x,simplify=FALSE)
        all.combos = c(all.combos, combos)
    }

    for(c in 1:length(all.combos)){
        combo = all.combos[[c]]
        if(length(combo) < length(markers)){
            all.combos[[c]] <- c(combo,paste0(setdiff(markers,combo),"-"))
        }
        all.combos[[c]] <- paste0(all.combos[[c]],collapse=",")
    }
    return(all.combos)
}

#' Extract non-subscript portion of cell type abbreviation
#' 
#' Extract non-subscript portion of cell type abbreviation
#' 
#' @param longAbbrev  character string representing complete form of cell type abbreviation
#' @param subscript   character string representing the part of logAbbrev that is to be 
#'                    subscripted whenever possible
getAbbrev <- function(longAbbrev, subscript){
    sapply(1:length(longAbbrev), function(x){
        gsub(subscript[x], "", longAbbrev[x])
    }) %>% unlist()
}

annotateCells <- function(annotatedCellsFile, dataDir = NULL, dataFiles = NULL, 
                          metaFiles = NULL, metaDataFile = NULL, numThreads = 1, 
                          filterExclusions = FALSE, controlMarker = "DAPI"){

    reclassify <- FALSE

    ###
    ### load data
    ###
    allDat <- NULL
    if(!fileDone(annotatedCellsFile)){
        
        annot <- loadStudyAnnotations(metaFiles = metaFiles, metaDataFile = metaDataFile)
        idMap <- annot$IDs

        ### load halo object analysis data if needed
        log_debug("Loading object analysis data...")
        dataFiles <- getFiles(path = dataDir, files = dataFiles, pattern=".rda")
        nThreads <- ifelse(length(dataFiles) < numThreads, length(dataFiles), numThreads)

        allDat <- loadAllHaloData(dataFiles = dataFiles,
                                  nThreads = nThreads,
                                  filterExclusions = TRUE, 
                                  controlMarker = controlMarker) %>%
                  joinIDs(idMap) %>%
                  spread(Marker, Value)

        log_debug(paste0("Starting with ",nrow(allDat)," ", controlMarker, "+ cells"))

    } else {
        allDat <- readRDS(annotatedCellsFile)
        reclassify <- TRUE
    }

    log_debug("Reading meta data...")
    ctFile    <- metaFiles[grep("_CellTypes.xlsx", metaFiles)]
    cellTypes <- openxlsx::read.xlsx(ctFile, 1, check.names = F) %>% as_tibble()

    mrkrFile  <- metaFiles[grep("_Markers.xlsx", metaFiles)]
    markerDesc <- openxlsx::read.xlsx(mrkrFile, 1, check.names = F) %>% as_tibble()

    ###
    ### add column containing comma-delimited list of markers positive in each cell
    ###
    log_debug("Adding positive marker columns...")
    annDat <- addCellPositiveMarkers(allDat, markerDesc)

    ###
    ### annotate cell level data with cell type, subtype, etc. 
    ###
    log_debug("Adding cell annotation...")
    annDat <- addCellTypes(annDat, cellTypes, markerDesc, reclassify=reclassify)

    ###
    ### organize and filter columns
    ###
    annDat <- annDat %>%
              select(UUID, CellDive_ID, Patient_ID, Sample_ID, Lesion_ID, FOV_ID, 
                     XMin, XMax, YMin, YMax, PositiveMarkers, CTMarkers, IDMarkers, 
                     Category, Cell_type, Subtype, Tag, Classifiers)

    log_info(paste0("Saving newly annotated cells in file: ", annotatedCellsFile))
    ## save data
    saveRDS(annDat, annotatedCellsFile)

    annDat

}
