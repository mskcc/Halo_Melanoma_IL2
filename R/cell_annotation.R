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
                                   curRow %>%
                                   mutate(Pos_markers = paste(x[!grepl("-",x)],
                                                              collapse=","),
                                          Neg_markers = paste(c(neg, 
                                                                gsub("\\-", "", x[grepl("-",x)])),
                                                                collapse=","))
                              }
                    ) %>%
                    bind_rows()
        empt <- which(is.null(expanded$Neg_markers) | expanded$Neg_markers == "")
        expanded$Neg_markers[!empt] <- paste(expanded$Neg_markers[!empt], curRow$Neg_markers, sep=",")
        expanded$Neg_markers[empt] <- curRow$Neg_markers

        labelUsing <- curRow %>% pull(Combo_label) %>% tolower()
        if(!labelUsing %in% c("tag","tags") && !is.na(labelUsing)){
            expanded$Tag <- labelCellTypes(expanded$Pos_markers, 
                                           labelUsing=labelUsing, 
                                           tags=expanded$Tag)
        }

        allExpanded <- bind_rows(allExpanded, expanded)

    }

    #### LATE MODIFICATIONS
    rw <- which(allExpanded$Cell_type == "B cell")
    allExpanded[rw, c("Subtype", "Abbrev_label_for_figures", "Subscript")] <- list("All", "B", "All")
    rw <- which(allExpanded$Cell_type == "Other leukocyte")
    allExpanded[rw, c("Abbrev_label_for_figures", "Subscript")] <- list("Leuk", "Other")
    rw <- which(allExpanded$Cell_type == "Tumor")
    allExpanded[rw, c("Subtype", "Abbrev_label_for_figures", "Subscript")] <- list("All", "Tumor", "All")
    rw <- which(allExpanded$Subtype == "Natural killer cell")
    allExpanded[rw, c("Abbrev_label_for_figures", "Subscript")] <- list("NK", "")
    rw <- which(allExpanded$Subtype == "Natural killer T cell")
    allExpanded[rw, c("Abbrev_label_for_figures", "Subscript")] <- list("NKT", "")      

 
    allCT <- tibble(Category = "Immune",
                    Cell_type = c("Immune", "T cell", "Natural killer cell overall",
                                 "Macrophage/monocyte"),
                    Subtype = "All",
                    Abbrev = c("Immune", "T", "NK", "MΦ"),
                    Subscript = c("All", "All", "All", "All"),
                    Classification_type = "type")

    allExpanded <- allExpanded %>%
                   mutate(Abbrev = getAbbrev(Abbrev_label_for_figures, Subscript)) %>%
                   select(-Abbrev_label_for_figures) %>%
                   bind_rows(allCT) 

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
#'
#' @return  vector of expanded class lists
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
#'
#' @return vector of all markers under Marker_name column of markers file
getAllMarkers <- function(markerDesc){
    markerDesc %>% pull(Marker_name)
}


#' Get vector of all 'cell type' markers used in study
#'
#' Get vector of all 'cell type' markers used in study
#'
#' @param markerDesc  read from Markers meta file
#'
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
#'
#' @return vector of markers classified as 'Identity' markers by the word 'Identity'
#'         in the Description column of markers file 
getIdentityMarkers <- function(markerDesc){
    markerDesc %>% filter(Description == "Identity") %>% pull(Marker_name)
}


#' Add columns to indicate all positive markers and positive cell identity markers
#'
#' Given a tibble where each row is a unique cell and columns for each marker have
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
#'
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
#'
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
#'
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
#' 
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
#'
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


#' Check that all group names exist in a set of column names 
#' 
#' Throw an error if not all column names to be used as grouping
#' variables exist in set of column names pulled from data
#'
#' @param cols     vector of tibble column names
#' @param groupBy  groupBy vector of group names
#'
#' @return nothing; throw error if not all group names are in col names
checkGrouping <- function(cols, groupBy){
    if(!all(groupBy %in% cols)){
        stop(paste0("Invalid groups: ",groupBy))
    }
}


#' Get vector of UUIDs for cells that fall into a certain class
#' 
#' Get vector of UUIDs for cells that fall into a certain class; this is
#' especially helpful for ensuring that there is no overlap/double counting when
#' comparing two groups of cells
#'
#' @param dat             tibble of *annotated* halo data containing one or more samples
#' @param class           class/cell type to find
#'
#' @return vector of UUIDs of cells in class
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
#'
#' @return  a character string to be used when grepping for class in Classifiers column of data
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
#' 
#' @return list where each element is a pair of vectors that contains one or more overlapping cells
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


#' Get vector of UUIDs that match one or more class/cell type
#' 
#' Get vector of UUIDs that match one or more class/cell type
#' 
#' @param dat        table of halo data where a row represents a single cell and contains,
#'                   at minimum, columns UUID and Classifiers
#' @param cellTypes  vector of one or more cell types to search for in Classifiers column
#'
#' @return  vector of all UUIDs that match one or more value in {cellTypes}
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
#'
#' @return  vector of all UUIDs that match ALL values in {cellTypes}
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
#'
#' @return  annDat table with columns for Category, CellType, Subtype and Tag
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
#'
#' @return  list of all combinations
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

#' Remove any markers determined to be problematic after analysis
#' 
#' Most exclusions are determined prior to rethresholding, but this function
#' allows users to exclude markers determined to be problematic during
#' QC or any other analysis step.
#'
#' @param dat    data tibble loaded from halo object analysis RDA files, joined
#'               with all sample and FOV IDs
#' @param annot  flat annotation data returned from loadStudyAnnotations(),
#'               including columns FOV_exclusion (see docs for 
#'               details on format)
#' @param return data tibble with any additional or markers removed 
postAnalysisMarkerExclusions <- function(dat, annot){
   
    if(all(is.na(annot$Marker_exclusion) | annot$Marker_exclusion == "")){ return(dat) }

    ## filter out markers
    mx <- lapply(annot$Marker_exclusion[!is.na(annot$Marker_exclusion) & annot$Marker_exclusion != ""], 
                 function(x) {
                   length(unlist(strsplit(x, ",")))
                 }
                ) %>% unlist %>% max
    sepCols <- paste0("Marker.", seq(1:mx))

    mExcl <- annot %>%
             filter(!is.na(Marker_exclusion), Marker_exclusion != "") %>%
             select(CellDive_ID, FOV_number, Marker_exclusion) %>%
             separate(Marker_exclusion, sepCols, sep=",") %>%
             gather(all_of(sepCols), key = "tmp", value = Marker) %>%
             select(-tmp) %>%
             filter(!is.na(Marker)) %>%
             unique %>%
             mutate(REMOVE = "YES", Marker = gsub(" ", "", Marker))

    tmp  <- dat %>% left_join(mExcl, by = intersect(names(.), names(mExcl)))

    ## log what is to be removed
    excl <- tmp %>%
            filter(REMOVE == "YES") %>%
            group_by(CellDive_ID, FOV_number, Marker) %>%
            summarize(Count = n()) %>%
            mutate(LOG = paste("Removed marker [", Marker, "] from", Count, 
                               "cells in FOV", FOV_number, "of sample", CellDive_ID))
    sapply(1:nrow(excl), function(x) log_debug(excl$LOG[x]))

    tmp %>%
    filter(is.na(REMOVE)) %>%
    select(-REMOVE)
}

#' Remove any FOVs determined to be problematic after analysis
#' 
#' Most exclusions are determined prior to rethresholding, but this function
#' allows users to exclude FOVs determined to be problematic during
#' QC or any other analysis step.
#'
#' @param dat    data tibble loaded from halo object analysis RDA files
#' @param annot  flat annotation data returned from loadStudyAnnotations(),
#'               including columns FOV_exclusion (see docs for 
#'               details on format)
#' @param return data tibble with any additional FOVs removed 
postAnalysisFOVexclusions <- function(dat, annot){

    totalCellsRemoved <- 0

    ## filter out FOVs
    fExcl <- annot %>%
             filter(!is.na(FOV_exclusion), FOV_exclusion != "") %>%
             select(CellDive_ID, FOV_number) %>%
             unique()

    excl <- dat
    if(nrow(fExcl) != 0){
        for(x in 1:nrow(fExcl)){
            rmv <- excl %>% 
                   filter(CellDive_ID == fExcl$CellDive_ID[x],
                          FOV_number == fExcl$FOV_number[x]) %>%
                   pull(UUID)
   
            log_debug(paste0("Removing ", length(rm), "cells",
                             " from FOV ", fExcl$FOV_number[x], 
                             " of sample ", fExcl$CellDive_ID[x]))  
        
            totalCellsRemoved <- totalCellsRemoved + length(rmv)

            excl <- excl %>% filter(!UUID %in% rmv)
        }
    }
    log_debug(paste0("Removed ", nrow(fExcl), " entire FOVs for a total of ", totalCellsRemoved, " cells"))
    
    excl
}

#' Add to cell-level tibble columns for positive markers and cell classifications
#' 
#' Given a tibble complete table of marker positivity for each cell, add columns
#' consolidating list of positive markers in each cell, positive cell identity markers
#' for each cell, and assign cell classes/types based on those positive markers
#' 
#' @param annotatedCellsFile   file that either already contains a cell-level tibble of
#'                             Halo data or the file to which said table will be saved
#' @param dataDir              directory of Halo files to be annotated; when NOT NULL, 
#'                             ALL RDA files in this directory will be loaded and annotated
#' @param dataFiles            vecotr of RDA files to be annotated
#' @param metaFiles            vector of meta data XLSX files, including sample annotation
#'                             marker information and cell type definitions files (default = NULL)
#' @param metaDataFile         RDA file of pre-compiled meta data (default = NULL)
#' @param numThreads           integer; number of threads
#' @param filterExclusions     logical indicating whether to remove cells marked with any text
#'                             in 'EXCLUDE' column; default = FALSE
#' @param controlMarker        marker whose negativity indicated the cell is not usable; these
#'                             cells are removed
#'
#' @return all annotated data
annotateCells <- function(annotatedCellsFile, dataDir = NULL, dataFiles = NULL, 
                          metaFiles = NULL, metaDataFile = NULL, numThreads = 1, 
                          filterExclusions = FALSE, controlMarker = "DAPI", 
                          forceReannotation = FALSE){

    reclassify <- FALSE

    ###
    ### load data
    ###
    allDat <- NULL
    if(!fileDone(annotatedCellsFile) || forceReannotation){
        
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
                  joinIDs(annot$IDs) %>%
                  postAnalysisMarkerExclusions(annot$flat) %>%
                  postAnalysisFOVexclusions(annot$flat) %>%
                  spread(Marker, Value)

        log_debug(paste0("Starting with ",nrow(allDat)," ", controlMarker, "+ cells"))

    } else {
        allDat <- readRDS(annotatedCellsFile)
        reclassify <- TRUE
    }

    log_debug("Reading meta data...")
    ctFile    <- metaFiles[grep("_CellTypes.xlsx", metaFiles)]
    cellTypes <- getCellTypes(ctFile)

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
    annDat <- addCellTypes(annDat, cellTypes, markerDesc)

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

#' Determine whether a cell meets the definition of a specific phenotype
#' 
#' Given a comma-delimited string representing a cell phenotype, in the form 
#' [Cell Type/Subtype/Tag],[Marker1],[Marker2],[Marker3], and the cells' cell type
#' classifiers and positive markers, return a 1 if the cell fits the phenotype and
#' a 0 if it does not.
#' 
#' @param phenStr        comma-delimited string where the first element contains a cell
#'                       type classifier (Cell_type, Subtype, Category or Tag) and the 
#'                       remaining elements are any combination of positive and negative 
#'                       markers
#' @param classifierStr  semi-colon delimited string containing all classifiers that 
#'                       apply to a cell
#' @param posMarkerStr   comma-delimited string of positive and negative markers
#' @param nameMap        list of out-dated classifiers named by the current ones used
#'                       in phenStr
cellFitsPhenotype <- function(phenStr, classifierStr, posMarkerStr, nameMap = NULL){

    if(length(classifierStr) != length(posMarkerStr)){
        stop("Vectors of classifier strings and positive marker strings differ in length.")
    }

    class <- gsub(",.*", "", phenStr)
    if(class %in% names(nameMap)){ class <- nameMap[[class]] }

    mrkrs <- unlist(strsplit(phenStr, ","))[-1]
    pos <- mrkrs[!grepl("\\-", mrkrs)]
    neg <- mrkrs[grepl("\\-", mrkrs)]

    lapply(1:length(classifierStr), function(x){


        if(!grepl(getClassifierPattern(class), classifierStr[x])){
            return(0)
        }

        if(length(pos) > 0){
            for(pm in pos){
                if(!grepl(getClassifierPattern(pos, delim = ","), posMarkerStr[x])){
                    return(0)
                }
            }
        }

        if(length(neg) > 0){
            if(grepl(getClassifierPattern(neg, delim = ","), posMarkerStr[x])){
                return(0)
            }
        }

        1

    }) %>% unlist()
}


