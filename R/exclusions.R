#' Remove points from data set that fall inside exclusion boundaries
#' 
#' Read Halo boundaries annotations file and remove from data set any cells
#' that fall inside the exclusion regions
#' 
#' @param ds     tibble containing at minimum, Sample, SPOT, UUID, X, Y, Marker, Value
#' @param excB   a list of excluded region tables as returned by cleanBoundaries;
#'               if NULL, must provide the original Halo annotations file in XML format
#' @param aFile  Halo boundaries annotation file in XML format; if NULL, must provide
#'               the list of excluded region tables as returned by cleanBoundaries()
#' @return tibble filtered to remove points that fall inside exclusion boundaries
getHaloExclusions <- function(fov, dd, excType, haloAnnotations=NULL, aFile=NULL, exc=NULL, 
                              boundaryColors=NULL, boundaryReassignmentFile=NULL){

    ds <- dd %>% filter(SPOT == fov)
    if(!all(c('X','Y') %in% colnames(ds))){
        ds <- convertCellMinMaxToMidpoints(ds)
    }

    annCodes <- list(Exc='excB', Epi='epiB', Gls='glsB', Tum='tumB')
    if(is.null(exc)){
        if(is.null(haloAnnotations)){
            boundaries <- readHaloAnnotations(aFile, boundaryColors=boundaryColors, 
                                              boundaryReassignmentFile=boundaryReassignmentFile)
            allBoundaries <- cleanBoundaries(boundaries) #remove exc regions completely surrounded by another
            excB <- allBoundaries[[annCodes[[excType]]]]
        } else {
            excB <- haloAnnotations[[as.character(fov)]][[annCodes[[excType]]]]
            #excB <- cleanBoundaries(split(excB, excB$RegionNum))
        }
        if(is.null(excB)){
            return(c())
        }
#        #excB <- split(excB, excB$RegionNum)
    }

    spCells <- ds %>% dplyr::select(UUID,X,Y)
    coordinates(spCells) <- ~X+Y

    excludedCellIdx <- NULL
    if(length(excB)>0){
        spExcB <- seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})
        for(jj in 1:length(spExcB)) {
            excCellJJ <- unlist(sp::over(spExcB[[jj]],geometry(spCells),returnList=T))
            excludedCellIdx <- union(excludedCellIdx,excCellJJ)
        }
    }
    fovHaloIdx <- which(dd$SPOT == fov)
    fovHaloIdx[excludedCellIdx]

}


joinFOVexclusions <- function(dat, sampAnn){
    dat %>%
    left_join(sampAnn %>% select(CellDive_ID, FOV_number, FOV_exclusion),  
              by = c("CellDive_ID", "FOV_number")) %>%
    mutate(FOV_exclusion = ifelse(FOV_exclusion == "X", "FOV_exclusion", NA))
}

joinMarkerExclusions <- function(dat, sampAnn){
    cols <- sapply(sampAnn$Marker_exclusion, function(x) length(unlist(strsplit(x, ","))) ) %>% unlist() %>% max()

    mExcl <- sampAnn %>% 
             select(CellDive_ID, FOV_number, Marker_exclusion) %>%
             separate(Marker_exclusion, paste0("M_", seq(cols))) %>%
             gather(paste0("M_", seq(cols)), key="TMP", value="Marker") %>%
             select(-TMP) %>%
             filter(!is.na(Marker), Marker != "") %>%
             mutate(Marker = trimws(Marker), Marker_exclusion = "Marker_exclusion")

    dat %>%
    left_join(mExcl, by = c("CellDive_ID", "FOV_number", "Marker")) 
}

#' Add column to halo data indicating whether a cell should be excluded due
#' to drift/loss percentage being greater than the set threshold
#'
#' Mark cells to be excluded because the percentage of drift/loss is
#' greater than the set threshold
#' 
#' @param dat       cell data
#' @param drift     table of drift/loss percentages
#' @param threshold maximum percent of drift/loss allowed
#' @return vector of dat indices that need to be marked for exclusion
joinDriftExclusions <- function(dat, drift, threshold){

    drft <- drift %>% 
            select(Image_Location=image_location, XMin=x_min, XMax=x_max,
                   YMin=y_min, YMax=y_max, drift_loss_pixel_pct) %>%
            filter(drift_loss_pixel_pct > threshold) %>%
            mutate(Drift_exclusion = paste0("Drift_", threshold))

    dat %>% 
    left_join(drft, by=c("Image_Location","XMin","XMax","YMin","YMax")) %>%
    select(-Image_Location, -drift_loss_pixel_pct)
}

#' Get indices of cell to be excluded because they lie outside the border padding 
#' 
#' Get indices of cell to be excluded because they lie outside the border padding 
#' 
#' @param dat        data tibble containing data to be marked, including columns
#'                   X and Y, the coordinates for each cell
#' @param bbFOV      list of four elements, X0, X1, Y0, Y1, representing the 
#'                   bounding box of a single FOV
#' @param borderPad  number in pixels to trim from data
joinBorderPaddingExclusions <- function(dat, bbFOV, borderPad){

    bb <- list(X0 = bbFOV$X0 + borderPad,
               X1 = bbFOV$X1 - borderPad,
               Y0 = bbFOV$Y0 + borderPad,
               Y1 = bbFOV$Y1 - borderPad)

    dat %>%
    mutate(Padding_exclusion = ifelse(X < bb$X0 | X > bb$X1 | Y < bb$Y0 | Y > bb$Y1, 
                                      paste0("Padding_", borderPad * pixel2um, "um"), NA))

}

markManualExclusions <- function(dat, exclusionBounds){
    excDat <- tibble()

    ## mark cells that fall inside Halo exclusion boundaries (exclusions, epidermis, glass) and also
    ## those that fall outside limits set by border pad 
    for(fov in unique(dat$FOV_number)){
        ds <- dat %>% filter(FOV_number == fov) %>% mutate(Manual_exclusion = "")

        exc <- exclusionBounds[[as.character(fov)]]
        exc <- exc[!names(exc) == "tumB"]
        
        if(is.null(exc) || length(exc) == 0){
            excDat <- bind_rows(excDat, ds)
            next
        } 

        bCode <- c("epiB" = "epidermis", "excB" = "exclusion", "glsB" = "glass")
        for(boundType in names(exc)){
            excB <- exc[[boundType]]
            spCells <- ds %>% dplyr::select(UUID,X,Y)
            coordinates(spCells) <- ~X+Y
            excludedCellIdx <- NULL
            spExcB <- seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})
            for(jj in 1:length(spExcB)) {
                excCellJJ <- unlist(sp::over(spExcB[[jj]],geometry(spCells),returnList=T))
                excludedCellIdx <- union(excludedCellIdx,excCellJJ)
            }
            ds$Manual_exclusion[excludedCellIdx] <- paste0(ds$Manual_exclusion[excludedCellIdx], 
                                                           rep(bCode[[boundType]], length(excludedCellIdx)))
        } 
        excDat <- bind_rows(excDat, ds)
    }
    excDat %>% mutate(Manual_exclusion = ifelse(Manual_exclusion == "", NA, Manual_exclusion))
}

#' Add exclusion reasons to a vector of existing exclusion reasons
#' 
#' Given an exclusion reason and the indices of cells to be excluded, 
#' properly add the reason to vector of existing exclusions
#' 
#' @param exclCol   vector of existing exclusions (EXCLUDE column from object analysis data)
#' @param reason    string to be added to exclude column indicating why a cell is to be excluded
#' @param idxs     vector of indices at which exclusion reason should be added
#' @return vector of exclusion reasons with new reason added
writeExclusionReasons <- function(exclCol,idxs,reason){
    for(idx in idxs){
        exclCol[idx] <- ifelse(exclCol[idx] == "", reason, paste0(exclCol[idx],",",reason)) 
    }
    exclCol
}

#' Given a tibble of object analysis data, add a column EXCLUDE to indicate
#' which cells should be excluded from further analysis
#'
#' Given a tibble of object analysis data, add a column EXCLUDE to indicate
#' which cells should be excluded from further analysis, according to various upstream
#' factors including drift/loss, border padding, and other technical reasons
#' determined by investigators
#'
#' @param dat             tibble containing Halo object analysis data, to which EXCLUDE column
#'                        will be added
#' @param drift           tibble of drift/loss summary 
#' @param fovAnn          FOV annotations for one sample
#' @param exclusionBounds parsed Halo boundaries in list form for a single sample
#' @param borderPad       number in microns indicating the minimum distance between a cell and the FOV
#'                        border in order for that cell NOT to be excluded
#' @param driftThreshold  maximum pixel percent determined to have drifted in order to NOT be excluded
#' @return tibble matching {dat} exactly, with an extra column, EXCLUDE, added to indicate which
#'         cells should not be analyzed
#' @export
markExclusions <- function(dat, drift, fovAnn, exclusionBounds, driftThreshold=0.1, borderPad = 20){

    borderPad_px <- borderPad/pixel2um

    dat %>%
    mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) %>%
    joinFOVexclusions(fovAnn) %>%
    joinMarkerExclusions(fovAnn) %>%
    joinDriftExclusions(drift, driftThreshold) %>%
    joinBorderPaddingExclusions(getConstantBoundingBox(dat$CellDive_ID[1]), borderPad_px) %>%
    markManualExclusions(exclusionBounds) %>%
    unite("EXCLUDE", 
          c(FOV_exclusion, Marker_exclusion, Drift_exclusion, Padding_exclusion, Manual_exclusion),
          sep = ";", na.rm = T)

}


#' Remove specific FOVs from tibble based on meta data criteria
#' 
#' Remove specific FOVs from tibble based on meta data criteria; this
#' function is generally used for MID-analysis FOV exclusions, ones
#' that were not anticipated before analysis started
#' 
#' @param exclusions list of meta criteria describing FOVs/samples to 
#'                   be excluded; each name must be a column name 
#'                   from meta data and values are vectors of meta
#'                   data field names corresponding to specific data
#'                   to be removed
#' @return filtered data
#' @export 
excludeFOVs <- function(dat, studyAnnotations, exclusions=NULL){
    if(!is.null(exclusions)){
        rmv <- list()
        for(ex in names(exclusions)){
            rmv[[length(rmv) + 1]] <- studyAnnotations %>%
                                      filter_at(vars(ex), any_vars(grepl(paste0(exclusions[[ex]], collapse="|"), .))) %>%
                                      select(FOV)
        }
        rmv <- bind_rows(rmv) %>% unique()
        log_info(paste0("Excluding the following FOVs: ", paste(rmv$FOV, collapse=", ")))
        dat <- dat %>%
               filter(!FOV %in% rmv$FOV)
    }
    postAnalysisRmv <- studyAnnotations %>%
                       select(FOV, FOV_exclusion_post_analysis) %>%
                       filter(FOV_exclusion_post_analysis == "X") %>%
                       pull(FOV)
    if(length(postAnalysisRmv) > 0){
        dat <- dat %>% filter(!FOV %in% postAnalysisRmv)
    }
    return(dat)
}



#' Exclude only certain markers (as opposed to entire cells or FOVs)
#' 
#' Exclude certain markers by either removing the row(s) containing them under the Marker column
#' OR removing them from columns PositiveMarkers & MarkerPosTag (depending on where in the process 
#' the data table comes from)
#' 
#' @param dat   table containing either 'PositiveMarkers' column or 'Marker' column
#' @param studyAnnotations  table of study annotations, including Marker_exclusion column in which an 'X' 
#'                          indicates that a marker should be excluded from an FOV
#' @return if 'PositiveMarkers' in dat, a tibble where excluded markers have been removed from this column where
#'         appropriate; if 'Marker' column in dat, a tibble where rows for excluded markers have been removed where 
#'         appropriate
excludeMarkers <- function(dat, studyAnnotations){
    excl <- studyAnnotations %>%
            filter(!is.na(Marker_exclusion), Marker_exclusion != "") %>%
            select(CellDive_ID, FOV_number, Marker_exclusion) %>%
            as.data.frame()

    exclToBeMade <- apply(excl, 1, function(x){
                                       etbm <- dat %>% filter(CellDive_ID == x["CellDive_ID"],
                                                              FOV_number == as.numeric(x["FOV_number"]), 
                                                              grepl(gsub(",", "|", gsub(" ","", x["Marker_exclusion"])), PositiveMarkers))
                                       if(nrow(etbm) > 0){
                                           return(tibble(CellDive_ID = x["CellDive_ID"], 
                                                         FOV_number = as.numeric(x["FOV_number"]), 
                                                         Marker_exclusion = x["Marker_exclusion"]))
                                       }
                                   }) %>%
                    bind_rows()

    if(is.null(exclToBeMade) || nrow(exclToBeMade) == 0){
        log_info("All marker exclusions already made.")
        cat("All marker exclusions already made.\n")
        return(dat)
    } else {
        log_info("Excluding additional markers...")
        cat("Excluding additional markers...\n")
    }

    if("Marker" %in% names(dat)){
        for(x in 1:nrow(excl)){
            mEx <- unlist(strsplit(excl[x,"Marker_exclusion"]), ",")
            for(m in mEx){
                dat <- dat %>% filter(!(FOV == excl$FOV[x] && Marker == m))
            }
        }
    } else if("PositiveMarkers" %in% names(dat)){
        for(x in 1:nrow(exclToBeMade)){
            mEx <- unlist(strsplit(exclToBeMade[x,] %>% pull(Marker_exclusion), ","))
            for(m in mEx){
                exUUIDs <- dat %>%
                           filter(CellDive_ID == exclToBeMade$CellDive_ID[x],
                                  FOV_number == exclToBeMade$FOV_number[x], 
                                  grepl(getClassifierPattern(m,delim=","),PositiveMarkers)) %>%
                           pull(UUID)
                if(length(exUUIDs) == 0){ next }
                log_info(paste0("    excluding marker ", m, " from CellDive_ID ",exclToBeMade$CellDive_ID[x]," FOV ",exclToBeMade$FOV_number[x]))
                dat$PositiveMarkers[dat$UUID %in% exUUIDs] <- sapply(dat$PositiveMarkers[dat$UUID %in% exUUIDs],
                                                                    function(x){
                                                                        markers <- unlist(strsplit(x,","))
                                                                        paste(markers[markers != m], collapse=",")
                                                                    })
            }
        }
    }

    ## check that all exclusions were made properly
    exclChk <- apply(exclToBeMade, 1, function(x){
                                     etbm <- dat %>% 
                                             filter(CellDive_ID == x["CellDive_ID"],
                                                    FOV_number == as.numeric(x["FOV_number"]), 
                                                    grepl(gsub(",", "|", gsub(" ","", x["Marker_exclusion"])), PositiveMarkers))
                                       if(nrow(etbm) > 0){
                                           return(tibble(CellDive_ID = x["CellDive_ID"], 
                                                         FOV_number = as.numeric(x["FOV_number"]), 
                                                         Marker_exclusion = x["Marker_exclusion"]))
                                       }
                                   })
    exclChk <- bind_rows(exclChk)
    if(nrow(exclChk) > 0){
        log_fatal("Something went wrong. Not all marker exclusions were properly made.")
        stop("Something went wrong. Not all marker exclusions were properly made.")
    }

    return(dat)
}


#' Remove FOVs marked for exclusion in study annotation
#'
#' From marker-level data (one row per marker per cell), filter out entire 
#' FOVs marked for exclusion in parsed study annotation
#'
#' @param dat      Halo data mega table that still includes a single row for marker
#'                 in every cell
#' @param sampAnn  table of sample annotation including columns for CellDive_ID, FOV_number
#'                 and FOV_exclusion
#' @return data table with entire FOVs removed as specified in study annotation
excludeMarkerLevelFOVs <- function(dat, sampAnn){

    totalCells <- length(unique(dat$UUID))

    fovExcl <- sampAnn %>% 
               select(Sample = CellDive_ID, SPOT = FOV_number, FOV_exclusion) %>%
               filter(FOV_exclusion == "X")

    for(x in 1:nrow(fovExcl)){
        log_debug(paste0("Excluding FOV ", fovExcl$SPOT[x], " from ", fovExcl$Sample[x]))
    }

    filtDat <- dat %>%
               left_join(fovExcl, by=c("Sample", "SPOT")) %>%
               filter(is.na(FOV_exclusion)) %>%
               select(-FOV_exclusion)

    cellsRemoved <- totalCells - length(unique(filtDat$UUID))
    log_info(paste0("Excluding ", cellsRemoved, " cells due to TOTAL FOV exclusions"))

    filtDat
                  
}


#excludeMarkers <- function(dat, sampAnn){

#    totalCells <- length(unique(dat$UUID))

#    markerExcl <- sampAnn %>% 
#                  select(CellDive_ID, FOV_number, Marker_exclusion) %>%
#                  filter(!is.na(Marker_exclusion))
#    markerExcl$Marker_exclusion <- gsub(" ","",markerExcl$Marker_exclusion)

#    for(x in 1:nrow(markerExcl)){
#        for(m in unlist(strsplit(markerExcl$Marker_exclusion[x], ","))){
#            log_debug(paste0("Excluding MARKER ", m, " from FOV ", markerExcl$FOV_number[x], " in ", markerExcl$CellDive_ID[x]))
#        }
#    }

#    maxNumMarkers <- lapply(strsplit(markerExcl$Marker_exclusion, ","), function(x){ length(x) }) %>% 
#                     unlist %>% 
#                     max    
 
#    expCols <- paste0("me", seq(1:maxNumMarkers))   
#    markerExcl <- markerExcl %>%
#                  separate(Marker_exclusion, into = expCols, sep=",") %>%
#                  gather(3:(2+length(expCols)), key="TMP", value="Marker") %>%
#                  select(-TMP) %>%
#                  mutate(Excl = "X")

#    filtDat <- dat %>%
#               left_join(markerExcl, by=c("CellDive_ID", "FOV_number", "Marker")) %>%
#               filter(is.na(Excl)) %>%
#               select(-Excl)

#    cellsRemoved <- totalCells - length(unique(filtDat$UUID))
#    log_info(paste0("Excluding ", cellsRemoved, " cells due to Marker exclusions (Should be zero!!)"))

#    filtDat
#}

