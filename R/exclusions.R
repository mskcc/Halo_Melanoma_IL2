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


getFOVexclusions <- function(dat, sampAnn){
    dat %>%
    left_join(sampAnn %>% select(CellDive_ID, FOV_number, FOV_exclusion),  
              by = c("CellDive_ID", "FOV_number"))
}

getMarkerExclusions <- function(dat, sampAnn){
    cols <- sapply(sampAnn$Marker_exclusion, function(x) length(unlist(strsplit(x, ","))) ) %>% unlist() %>% max()

    mExcl <- sampAnn %>% 
             select(CellDive_ID, FOV_number, Marker_exclusion) %>%
             separate(Marker_exclusion, paste0("M_", seq(cols))) %>%
             gather(paste0("M_", seq(cols)), key="TMP", value="Marker") %>%
             select(-TMP) %>%
             filter(!is.na(Marker), Marker != "") %>%
             mutate(Marker = trimws(Marker), Marker_exclusion = "X")

    dat %>%
    left_join(mExcl, by = c("CellDive_ID", "FOV_number", "Marker")) 
}

#' Get indices of cells to be excluded because the percentage of drift/loss is
#' greater than the set threshold
#'
#' Get indices of cells to be excluded because the percentage of drift/loss is
#' greater than the set threshold
#' 
#' @param drift     table of drift/loss percentages
#' @param dat       cell data
#' @param threshold maximum percent of drift/loss allowed
#' @return vector of dat indices that need to be marked for exclusion
getDriftExclusions <- function(drift, dat, threshold){

    excl <- c()
    drft <- drift %>% select(Image_Location=image_location, XMin=x_min, XMax=x_max,
                             YMin=y_min, YMax=y_max, drift_loss_pixel_pct)
    tmp <- full_join(dat, drft, by=c("Image_Location","XMin","XMax","YMin","YMax"))
    excl <- which(!is.na(tmp$drift_loss_pixel_pct) & tmp$drift_loss_pixel_pct > threshold) 
    excl
}

#' Get indices of cell to be excluded because they lie outside the border padding 
#' 
#' Get indices of cell to be excluded because they lie outside the border padding 
#' 
#' @param fov        fov to be marked for exclusions
#' @param samp       sample to be marked
#' @param dat        data tibble containing data to be marked
#' @param borderPad  number in pixels to trim from data
getBorderPaddingExclusions <- function(fov, samp, dat, borderPad){

    fovDat <- dat %>% filter(Sample==samp, SPOT==fov)
    
    #bbFOV <- list(X0=1,Y0=-3375,X1=5363,Y1=1)
    bb <- list(X0 = bbFOV$X0 + borderPad,
               X1 = bbFOV$X1 - borderPad,
               Y0 = bbFOV$Y0 + borderPad,
               Y1 = bbFOV$Y1 - borderPad)

    which( dat$Sample == samp & dat$SPOT == fov &
         ( dat$X < bb$X0 | dat$X > bb$X1 | dat$Y < bb$Y0 | dat$Y > bb$Y1)) 
    
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
#' @param sampAnn         FOV annotations for one sample
#' @param borderPad       number in microns indicating the minimum distance between a cell and the FOV
#'                        border in order for that cell NOT to be excluded
#' @param driftThreshold  maximum pixel percent determined to have drifted in order to NOT be excluded
#' @return tibble matching {dat} exactly, with an extra column, EXCLUDE, added to indicate which
#'         cells should not be analyzed
#' @export
markExclusions <- function(dat, drift, fovAnn, driftThreshold=0.1){

    borderPad_px <- borderPad/pixel2um

    exFiles <- epFiles <- gFiles <- NULL

    log_debug("Converting min/max coordinates to midpoint coordinates")
    dd <- dat %>% mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) ## need to keep X|Y min|max for drift loss 
    dd$EXCLUDE <- ""

    log_info(paste("\t","Sample","FOV","Reason","Num.Excl.Cells",sep="\t"))


    ## get exclusions determined by lab personel, found in FOV annotations meta data spreadsheet   
    dd <- dd %>%
          getFOVexclusions(sampAnn) %>%
          getMarkerExclusions(sampAnn)
    

    ######### START HERE ##########


    ## get drift/loss exclusions
    if(!is.null(drift)){
        driftExcl <- getDriftExclusions(drift, dd, driftThreshold)
        if(length(driftExcl) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, driftExcl, paste0("DRIFT_",driftThreshold*100))
            tmp <- dd %>% 
               filter(grepl("DRIFT",EXCLUDE)) %>% 
               group_by(Sample, SPOT, EXCLUDE) %>%
               summarize(`Num.Excl.Cells`=n())
            for(row in 1:nrow(tmp)){
               log_info(paste("\t",paste(tmp[row,], collapse="\t")))
            }
        } else {
            log_info(paste("\t",samp,"ALL",paste0("DRIFT_",driftThreshold*100),0,sep="\t"))
        }
    } else {
        log_info(paste("\t",samp,"ALL",paste0("DRIFT_",driftThreshold*100),0,sep="\t"))
    }


    ## mark cells that fall inside Halo exclusion boundaries (exclusions, epidermis, glass) and also
    ## those that fall outside limits set by border pad 
    for(fov in unique(dd$SPOT)){ 
        aFile <- epFile <- gFile <- iFile <- NULL

        ## exclude points that fall within outside border
        borderExcl <- getBorderPaddingExclusions(fov, samp, dd, borderPad_px)
        if(length(borderExcl) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, borderExcl, paste0("PADDED_",borderPad,"_um"))
        }
        log_info(paste("\t", samp, fov, paste0("PADDED_",borderPad,"_um"), length(borderExcl), sep="\t"))

        ## if parsed halo annotation is not given, get all halo annotations files for this FOV
        if(is.null(haloAnn)){
            sampAnns <- aFiles[grep(paste0(samp,"_"),aFiles)]
            aFile <- sampAnns[grep(paste0("Spot",fov,".annotations"),sampAnns)]
            sampAnns <- epFiles[grep(paste0(samp,"_"),epFiles)]
            epFile <- sampAnns[grep(paste0("Spot",fov,".annotations"),sampAnns)]
            sampAnns <- gFiles[grep(paste0(samp,"_"),gFiles)]
            gFile <- sampAnns[grep(paste0("Spot",fov,".annotations"),sampAnns)]
            sampAnns <- iFiles[grep(paste0(samp,"-"),iFiles)]
            iFile <- sampAnns[grep(paste0("Spot",fov,".annotations"),sampAnns)]
        }

        ## get halo exclusions 
        haloExcl <- getHaloExclusions(fov, dd, "Exc", haloAnnotations=haloAnn, aFile=aFile, 
                                      boundaryColors=boundaryColors, boundaryReassignmentFile=boundaryReassignmentFile)
        if(length(haloExcl) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, haloExcl, "HALOExclusion")
        }
        log_info(paste("\t", samp, fov, "HALOExclusion", length(haloExcl), sep="\t"))

        ## get halo epidermis exclusions 
        haloEpi <- getHaloExclusions(fov, dd, "Epi", haloAnnotations=haloAnn, aFile=epFile, 
                                     boundaryColors=boundaryColors, boundaryReassignmentFile=boundaryReassignmentFile)
        if(length(haloEpi) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, haloEpi, "HALOEpidermis")
        }
        log_info(paste("\t", samp, fov, "HALOEpidermis", length(haloEpi), sep="\t"))

        ## get halo glass exclusions 
        haloGls <- getHaloExclusions(fov, dd, "Gls", haloAnnotations=haloAnn, aFile=gFile, 
                                     boundaryColors=boundaryColors, boundaryReassignmentFile=boundaryReassignmentFile)
        if(length(haloGls) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, haloGls, "HALOGlass")
        }
        log_info(paste("\t", samp, fov, "HALOGlass", length(haloGls), sep="\t"))

        ## debug with plots
        if(printPlots){
            log_info("Printing plots for debugging")
            plotExclusions(dd[which(dd$SPOT==fov),],haloAnnotations=haloAnn,iFiles,aFiles,epFiles,gFiles,outDir=debugDir,
                           boundaryReassignmentFile=boundaryReassignmentFile)
        }
    }
    dd %>% select(-(X:Y))
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

#' Remove FOVs marked for exclusion in study annotation
#' 
#' From cell-level data (one row per cell), filter out entire FOVs marked for
#' exclusion in study annotation.
#' 
#' @param dat      cell level data table including columns for Sample_ID, FOV_ID,
#' @param sampAnn  a table of sample annotation including columns for Sample_ID,
#'                 FOV_ID and FOV_exclusion
#excludeCellLevelFOVs <- function(dat, sampAnn){

#    totalCells <- length(unique(dat$UUID))

#    fovExcl <- sampAnn %>%
#               select(CellDive_ID
    
#}


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

