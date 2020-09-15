options(dplyr.summarise.inform = FALSE)


loadRawNeighborhoodData <- function(nbhdDirs, nbhdAnalyses, fovs, nameMap = NULL){
    
    centers <- unique(c(nbhdAnalyses$nfracs$`Center Population A`,
                        nbhdAnalyses$nfracs$`Center Population B`,
                        nbhdAnalyses$navgcounts$Center))

    mainCTs <- unique(unlist(lapply(strsplit(centers,","), function(x){ x[1] })))

    oldCenters <- centers
    if(!is.null(nameMap)){
        for(nm in names(nameMap)){
            oldCenters <- gsub(paste0("^",nm), nameMap[[nm]], oldCenters)
        }
    }

    nbhdFiles <- lapply(nbhdDirs, function(x) file.path(x, dir(x))) %>% unlist

    nbhdDat <- tibble()
    for(ct in mainCTs){
        oldC <- ifelse(ct %in% names(nameMap), nameMap[[ct]], ct)

        log_debug("**************************************")
        log_debug(ct)
        log_debug("**************************************")

        nFile <- nbhdFiles[grepl(gsub("/","_",paste0("\\-",oldC,"____")),nbhdFiles)]
        log_info(paste0("Loading unformatted neighborhood counts from file: ", nFile))
        nbhdDat <- nbhdDat %>% 
                   bind_rows(readRDS(nFile) %>% filter(FOV %in% fovs) %>% select(-any_of(c("C2", "C3")))) %>%
                   unique()
    }

    nbhdDat 

}


filterForNeighborhoodCenterCells <- function(rawNbhds, center, centerCol = "C.Classifiers",
                                             markerCol = "C.FullPosMarkerStr", nameMap = NULL){

    cls    <- gsub(",.*", "", ct)
    oclass <- ifelse(cls %in% names(nameMap), nameMap[[cls]], cls)
    mrkrs  <- unlist(strsplit(ct, ","))[-1]

    ctDat <- rawNbhds %>%
             filter(grepl(getClassifierPattern(oclass), !!as.name(centerCol)))

    if(length(mrkrs) > 0){
        pos <- mrkrs[!grepl("\\-", mrkrs)]

        neg <- gsub("\\-", "", mrkrs[grepl("\\-", mrkrs)])
 
        if(length(neg) > 0){
            ctDat <- ctDat %>% 
                     filter(!grepl(paste(getClassifierPattern(neg, delim = ","), collapse = "|"), 
                                   !!as.name(markerCol)))
        }
        for(m in pos){
            ctDat <- ctDat %>%
                     filter(grepl(getClassifierPattern(m, delim = ","), !!as.name(markerCol)))
        }           
    }

    ctDat

}

getNeighborhoodCounts <- function(nbhdDirs, nbhdAnalyses, nbhdCountsDir, annCells, 
                                   fovs, markers){

    nameMap <- c("MHCIIpos_macro" = "M1",
                 "MHCIIneg_macro" = "M2")
 
    centers <- unique(c(nbhdAnalyses$nfracs$`Center Population A`,
                        nbhdAnalyses$nfracs$`Center Population B`,
                        nbhdAnalyses$navgcounts$Center))

    nbhds <- c(nbhdAnalyses$navgcounts %>% 
                 filter(Center %in% centers) %>% 
                 pull(Neighborhood),
               nbhdAnalyses$nfracs %>%
                 filter(`Center Population A` %in% centers) %>%
                 pull(`Neighborhood Population A`),
               nbhdAnalyses$nfracs %>%
                 filter(`Center Population B` %in% centers) %>%
                 pull(`Neighborhood Population B`)) %>% 
             unique()


    rawNbhds <- loadRawNeighborhoodData(nbhdDirs, nbhdAnalyses, fovs, nameMap = nameMap)

    ncounts <- tibble()
    for(ct in centers){
        log_debug(paste0("Center: ", ct))
        centerDat <- filterForNeighborhoodCenterCells(rawNbhds, ct, 
                                                      centerCol = "C.Classifiers", 
                                                      markerCol = "C.FullPosMarkerStr", 
                                                      nameMap = nameMap)                     
        
        for(nbhd in nbhds){
            log_debug(paste0("    ", nbhd))
            oldNbhd <- nbhd
            if(gsub(",.*", "", oldNbhd) %in% names(nameMap)){
                oldNbhd <- paste0(nameMap[[gsub(",.*", "", oldNbhd)]], 
                                  gsub("^.*,", ",", oldNbhd))      
            }

            ncounts <- ncounts %>%
                       bind_rows(centerDat %>%
                                 filterForNeighborhood(oldNbhd, markers) %>%
                                 mutate(CenterCellType = ct, NeighborhoodCellType = nbhd) %>%
                                 unique() %>%
                                 group_by(FOV, C.UUID, CenterCellType, NeighborhoodCellType) %>%
                                 summarize(N.Count = n()))
        }                
    }

    cellDat <- annCells %>%
               select(any_of(c("CellDive_ID", "Sample_ID", "FOV_ID", "Band", "UUID", "PositiveMarkers"))) %>%
               rename(C.UUID = UUID, C.PositiveMarkers = PositiveMarkers)

    ncounts %>% 
    addZeroNeighborhoodCounts() %>%
    left_join(cellDat, by = intersect(names(.), names(cellDat)))
}


addZeroNeighborhoodCounts <- function(ncounts){

       
    fovs <- unique(ncounts$FOV)
    nCellTypes <- unique(ncounts$NeighborhoodCellType)

    ## fill in missing neighborhood cell types in all FOVs
    allNcounts <- lapply(fovs, function(fov){ 

        nCts <- ncounts %>%
                filter(FOV == fov) %>%
                spread(NeighborhoodCellType, N.Count, fill=0)

        for(nct in setdiff(nCellTypes, names(nCts))){
            nCts[[nct]] <- 0
        }       

        nCts %>%
        gather(all_of(nCellTypes), 
               key="NeighborhoodCellType",
               value="N.Count") %>%
        rename(FOV_ID = FOV)
        
    }) %>%
    bind_rows()

    errs <- allNcounts %>%
            group_by(FOV_ID, C.UUID, CenterCellType) %>%
            summarize(Count = n()) %>%
            filter(Count != length(nCellTypes))

    if(nrow(errs) > 0){
        stop("Center cells do not all have the same number of neighborhood cell types!!")
    }

    allNcounts
}


#### COUNTS WILL BE LOADED/GENERATED FOR ALL CELLDIVE IDS IN cellDat
loadNeighborhoodCounts <- function(nbhdDirs, nbhdCountsDir, cellDat, markers){ 

    countsFiles <- file.path(nbhdCountsDir, 
                             paste0(unique(cellDat$CellDive_ID), "_nbhd_counts.rda"))

    ncounts <- lapply(countsFiles, function(x){
                   if(fileDone(x)){
                       readRDS(x)
                   } #else {
#                       getNeighborhoodCounts(nbhdDirs, 
#                                             nbhdAnalyses,  
#                                             nbhdCountsDir,  
#                                             cellDat, 
#                                             unique(cellDat$FOV_ID),  
#                                             markers)
#                   }
               }) %>%
               bind_rows()

    ncounts
}

        
     
    
#for(cdid in unique(fin$CellDive_ID)){
#    tmp <- fin %>%
#           filter(CellDive_ID == cdid) %>%
#           select(CellDive_ID, Sample_ID, FOV_ID, C.UUID, everything())

#    saveRDS(tmp, file.path(nbhdCountsDir, paste0(cdid, "_nbhd_counts.rda")))
#}


