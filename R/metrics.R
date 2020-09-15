tableComplete <- function(dat, uniqueGrpCols, colToCheck){
    res <- dat %>%
           select_at(c(uniqueGrpCols, colToCheck)) %>%
           unique() %>%
           group_by_at(uniqueGrpCols) %>%
           summarize(Count = n()) %>%
           pull(Count) %>%
           unique()

    if(length(res) > 1 || res != nrow(dat %>% select(colToCheck) %>% unique())){
        return(FALSE)
    }
    return(TRUE)
}


areaBB<-function(bb){
    (bb$X1-bb$X0)*(bb$Y1-bb$Y0)
}

pointsInPolygon <- function(pts,poly){
    point.in.polygon(pts$X,pts$Y,poly$X,poly$Y)==1
}

pointsInsidePolygonList <- function(pts,bList) {

    inside=bList %>%
        map(function(bn){point.in.polygon(pts$X,pts$Y,bn$X,bn$Y)>0}) %>%
        as.data.frame %>%
        apply(.,1,any)

    inside

}

#' Generate data frame of random points (is it random???)
#' 
#' Generate a data frame of random points that fall within the plot area
#' 
#' @param   nGrid   number of points to generate
#' @param   bbG     a list containing X0,X1,Y0,Y1 representing the boundary
#'                  of the plot area
#' @return data frame of X and Y values of the random points
generateSobolGrid <- function(nGrid,bbG) {
    gg=sobol(nGrid,2,scrambling=2)
    data.frame(
        X=gg[,1]*(bbG$X1-bbG$X0)+bbG$X0,
        Y=gg[,2]*(bbG$Y1-bbG$Y0)+bbG$Y0
        )
}

#' Add columns for band area to existing tibble with a column of bands
#'
#' Add columns for band area to existing tibble with a column of bands
#'
#' @param bdat      tibble of band data, including at least a column named Band
#' @param bandArea  a data frame or tibble with an area value for each band in bdat
#' @return tibble containing all data in bdat and bandArea
joinBandArea <- function(bdat,bandArea) {

    bdat$Band <- as.character(bdat$Band)
    log_debug("adding band area to band counts\n")
    xx        <- full_join(bdat,bandArea)
    log_debug("filling in NAs with 0s\n")
    xx$n[is.na(xx$n)] <- 0
    xx$Band <- factor(xx$Band,levels=bandArea$Band)
    log_debug("arranging by Band\n")
    xx %<>% arrange(Band)
    xx

}

#' Remove points from data set that fall inside exclusion boundaries
#' 
#' Read Halo boundaries annotations file and remove from data set any cells
#' that fall inside those regions
#' 
#' @param ds     tibble containing at minimum, Sample, SPOT, UUID, X, Y, Marker, Value
#' @param excB   a list of excluded region tables as returned by cleanBoundaries;
#'               if NULL, must provide the original Halo annotations file in XML format
#' @param aFile  Halo boundaries annotation file in XML format; if NULL, must provide
#'               the list of excluded region tables as returned by cleanBoundaries()
#' @return tibble filtered to remove points that fall inside exclusion boundaries
removeExcludedPoints <- function(ds,aFile=NULL,excB=NULL) {

    if(!all(c('X','Y') %in% colnames(ds))){
        ds <- convertCellMinMaxToMidpoints(ds)
    }

    if(is.null(excB)){
        boundaries <- readHaloAnnotations(aFile)
        allBoundaries <- cleanBoundaries(boundaries)
        excB <- c(allBoundaries$excB,allBoundaries$glsB,allBoundaries$epiB)
    }

    spCells <- ds %>% dplyr::select(UUID,X,Y)
    coordinates(spCells) <- ~X+Y

    if(length(excB)>0){
        spExcB <- seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})
        excludedCellIdx <- NULL
        for(jj in length(spExcB)) {
            excCellJJ <- unlist(over(spExcB[[jj]],geometry(spCells),returnList=T))
            excludedCellIdx <- union(excludedCellIdx,excCellJJ)
        }
        if(length(excludedCellIdx)>0){
            ds <- ds[-excludedCellIdx,]
        }
    }

    return(ds)

}

#' Calculate total area in a single FOV
#'
#' Calculate total area in a single FOV
#' 
#' @param dat              tibble of all halo data for a single FOV
#' @param boundaries       list of pre-parsed halo annotations for single FOV
#' @param maxG             ???? default=5
#' @return total FOV area in pixels/mm^2 (????) 
#' @export
calculateAreaTotalFOV <- function(fovDat, boundaries, maxG=5){

    s          <- unique(fovDat$CellDive_ID)
    excB       <- c(boundaries$excB, boundaries$epiB, boundaries$glsB)
    bb         <- getConstantBoundingBox(s)
    gg         <- generateSobolGrid(10^maxG, bb)
    if(length(excB) > 0){
        excludedPts <- pointsInsidePolygonList(gg, excB)
        gg <- gg[!excludedPts,]
    }

    p2tomm * areaBB(bb) * nrow(gg)/10^maxG

}

#' Calculate all FOV areas for a single sample
#'
#' Given a complete table of all cell coordinates for a single 
#' sample, represented by a CellDive_ID, and all exclusion
#' boundary coordinates, calculate the total tissue area
#' for each FOV
#' 
#' @param cdid        CellDive_ID of sample for which FOV areas will
#'                    be calculated
#' @param dataFile    RDA file containing all sample cell coordinates
#' @param boundaries  nested list of boundary coordinates including
#'                    all exclusion boundaries
#' @param maxG
#' 
#' @return tibble with one row per FOV and FOVArea for each
getFOVAreas <- function(cdid, dataFile, boundaries, maxG = 5){ 
    sdat <- loadHaloDataFile(dataFile, filterExclusions = TRUE)
    areas <- tibble()
    for(fov in unique(sdat %>% filter(CellDive_ID == cdid) %>% pull(FOV_number))){
        fdat <- sdat %>% filter(FOV_number == fov)
        fbounds <- boundaries[[cdid]][[as.character(fov)]]    
        areas <- areas %>% bind_rows(tibble(CellDive_ID = cdid,    
                                            FOV_number = fov,
                                            FOVArea = calculateAreaTotalFOV(fdat, fbounds, maxG=maxG)))
    }
    return(areas)
}


#' Compute band areas for each interface bin
#' 
#' Given a list of distances from a tumor boundary, or 'interface bins', calculate area
#' for each bin using a sobol grid
#' 
#' @param tumorBoundariesList   a list where each element is a tibble containing tumor boundaries
#'                              returned from readHaloAnnotations()
#' @param bb                    a list containing X0,X1,Y0,Y1 representing the boundary
#'                              of the trimmed FOV
#' @param interfaceBins         vector of distances (integers) that define each band; default = (-20:20)*10 
#' @param nGrid                 number of random points to use for area calculations; default=100
#' @param maxSeg                maximum segments to divide long boundary segments into; default=1000
#' @param exclusionsB           a list where each element is a tibble containing exclusion boundaries
#'                              returned from readHaloAnnotations()
computeBandAreasMC <- function(tumorBoundariesList,bb,interfaceBins,nGrid=100,
                               maxSeg=1000,exclusionB) {

    ## generate random points all across the plot
    log_debug("generating sobol grid\n")
    gg <- generateSobolGrid(nGrid,bb)
    if(length(exclusionB)>0) {
        log_debug("There are %s exclusion boundaries. Removing excluded points",length(exclusionB),"\n")
        excludedPts <- pointsInsidePolygonList(gg,exclusionB)
        gg          <- gg[!excludedPts,]
    }

    log_debug("trimming boundaries\n")
    tumBTrim <- tumorBoundariesList %>%
                map(trimDFToBB,bb) %>%
                map(subDivideLongSegments,maxSeg)

    log_debug("getting X,Y points for tumor boundaries\n")
    tumorBoundariesMergeXY <- tumBTrim %>% bind_rows %>% dplyr::select(X,Y) %>% as.data.frame

    log_debug("finding distance from points to interface\n")
    gg$Z <- findDistanceFromPointsToInterfacePoint(gg,tumorBoundariesMergeXY)

    log_debug("assign negative distance to points inside tumor\n")
    insideTumor       <- pointsInsidePolygonList(gg,tumorBoundariesList)
    gg$Z[insideTumor] <- -gg$Z[insideTumor]

    log_debug("getting area table\n")
    areaTable <- table(cut(gg$Z*pixel2um,interfaceBins))

    log_debug("making bandAreas dataframe\n")
    bandAreas <- as.data.frame.table(p2tomm*areaBB(bb)*areaTable/nGrid)
    colnames(bandAreas) <- c("Band","Area")
    return(bandAreas)
}

#' Get counts of cells that fall within given distances of
#' tumor boundaries for a single FOV
#'
#' Given a set of interface bands, or set distances from a tumor 
#' boundary, count the numbers of cells that fall within those bands
#' 
#' @param  ds             a tibble containing data for one FOV
#' @param bb              a list containing X0,X1,Y0,Y1 representing the boundary
#'                        of the trimmed FOV
#' @param allBoundaries   a list returned from readHaloAnnotations
#' @param interfaceBins   vector of distances (integers) that define each band; default = (-20:20)*10 
#' @return tibble
#' @export
getBandAssignments <- function(ds,bb,allBoundaries,interfaceBins=(-36:36)*10){

    log_debug("getting boundaries\n")
    tumB <- allBoundaries$tumB
    if(is.null(tumB) || length(tumB) == 0){
        log_info("No TUMOR boundaries found in file %s. Skipping.")
        return()
    }

    log_debug("getting midpoint of each cell\n")
    ## get midpoint of each cell
    if(!all(c("X","Y") %in% names(ds))){
        ds <- ds %>%
              mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2)
    }

    ## create a SpatialPointsDataFrame
    spCells <- ds %>% dplyr::select(UUID,X,Y)
    coordinates(spCells) <- ~X+Y

    tumorContour <- tumB %>% map(trimDFToBB,bb)
    log_debug("pulling out X,Y coords\n")
    ## get x and y coords of all remaining tumor cells
    tumorBoundariesMergeXY <- tumorContour %>%
                              bind_rows %>%
                              dplyr::select(X,Y) %>%
                              as.data.frame

    log_debug("finding distances closest to boundary\n")
    ## for each point, find distance to closes boundary point
    dsPt   <- ds %>% dplyr::select(X,Y) %>% as.data.frame
    dTumor <- findDistanceFromPointsToInterfacePoint(
                  dsPt,
                  tumorBoundariesMergeXY
               )
    ## list of true/false indicating whether point is in a tumor
    tumorPoints <- pointsInsidePolygonList(dsPt,tumB)
    dTumor      <- ifelse(tumorPoints,-dTumor,dTumor)

    ## convert to um
    dTumor      <- dTumor*pixel2um
    ds$Distance <- dTumor
    ds$Band     <- cut(dTumor,interfaceBins) ## assign each point to a bin

    return(ds)

}


#' Calculate area of interface bands, each defined as the collection of points that are located
#' between X and Y microns of the tumor boundary [ TO DO: REWORD THIS? ]
#' 
#' Calculate area of each band, defined as the collection of cells that are located
#' between X and Y microns of the tumor boundary [ TO DO: REWORD THIS? ]
#' **Do this for 10^maxG and for 10^maxG-1
#' 
#' @param allBoundaries  halo boundaries after removal of contained boundaries (list returned by 
#'                       cleanBoundaries())
#' @param maxG           maximum factor of 10 to use for generating random points on grid [ REWORD??? ]
#' @param bb             a list containing X0,X1,Y0,Y1 representing the boundary
#'                       of the trimmed FOV
#' @param interfaceBins  vector of distances (integers) that define each band; default = (-20:20)*10
#' @return a dataframe with a row for each interface bin and a column for each area calculation (one for 
#'         10^maxG and one for 10^maxG-1)
calculateBandAreas <- function(allBoundaries,maxG=5,bb,interfaceBins) {

    if(is.null(interfaceBins)){
        interfaceBins <- (-20:20)*10
    }
    tumB <- allBoundaries$tumB
    excB <- c(allBoundaries$excB,allBoundaries$epiB,allBoundaries$glsB)

    bandArea <- list()
    for(nGrid in 10^((maxG-1):maxG)) {
        maxSeg <- max(unlist(bb))/sqrt(10*nGrid)
        log_debug("nGrid = %s  maxSeg = %s ",nGrid,maxSeg,"\n")
        bandArea[[as.character(nGrid)]] <- computeBandAreasMC(tumB,bb,interfaceBins,
                                                               nGrid,maxSeg,excB)
    }

    aa           <- bind_cols(bandArea) %>%
                    dplyr::select(dplyr::matches("Area"))
    rownames(aa) <- bandArea[[1]][,1]

    aa

}


#' Calculate area of each X-micron band around a tumor boundary 
#' 
#' Calculate area of each X-micron band around a tumor boundary
#'
#' @param dat                 Halo data for one or more samples, loaded from *.rda file  
#' @param haloAnnotations     list of halo boundary annotations
#' @param maxG                the maximum factor of 10 to use for generating random points for area calculation
#' @param interfaceBins       a vector of distances from tumor interface in which cells will be binned; 
#'                            default=(-36:36)*10
#' @export
calculateInterfaceArea <- function(dat, haloAnnotations, maxG=5,
                                   interfaceBins=(-36:36)*10,
                                   numThreads=6){

    cl <- makeCluster(numThreads, type = "FORK", outfile = "")
    clusterExport(cl,
                  c("dat", "haloAnnotations", "maxG", "interfaceBins"),
                  envir = environment())

    allAreas <- parLapply(cl, unique(dat$CellDive_ID), function(s){
        fovs <- dat %>% filter(CellDive_ID == s) %>% pull(FOV_number) %>% unique()
        areas <- list()
        for(fov in fovs){
            log_debug(paste0("getting tumor boundary interval areas for\tSample ",s, "\tFOV ",fov))
            boundaries <- haloAnnotations[[s]][[as.character(fov)]]
            if(is.null(boundaries) || length(boundaries) == 0 || length(boundaries$tumB) == 0){
                log_debug("no tumor boundaries found. skipping\n")
                next
            }
            boundaries$tumB <- boundaries$tumB[!is.na(names(boundaries$tumB))]
            bb <- getConstantBoundingBox(s)
            ba <- calculateBandAreas(boundaries, maxG=maxG, bb, interfaceBins)
            areas[[length(areas) + 1]] <- tibble(CellDive_ID = s,
                                                 FOV_number = fov,
                                                 Band = rownames(ba),
                                                 Area = ba[,1])
        }
        if(length(areas) > 0){ areas <- bind_rows(areas) }
        areas
    })
    stopCluster(cl)

    if(length(allAreas) == 0){ return(NULL) }
    bind_rows(allAreas)

}


#' Calculate distances from each cell to nearest point on tumor boundary 
#' 
#' Calculate distances from each cell to nearest point on tumor boundary and assign
#' each cell to an interface band/bin
#'
#' @param dat                 Halo data for one sample or more samples, loaded from *.rda file(s)  
#' @param haloAnnotations     list of pre-parsed halo boundary annotations
#' @param interfaceBins       a vector of distances from tumor interface in which cells will be binned; 
#'                            default=(-36:36)*10
#' @return table of all cells including distance to nearest tumor boundary and band assignments 
#' @export
addInterfaceBandInfo <- function(dat, haloAnnotations=NULL,
                                 interfaceBins=(-36:36)*10, numThreads=6){

    cl <- makeCluster(numThreads, type="FORK", outfile="")
    clusterExport(cl, c("dat","haloAnnotations","interfaceBins"), envir=environment())
    allBa <- parLapply(cl, unique(dat$CellDive_ID), function(s){
        dd <- dat %>% filter(CellDive_ID == s)
        bandAssignments <- list()
        for(fov in unique(dd$FOV_number)){
            log_debug(paste0("Getting band assignments for sample ",s, " FOV ",fov))
            ds <- dd %>% filter(FOV_number == fov)
            boundaries <- haloAnnotations[[s]][[as.character(fov)]]
            if(is.null(boundaries) || length(boundaries) == 0 || length(boundaries$tumB) == 0){
                next
            }
            boundaries$tumB <- boundaries$tumB[!is.na(names(boundaries$tumB))]
            bb <- getConstantBoundingBox(s)
            ba <- getBandAssignments(ds, bb, boundaries, interfaceBins)
            bandAssignments[[length(bandAssignments) + 1]] <- ba
        }
        bind_rows(bandAssignments)
    }) %>%
    bind_rows()
    stopCluster(cl)

    return(allBa)
}


#' Get full path to appropriate metrics directory depending on cell region 
#' of interest
#' 
#' Get full path to appropriate metrics directory depending on cell region 
#' of interest
#'
#' @param cellRegion   cell region of interest [fov|interface inside|interface outside]
#' @param metricsDir   root metrics directory
#'
#' @return full path to directory containing the metrics of a specific cell region
getMetricsSubDir <- function(cellRegion, metricsDir){
    metrics_subdir <- switch(cellRegion,
                             "fov" = "fovs",
                             "interface inside" = "inside_interface",
                             "interface outside" = "outside_interface",
                             "interface" = "infiltration",
                             "neighborhood" = "tumor_neighborhood",
                             NULL)
    if(is.null(metrics_subdir)){ return(NULL) }
    file.path(metricsDir,metrics_subdir)
}

#' Calculate and save all metrics for one cell region and one calculation unit
#'
#' Precalculate all metrics specified in a list of analyses (fractions, densities, etc.)
#' focusing on a single cell region and one type of calculation unit. A single calculation
#' will be made for each calculation type (fracs, dens, etc.) and each unique calculation
#' unit (e.g., FOV, Sample, FOV+Band, etc.). Results will be save to the specified directory.
#' 
#' @param annCells         tibble of annotated cell data where each row represents a single 
#'                         unique cell
#' @param sampAnn          tibble of flattened study sample annotation
#' @param markers          vector of all markers used in study
#' @param nbhdCounts       tibble of macrophage neighborhood counts
#' @param tumorNbhdCells   vector of UUIDs of cells that are within 30 microns of at least
#'                         one tumor cell
#' @param analysisList     list of tibbles where each table describes all analyses to run
#'                         for a single calculation type (fractions, densities, macrophage neighborhoods)
#' @param metricsDir       directory where results will be written
#' @param cellRegion       cell region to focus on; [fov|interface|interface inside|interface outside|neighborhood]
#' @param calcUnitCols     column name(s) of calculation units
#' @param numThreads       number of threads 
#' 
#' @return nothing
precalculateMetrics <- function(annCells, sampAnn, markers, nbhdCounts, tumorNbhdCells, analysisList,
                                metricsDir, cellRegion, calcUnitCols, cellDiveID = "All", 
                                numThreads = 1){

    cr <- getCellRegion(cellRegion)
    outDir <- getMetricsSubDir(cr, metricsDir)
    calcUnit <- getCalcUnit(calcUnitCols)

    log_info(paste0("Calculating metrics for [", calcUnit, "], cell region [", cr, "]"))

    analyses <- analysisList
    tme <- NULL
    nbhd <- NULL
    areas <- NULL

    if(cr == "neighborhood"){
        analyses <- analysisList[names(analysisList) == "fractions"]
        tme      <- tumorNbhdCells
    } else {
        if(any(c("navgcounts", "nfracs") %in% names(analysisList))){
            nbhd <- nbhdCounts %>%
                    left_join(sampAnn %>% select(FOV_ID, Sample_ID), 
                              by = intersect(names(nbhdCounts), c("FOV_ID", "Sample_ID"))) %>%
                    mutate(UUID = C.UUID) %>%
                    filterDataForCellRegion(cr) %>%
                    select(-UUID) %>%
                    transformCalcUnit(calcUnitCols)
        }
        if("densities" %in% names(analysisList)){
            areas <- filterForAreas(annCells, cr, groupBy = calcUnitCols) %>%
                     transformCalcUnit(calcUnitCols)
        }
    }

    dat <- annCells %>%
           filterDataForCellRegion(cr, tumorNbhdCells = tme) %>%
           transformCalcUnit(calcUnitCols)

    popDat <- getPopulationMetrics(dat,
                                   analyses,
                                   markers,
                                   areas = areas,
                                   calcUnit = calcUnit,
                                   nbhdCounts = nbhd,
                                   outDir = outDir,
                                   numThreads = numThreads, 
                                   cellDiveID = cellDiveID)
    return(NULL)
}

#' Get numbers of cells belonging to a list of populations of interest
#' 
#' For each population of interest, count the number of cells determined to be in that population
#' 
#' @param dat           table of annotated cells where one row corresponds to a unique cell and columns
#'                      include at least Classifiers and PositiveMarkers
#' @param populations   vector of comma-delimited character strings describing the cell populations to be counted;
#'                      each element of a single delimited string may be a cell Classifier or a positive or
#'                      negative marker; a positive marker is indicated simply by the marker name and a negative
#'                      marker is indicated by the marker name with '-' appended (e.g., "PD1-")
#' @param markers       vector of all marker names, used to distinguish a marker from a classifier in a
#'                      population string
#' @param calcUnit      data column name containing values for which a single count should be made
#'                      (default: FOV_ID)
#' @param numThreads    max number of threads to use when available
#' @param outFile       path to file that contains pre-counted population data OR the file to which data
#'                      should be saved
#'
#' @return table with columns 'Population', [calcUnit], and Count
getPopulationCounts <- function(dat, populations, markers, calcUnit = "FOV_ID", numThreads = 1, outFile = NULL){

    if(fileDone(outFile)){
        log_info(paste0("Loading pre-computed counts from file ",outFile))
        return(readRDS(outFile))
    }

    popCounts <- tibble()

    log_info("No pre-computed counts file found. Generating now...")
    cl <- makeCluster(numThreads, type="FORK", outfile="")
    clusterExport(cl, c("dat","markers","calcUnit","populations"), envir=environment())
    popCounts <- parLapply(cl, 1:length(populations), function(x){
                        log_debug( paste0(" [", x ,"] ", populations[x]) )
                        filterDataForCondition(dat, populations[x], markers) %>%
                        group_by_at(calcUnit) %>%
                        summarize(Count = n()) %>%
                        mutate(Population = populations[x])
                 }) %>%
                 bind_rows()
    stopCluster(cl)

    cuStr <- paste0(paste(calcUnit, collapse = "+"),
                    ifelse(length(calcUnit) > 1, " combo", ""))
    numUnits <- popCounts %>% select_at(calcUnit) %>% unique() %>% nrow()
    log_debug(paste0("Counted populations in ", numUnits, " ", cuStr, "s"))

    ## make sure all populations are included in data
    ## do NOT use fill in spread because some NAs are valid 
    ## and should not be filled (other leuk)
    popCounts <- popCounts %>% spread(Population, Count)
    for(pop in populations[!populations %in% names(popCounts)]){
        popCounts <- popCounts %>% mutate(!!pop := NA)
    }
    popCounts <- popCounts %>% gather(populations, key = "Population", value = "Count")

    ## NAs are converted to strings when data is spread so 
    ## need to be converted back
    for(cu in calcUnit){ popCounts[[cu]][popCounts[[cu]] == "<NA>"] <- NA }

    ## Other_leuk is the only cell type that may be NA 
    ## (CD45 was not included in all samples)
    ## so set all other NA to 0.
    popCounts$Count[intersect(which(is.na(popCounts$Count)),
                      which(!grepl("Other_leuk", popCounts$Population)))] <- 0

    ## make sure there is a value for every 
    chkTbl <- popCounts %>% filter(Population != "Other_leuk") %>% ungroup()
    if(!tableComplete(chkTbl, calcUnit, "Population")){
        msg <- paste0("Some population counts missing for some ", cuStr, "s")
        log_error(msg)
        stop(msg)
    }

    if(!is.null(outFile)){
        saveRDS(popCounts, file = outFile)
    }

    popCounts

}

#' Get fractions of cell populations
#'
#' Given a tibble of cell population counts and a set of population fractions to calculate,
#' generate a table of said fractions 
#'
#' @param counts      tibble of counts for all cell populations included in list of fractions
#' @param conditions  tibble describing population fractions to calculate, including at minimum
#'                    columns Condition and Population, where Population is a character string
#'                    describing the overall cell population (denominator) and Condition is a
#'                    character string describing the subpopulation (numerator)
#' @param calcUnit    data column name containing values for which a single fraction
#'                    should be calculated (default: FOV_ID)
#' @param outFile     optional; path to RDA file where results should be saved
#'
#' @return tibble of cell population fractions
getPopulationFractions <- function(counts, conditions, calcUnit = "FOV_ID", outFile = NULL){

    if(fileDone(outFile)){
        log_info(paste0("Loading pre-computed fractions from file ",outFile))
        return(readRDS(outFile))
    }

    fracs <- tibble()

    condCounts <- counts %>% select(Condition = Population, CondCount = Count, everything())
    popCounts  <- counts %>% select(PopCount = Count, everything())

    log_info("No pre-computed fractions file found. Generating now...")
    fracs <- conditions %>%
             left_join(condCounts, by = intersect(names(conditions), names(condCounts)))
    fracs <- fracs %>%
             left_join(popCounts, by = intersect(names(fracs), names(popCounts))) %>%
             mutate(Fraction = CondCount/PopCount)

    if(!is.null(outFile)){
        saveRDS(fracs, file = outFile)
    }

    fracs
}

#' Get densities of cell populations
#' 
#' Given a tibble of cell population counts and one of areas per calculation unit, calculate
#' density for each population 
#' 
#' @param counts      tibble of counts for all cell populations for which densities are 
#'                    to be calculated
#' @param areas       tibble of areas including a single value for each calculation unit
#' @param calcUnit    data column name containing values for which a single density value
#'                    should be calculated (default: FOV_ID)
#' @param outfile     optional; path to RDA file where results should be saved
#'
#' @return tibble of cell population densities
getPopulationDensities <- function(counts, areas, calcUnit = "FOV_ID", outFile = NULL){

    if(fileDone(outFile)){
        log_info(paste0("Loading pre-computed densities from file ",outFile))
        return(readRDS(outFile))
    }

    dens <- tibble()

    log_info("No pre-computed densities file found. Generating now...")
    areasPer <- areas %>%
                group_by_at(calcUnit) %>%
                summarize(TotalArea = sum(Area, na.rm = T))

    dens <- counts %>%
            group_by_at(c(calcUnit,"Population")) %>%
            mutate(TotalCounts = sum(Count, na.rm = T)) %>%
            left_join(areasPer, by=c(calcUnit)) %>%
            mutate(Density = TotalCounts/TotalArea)

    if(!is.null(outFile)){
        saveRDS(dens, file = outFile)
    }

    dens
}

#' Get counts of neighborhood cells by type for each center cell type
#' 
#' Get full list of neighborhood cells by filtering for specific center cell types
#' including any functional markers and appending to original list
#' 
#' @param nbhdCounts   formatted tibble of counts including at minimum columns:
#'                     C.UUID - UUID of a single center cell
#'                     CenterCellType - type or subtype of center cell
#'                     C.PositiveMarkers - comma delimited character string of all
#'                                         markers positive in center cell
#'                     NeighborhoodCellType - type or subtype of cells within
#'                                            the neighborhood of the center cell 
#'                                            (default: 30um)
#'                     N.Count - number of cells of NeighborhoodCellType existing
#'                               in neighborhood of center cell (default: 30um)
#'                     Band - interface band in which the center cell lies
#' @param allCenterPopulations  vector of all center populations to be filtered out
#'                              and if not already existing in nbhdCounts, appended to
#'                              nbhdCounts
#' @return  expanded neighborhood counts table including subsets corresponding to
#'          exact center cell types
getNeighborhoodCounts <- function(nbhdCounts, allCenterPopulations, calcUnit){
    allCounts <- tibble()
    for(cct in allCenterPopulations){
        cctCounts <- getCenterCellTypeSubset(nbhdCounts, cct)
        allCounts <- allCounts %>% bind_rows(cctCounts)
    }
    if(nrow(allCounts) == 0){
        return(allCounts)
    }

    ## add missing neighborhood types
    NCTs <- unique(allCounts$NeighborhoodCellType)
    FOV_IDs <- unique(allCounts$FOV_ID)

    allCounts <- allCounts %>%
                 spread(NeighborhoodCellType, N.Count, fill=0) %>%
                 gather(NCTs, key="NeighborhoodCellType", value="N.Count")

    ## add missing FOVs
    tmp  <- allCounts %>%
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

    if(nrow(tmp) > 0){
        allCounts <- allCounts %>% bind_rows(tmp)
    }

    allCounts
}


#' Get mean fraction of a cell subtype in the neighborhood of a center cell type
#' Given a center cell type (e.g., macrophage), a neighborhood cell type (e.g., Tconv8) and
#' a neighborhood cell subpopulation (e.g., Tconv8,PD1,LAG3,TIM3), calculate for each cell the fraction
#' of neighborhood subpopulation out of all cells of neighborhood type. i.e., for this example,
#' calculate in each cell the fraction of Tconv8 cells in the cell's neighborhood (typically 
#' within 30 microns) that are exhausted. Calculate the mean fraction across all macrophage cells
#' in each calculation unit (typically Sample_ID or FOV_ID) 
#' 
#' @param neighborhoodCounts  table of center cell-neighborhood cell type pairs, including at minimum 
#'                            columns: 
#'                              - C.UUID                UUID of a single center cell
#'                              - CenterCellType        cell type or subtype of the cell identified by 
#'                                                      C.UUID
#'                              - C.PositiveMarkers     all markers positive in center cell
#'                              - NeighborhoodCellType  a single cell type/state within the neighborhood 
#'                                                      (30 microns) of center cell
#'                              - N.Count               the number of cells of type NeighborhoodCellType
#'                                                      in neighborhood of center cell
#' @param fracsToCalculate    table describing fractions to be calculated, including columns:
#'                              - `Center Population A`       center cell type of numerator 
#'                              - `Neighborhood Popluation A` neighborhood cell type of numerator
#'                              - `Center Population B`       center cell type of denominator 
#'                              - `Neighborhood Population B` neighborhood cell type of denominator
#' @param calcUnit            data column name containing values for which a single fraction summary
#'                            should be calculated (default: Sample)
#' @param outFile             file that either contains fractions from which data should be
#'                            loaded or to which data should be saved
#' @return table of fraction summaries for each calcUnit
getNeighborhoodFractions <- function(neighborhoodCounts, fracsToCalculate,
                                     calcUnit = "Sample", outFile=NULL){
    if(fileDone(outFile)){
        log_info(paste0("Loading pre-computed neighborhood fractions from file ",outFile))
        return( readRDS(outFile) )
    }

    log_info("No pre-computed neighborhood fractions file found. Generating now...")

    ftc <- fracsToCalculate
    allFracs <- tibble()

    ## check for center cell types that need to be subset (indicated by presence of a comma)
    ## anything after the comma SHOULD be a marker combo, otherwise it would already
    ## exist in the list of 'center' cell types
    allCenterPops <- unique(c(ftc$`Center Population A`, ftc$`Center Population B`))
    for(cct in allCenterPops){
        log_debug(paste0("Adding rows for center cell type ",cct))
        subDat <- getCenterCellTypeSubset(neighborhoodCounts, cct)
        subftc <- ftc %>%
                  filter(`Center Population A` == cct | `Center Population B` == cct)

        ## pull out counts for numerator condition
        num <- subftc %>%
               select(`Cell State ID`, dplyr::matches("Population A")) %>%
               left_join(subDat %>%
                         select(`Center Population A` = CenterCellType,
                                `Neighborhood Population A` = NeighborhoodCellType,
                                everything()),
                         by=c("Center Population A",
                              "Neighborhood Population A")) %>%
               dplyr::rename(NCountA=N.Count) %>%
               filter(!is.na(`Cell State ID`)) %>%
               select(`Cell State ID`,
                      `Center Population A`,
                      `Neighborhood Population A`,
                      FOV_ID,
                      C.UUID,
                      NCountA)

        ## pull out counts for denominator condition
        dnm <- subftc %>%
               select(`Cell State ID`, dplyr::matches("Population B")) %>%
               left_join(subDat %>%
                         select(`Center Population B` = CenterCellType,
                                `Neighborhood Population B` = NeighborhoodCellType,
                                everything()),
                         by=c("Center Population B",
                              "Neighborhood Population B")) %>%
               dplyr::rename(NCountB=N.Count) %>%
               filter(!is.na(`Cell State ID`))

        ## join numerator and denominator
        fracs <- num %>%
                 full_join(dnm, by = intersect(names(num), names(dnm))) %>%
                 filter(NCountB > 0) %>%
                 filter(`Center Population A` != `Center Population B` |
                        `Neighborhood Population A` != `Neighborhood Population B`) %>%
                 mutate(Fraction = NCountA/NCountB) %>%
                 group_by_at(c("Cell State ID",
                               "Center Population A", "Neighborhood Population A",
                               "Center Population B", "Neighborhood Population B",
                             calcUnit)) %>%
                 summarize(MeanFraction = mean(Fraction, na.rm = T))

        allFracs <- allFracs %>% bind_rows(fracs)
    }
    if(length(outFile) > 0 && !is.null(outFile)){
        saveRDS(allFracs, outFile)
    }
    
    allFracs

}

#' Get average number of a certain cell type within the neighborhood of another cell type
#' 
#' Given a table of neighborhood counts, where a row contains a single 'center cell' of 
#' type X and the number of cells of type Y that are within its neighborhood (typically
#' within 30 microns)
#' 
#' @param neighborhoodCounts  table of center/neighborhood cell pairs including at least columns 
#'                            'C.UUID', 'CenterCellType', 'C.PositiveMarkers', 'NeighborhoodCellType', 
#'                            'N.Count'          
#' @param avgsToCalculate     table describing neighborhoods to be averaged
#' @param calcUnit            data column name containing values for which a single fraction summary
#'                                should be calculated (default: Sample)
#' @param outFile             file that either contains averages from which data should be
#'                            loaded or to which data should be saved
#' @return table of neighborhood averages for each calcUnit
getNeighborhoodAverageCounts <- function(neighborhoodCounts, avgsToCalculate,
                                         calcUnit="FOV_ID", outFile=NULL){

    if(fileDone(outFile)){
        log_info(paste0("Loading pre-computed neighborhood averages from file ",outFile))
        return( readRDS(outFile) )
    }

    log_info("No pre-computed neighborhood averages file found. Generating now...")

    calcAvgs <- function(dat, calcUnit){
        dat %>%
        group_by_at(c("CenterCellType", "NeighborhoodCellType", calcUnit)) %>%
        summarize(AvgNbhdCellTypeCount = mean(N.Count, na.rm = T))
    }

    avgs <- tibble()

    subs <- unique(avgsToCalculate$Center)
    for(cct in subs){
        log_debug(paste0("Getting average counts by neighborhood population for center cell type ",cct))
        nbhds <- avgsToCalculate %>%
                 filter(Center == cct) %>%
                 pull(Neighborhood) %>%
                 unique()
        subAvgs <- getCenterCellTypeSubset(neighborhoodCounts, cct) %>%
                   filter(NeighborhoodCellType %in% nbhds) %>%
                   calcAvgs(calcUnit)
        avgs <- avgs %>% bind_rows(subAvgs)
    }

    if(length(outFile) > 0 && !is.null(outFile)){
        saveRDS(avgs, outFile)
    }
    return( avgs )

}


#' Generate and/or load all population metrics of interest for all or a subset of all 
#' calculation units
#'
#' @param dat           tibble of cell level annotated data
#' @param allAnalyses   list of all analyses to be run/loaded where each list contains a tibble
#'                      describing analyses for one calculation type (fraction, density, etc.)
#' @param markers       vector of all markers involved in study
#' @param areas         optional; tibble containing one area value per calculation unit to be
#'                      used in density calculations if included in analyses (default = NULL)
#' @param calcUnit      data column name containing values for which a single value
#'                      should be calculated (default: FOV_ID)
#' @param include       optional; vector of calculation units to include in results; when NULL, 
#'                      all calc units will be included
#' @param nbhdCounts    optional; required when neighborhood analyses included in analysisList; 
#'                      tibble of macrophage neighborhood counts (or other neighborhood counts
#'                      if other center cells included in analysis list); default = NULL
#' @param outdir        optional; if provided and dataStored is TRUE, will load pre-computed 
#'                      calculations from files in this directory if available and if not, will
#'                      write them to this directory 
#' @param numThreads    number of threads
#' @param dataStored    when true, any data already saved in outDir will be loaded and calculations
#'                      will NOT be rerun; also, any new calculations for calc types not already
#'                      run will be saved; WARNING: if you are not sure what data has been previously 
#'                      calculated, it's best to set this to FALSE 
#' @return list of metrics for specified set of calculation units where each element of list
#'         contains data for a single calculation type 
getPopulationMetrics <- function(dat, allAnalyses, markers, areas = NULL, calcUnit = "FOV_ID",
                               include = NULL, nbhdCounts = NULL, outDir = NULL, numThreads = 1,
                               dataStored = TRUE, cellDiveID = "All"){

    countFiles <- fracFiles <- densFiles <- navgFiles <- nfracFiles <- NULL

    if(!is.null(outDir) && dataStored){
        mkdir(outDir)
        if(tolower(cellDiveID) == "all"){
            cdids <- dat %>% pull("CellDive_ID") %>% unique
        } else {
            cdids <- cellDiveID
        }
        countFiles <- file.path(outDir, paste0(cdids, "_counts_per_", calcUnit, ".rda"))
        fracFiles  <- file.path(outDir, paste0(cdids, "_fractions_per_", calcUnit, ".rda"))
        densFiles  <- file.path(outDir, paste0(cdids, "_densities_per_", calcUnit, ".rda"))
        navgFiles  <- file.path(outDir, paste0(cdids, "_neighborhood_averages_per_", calcUnit, ".rda"))
        nfracFiles <- file.path(outDir, paste0(cdids, "_neighborhood_mean_fractions_per_", calcUnit, ".rda"))
    }


    ## get full fractions/densities population list
    allPops <- unique(c(allAnalyses$fractions$Condition,
                        allAnalyses$fractions$Population,
                        allAnalyses$densities$Condition))
    popDat <- list()

    if(any(c("fractions","densities") %in% names(allAnalyses))){
        popDat$counts <- tibble()
        if(is.null(countFiles) || length(countFiles) == 0){
            popDat$counts <- getPopulationCounts(dat, allPops, markers, 
                                                 calcUnit = calcUnit,
                                                 numThreads = numThreads)
        } else {
            for(countFile in countFiles){
                popDat$counts <- popDat$counts %>%
                                 bind_rows(getPopulationCounts(dat, allPops, markers, 
                                                               calcUnit = calcUnit,
                                                               numThreads = numThreads, 
                                                               outFile = countFile)) 
            }
        }
        popDat$counts <- popDat$counts %>%
                         filterForCalculationUnits(calcUnit = calcUnit, include = include)
    }

    ### fractions
    if("fractions" %in% names(allAnalyses)){
        popDat$fractions <- tibble()
        if(is.null(fracFiles) || length(fracFiles) == 0){
            popDat$fractions <- getPopulationFractions(popDat$counts,
                                                       allAnalyses$fractions %>%
                                                         select(Condition, Population),
                                                       calcUnit = calcUnit)
        } else {
            for(fracFile in fracFiles){
                popDat$fractions <- popDat$fractions %>%
                                    bind_rows(getPopulationFractions(popDat$counts,
                                                                     allAnalyses$fractions %>% 
                                                                       select(Condition, Population),
                                                                     calcUnit = calcUnit,
                                                                     outFile = fracFile))
            }
        }
        popDat$fractions <- popDat$fractions %>%
                            filterForCalculationUnits(calcUnit = calcUnit, include = include) %>%
                            full_join(allAnalyses$fractions %>% 
                                        select(`Cell State ID`, Condition, Population),
                                      by = c("Condition", "Population"))
    }

    ### densities
    if("densities" %in% names(allAnalyses)){
        if(is.null(areas) || nrow(areas) == 0 || !calcUnit %in% names(areas)){
            log_warn(paste0("Invalid or empty table of ", calcUnit,
                            " areas. If provided and not empty, make sure column ",
                            calcUnit, " exists in table."))
            popDat$densities <- NULL
        } else {
            if(is.null(densFiles) || length(densFiles) == 0){
                popDat$densities <- getPopulationDensities(popDat$counts,
                                                           areas,
                                                           calcUnit = calcUnit)
            } else {
                popDat$densities <- tibble()
                for(densFile in densFiles){
                    popDat$densities <- popDat$densities %>%
                                        bind_rows(getPopulationDensities(popDat$counts,
                                                                         areas,
                                                                         calcUnit = calcUnit,
                                                                         outFile = densFile))
                }
            }
            popDat$densities <- popDat$densities %>%
                                filterForCalculationUnits(calcUnit = calcUnit, include = include) %>%
                                full_join(allAnalyses$densities %>% 
                                            select(`Cell State ID`, Population = Condition),
                                          by = c("Population"))
        }
    }

    ### make sure ALL conditions exist in neighborhood data
    if(any(c("navgcounts", "nfracs") %in% names(allAnalyses))){
        allCenters <- c(allAnalyses$nfracs$`Center Population A`,
                        allAnalyses$nfracs$`Center Population B`,
                        allAnalyses$navgcounts$`Center`) %>%
                      unique()

        allCounts <- nbhdCounts %>%
                     filter(C.UUID %in% dat$UUID) %>%
                     filterForCalculationUnits(calcUnit = calcUnit, include = include)

        popDat$ncounts <- nbhdCounts %>% #getNeighborhoodCounts(allCounts, allCenters, calcUnit) %>%
                          filter(!is.na(!!as.name(calcUnit)))
    }

    ### neighborhood average counts
    if("navgcounts" %in% names(allAnalyses)){
        if(is.null(navgFiles) || length(navgFiles) == 0){
            popDat$navgcounts <- getNeighborhoodAverageCounts(popDat$ncounts,
                                                              allAnalyses$navgcounts,
                                                              calcUnit = calcUnit)
        } else {
            popDat$navgcounts <- tibble()
            for(navgFile in navgFiles){
                popDat$navgcounts <- popDat$navgcounts %>%
                                     bind_rows(getNeighborhoodAverageCounts(popDat$ncounts,
                                                                            allAnalyses$navgcounts,
                                                                            calcUnit = calcUnit,
                                                                           outFile = navgFile))
            }
            popDat$navgcounts <- popDat$navgcounts %>%
                                 filterForCalculationUnits(calcUnit = calcUnit, include = include) %>%
                                 full_join(allAnalyses$navgcounts %>%
                                             select(CenterCellType = Center,
                                                    NeighborhoodCellType = Neighborhood,
                                                    everything()),
                                           by = c("CenterCellType", "NeighborhoodCellType")) %>%
                                               filter(!is.na(Cell_State))
        }
    }

    ### neighborhood fractions
    if("nfracs" %in% names(allAnalyses)){
        if(is.null(nfracFiles) || length(nfracFiles) == 0){
            popDat$nfracs <- getNeighborhoodFractions(nbhdCounts,
                                                      allAnalyses$nfracs,
                                                      calcUnit = calcUnit)
        } else {
            popDat$nfracs <- tibble()
            for(nfracFile in nfracFiles){
                popDat$nfracs <- popDat$nfracs %>%
                                 bind_rows(getNeighborhoodFractions(nbhdCounts,
                                                                    allAnalyses$nfracs,
                                                                    calcUnit = calcUnit,
                                                                    outFile = nfracFile))
            }
            popDat$nfracs <- popDat$nfracs %>%
                             filterForCalculationUnits(calcUnit = calcUnit, include = include) %>%
                             full_join(allAnalyses$nfracs,
                                         by = intersect(names(popDat$nfracs),
                                                        names(allAnalyses$nfracs))) %>%
                             filter(!is.na(Cell_State))
        }
    }
    popDat
}


getTMEfractions <- function(tme, cellStateID, center, subpop, population, markers){

    fracs <- tibble()

    for(p in c("A", "B")){
        col <- paste0("Neighborhood Population ", p)
        pop <- ifelse(p == "A", subpop, population)
        popDat <- tme %>%
                  filterDataForCondition(pop, markers) %>%
                  mutate(`Cell State ID` = cellStateID,
                         `Center Population A` = center,
                         `Center Population B` = center,
                         !!as.name(col) := pop) %>%
                  group_by_at(names(.)[!grepl("Marker|Count|Classifiers", names(.))]) %>%
                  summarize(!!as.name(paste0("N.pop.",p,"_Count")) := sum(Count))
        if(p == "A"){
            fracs <- popDat
        } else {
            fracs <- fracs %>%
                     full_join(popDat, by = intersect(names(.), names(popDat)))

            naPop <- which(is.na(fracs$`Neighborhood Population A`))
            fracs$`Neighborhood Population A`[naPop]        <- subpop
            fracs$N.pop.A_Count[is.na(fracs$N.pop.A_Count)] <- 0
        }
    }
    calcUnit <- "C.UUID"
    fracs <- fracs %>%
             mutate(Fraction = N.pop.A_Count/N.pop.B_Count) %>%
             filter(!!as.name(groupVar) %in% c("Pos.Env","Neg.Env"))
    fracs

}

