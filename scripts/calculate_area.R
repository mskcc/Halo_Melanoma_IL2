suppressMessages(library(logger))
suppressMessages(library(R.utils))
suppressMessages(library(funr))
suppressMessages(library(yaml))
suppressMessages(library(tidyverse))
suppressMessages(library(xlsx))
suppressMessages(library(randtoolbox))
suppressMessages(library(rgeos))
suppressMessages(library(sp))
suppressMessages(library(parallel))
suppressMessages(library(SearchTrees))

log_threshold(DEBUG)

#####################################
#####        SOURCE CODE        #####
#####################################

sdir <- dirname(get_script_path())
log_debug(paste0("loading source files from: ",sdir))

source(file.path(sdir, "R/constants.R"))
source(file.path(sdir, "R/script_utils.R"))
source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/metrics.R"))
source(file.path(sdir, "R/load_data.R"))
source(file.path(sdir, "R/validation.R"))
source(file.path(sdir, "R/data_utils.R"))

usage <- function(){

    cat("\nUsage:  Rscript calculate_area.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --halo_boundaries_file  path to subdirectories of Halo XML files of boundary annotations
            --fov_area_dir          output directory for FOV area files (one per sample)
            --band_dir              output directory infiltration band assignments (one per sample) 
            --infiltration_area_dir output directory for infiltration band area files (one per sample)
            --focus                 character indicating which area level to focus on [fov|interface]

          [OPTIONAL]
            --manifest                      YAML file containing one or more parameter; 
                                            NOTE: arguments on comman line override manifest arguments!!!         
            --data_files                    halo data RDA file; must be specified if data_dir is not
            --data_dir                      directory containing all halo data RDA files for which area should 
                                            be calculated
            --meta_dir                      path to meta files in XLSX format, required IF meta_data_file is NULL
            --meta_data_file                RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                                            meta_files are NULL
            --max_g                         default = 5 [ STILL NOT SURE WHAT THIS IS ]
            --band_width                    width of each interface interval for which area should be calculated
            --max_distance_from_interface   max width in microns of interface, inside and 
                                            outside (total width = 2*this value)
            --number_threads                number of threads to use for parallel processes
        \n"
    )
}

minReq <- list(c("meta_dir","meta_files","meta_data_file"),
               "halo_boundaries_file",
               "fov_area_dir",
               "band_dir",
               "infiltration_area_dir",
               c("data_files", "data_dir"),
               "number_threads",
               "focus")

used <- c("data_files", "data_dir", "meta_dir","meta_files","meta_data_file","halo_boundaries_file",
          "fov_area_dir", "band_dir", "infiltration_area_dir", "cell_dive_id","focus",
          "number_threads","max_g","band_width","max_distance_from_interface")

defaults <- list(focus = "fov", max_g = 5, band_width = 10, max_distance_from_interface = 360, number_threads = 6)

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- processCMD(read_yaml("input/config/study_config.yaml"), defaults, minReq, usage)
    args$data_files <- "preprocessing/02_rethresholded/mel_1_MegaTableV5b_Excl_SOX10reThres_Rule2__TnullFix.rda"
    args$cell_dive_id <- "mel_1"
}

logParams(args, used)

### load all previously created data
annot <- loadStudyAnnotations(metaFiles = getFiles(path = args$meta_dir, pattern = ".xlsx"), 
                              metaDataFile = args$meta_data_file)$flat

boundaries <- readRDS(args$halo_boundaries_file)
dataFiles  <- getFiles(path = args$data_dir, files = args$data_files) 
threads    <- min(args$number_threads, length(dataFiles))

###
### FOV area
###
if(args$focus == 'fov'){
    getSampleFOVAreas <- function(cdid, dataFile, boundaries, maxG = 5){
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

    if(length(dataFiles) > 1 && threads > 1){
    
        cl <- makeCluster(threads, type="FORK", outfile="")
        clusterExport(cl, c("dataFiles","boundaries", "maxG", args$fov_area_dir), envir=environment())

        x <- parLapply(cl, unique(annot$CellDive_ID), function(cdid){
                 fPat <- paste0("^", cdid, "_")
                 if(!any(grepl(fPat, basename(dataFiles)))){ return(NULL) }
                 haloFile <- dataFiles[grepl(fPat, basename(dataFiles))]
                 log_debug(paste0("Reading from...",haloFile))
                 fAreas <- getSampleFOVAreas(cdid, haloFile, boundaries, maxG = args$max_g)
                 log_debug(paste0("Saving to...",file.path(args$fov_area_dir, paste0(cdid, "_fov_areas.rda"))))
                 saveRDS(fAreas, file.path(args$fov_area_dir, paste0(cdid, "_fov_areas.rda")))
             })

        stopCluster(cl)    

    } else {

        for(cdid in unique(annot$CellDive_ID)){
            fPat <- paste0("^", cdid, "_")
            if(!any(grepl(fPat, basename(dataFiles)))){ next }
            haloFile <- dataFiles[grepl(fPat, basename(dataFiles))]
            log_debug(haloFile)
            fAreas <- getSampleFOVAreas(cdid, haloFile, boundaries, maxG = args$max_g)
            saveRDS(fAreas, file.path(args$fov_area_dir, paste0(cdid, "_fov_areas.rda")))
        }

    }
} else if(args$focus == 'interface'){

    ###
    ### Infiltration area
    ###
    numBands <- args$max_distance_from_interface/args$band_width
    bins <- (-numBands:numBands)*args$band_width

    for(cdid in unique(annot$CellDive_ID)){
        
        fPat         <- paste0("^", cdid, "_")
        bandAreaFile <- file.path(args$infiltration_area_dir, paste0(cdid, "_infiltration_band_areas.rda")) 
        bandFile     <- file.path(args$band_dir, paste0(cdid, "_infiltration_bands.rda"))

        if(!any(grepl(fPat, basename(dataFiles)))){ next }

        haloFile <- dataFiles[grepl(fPat, basename(dataFiles))]
        log_debug(paste0("Reading from...",haloFile))
        sbands <- addInterfaceBandInfo(loadHaloDataFile(haloFile, filterExclusions = TRUE), 
                                       haloAnnotations = boundaries,
                                       interfaceBins=bins, 
                                       numThreads = args$number_threads) 

        if(!is.null(sbands) && nrow(sbands) > 0){
            sbands <- sbands %>%
                      select(UUID, CellDive_ID, FOV_number, X, Y, Distance, Band) %>%
                      unique()
    
            sArea  <- calculateInterfaceArea(sbands, haloAnnotations = boundaries, 
                                             maxG = args$max_g, interfaceBins = bins,
                                             numThreads = args$number_threads)
            bandDat <- sbands %>%
                       select(CellDive_ID, FOV_number, UUID, X, Y, Distance, Band) %>%
                       left_join(sArea, by=c("CellDive_ID", "FOV_number", "Band")) %>%
                       rename(FOVBandArea = Area)
        } else {
            sdat  <- loadHaloDataFile(haloFile, filterExclusions = TRUE) 
            sArea <- sdat %>%
                     select(CellDive_ID, FOV_number) %>%
                     unique() %>%
                     mutate(Band = NA, Area = NA)
            bandDat <- sdat %>%
                       convertCellMinMaxToMidpoints() %>% 
                       select(CellDive_ID, FOV_number, UUID, X, Y) %>%
                       unique() %>%
                       mutate(Distance = NA, Band = NA, FOVBandArea = NA)
        }
        log_debug(paste0("Saving to...", bandAreaFile))
        saveRDS(sArea, bandAreaFile) 

        log_debug(paste0("Saving to...", bandFile))
        saveRDS(bandDat, bandFile)  

    }
}
