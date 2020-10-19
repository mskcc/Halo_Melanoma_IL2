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
            --band_dir              output directory infiltration band assignments (one per sample) 
            --data_dir              directory containing all halo data RDA files for which area should 
                                    be calculated; alternatively, a specific vector of files can be
                                    passed to --data_files
            --halo_boundaries_file  path to subdirectories of Halo XML files of boundary annotations
            --infiltration_area_dir output directory for infiltration band area files (one per sample)
            --meta_dir              path to meta files in XLSX format, required IF meta_data_file is NULL 

          [OPTIONAL]
            --cell_dive_id                  CellDive_ID of sample to be processed; to be used only if a 
                                            single data file is provided
            --band_width                    width of each interface interval for which area should be calculated;
                                            default = 10
            --data_files                    halo data RDA file; must be specified if data_dir is not
            --manifest                      YAML file containing one or more parameter;
                                            NOTE: arguments on comman line override manifest arguments!!!
            --max_distance_from_interface   max width in microns of interface, inside and
                                            outside (total width = 2*this value)
                                            be calculated
            --max_g                         default = 5 [ STILL NOT SURE WHAT THIS IS ]
            --meta_data_file                RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                                            meta_files are NULL
            --number_threads                number of threads to use for parallel processes
        \n"
    )
}

minReq <- list(c("band_dir","fov_area_dir"),
               c("data_files", "data_dir"),
               c("meta_dir","meta_files","meta_data_file"),
               "halo_boundaries_file")

defaults <- list(max_g = 5, 
                 band_width = 10, 
                 max_distance_from_interface = 360, 
                 number_threads = 6)

###########################
###    CONFIGURATION    ###
###########################

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

boundaries <- readRDS(args$halo_boundaries_file)
dataFiles  <- getFiles(path = args$data_dir, files = args$data_files) 
threads    <- min(args$number_threads, length(dataFiles))
numBands   <- args$max_distance_from_interface/args$band_width
bins       <- (-numBands:numBands)*args$band_width

logParams(args, names(args)[!names(args) %in% c("no-restore", "slave", "args", "file")])

### load all previously created data
annot <- loadStudyAnnotations(metaFiles = getFiles(path = args$meta_dir, pattern = ".xlsx"),
                              metaDataFile = args$meta_data_file)$flat

#
## for each cell dive ID, measure distances between each cell and the 
## closest point on a tumor boundary in order to bin cells into bands
## around said boundary; for each band in each FOV, calculate area. 
## for each sample, two files are created: one is all band info associated
## with every cell and one is a pared down table containing only area values
#
for(cdid in unique(annot$CellDive_ID)){
        
    fPat     <- paste0("^", cdid, "_")
    haloFile <- dataFiles[grepl(fPat, basename(dataFiles))]
    if(length(haloFile) == 0){ next }

    log_info("Reading Halo data from... ")
    log_info(paste0("  ",haloFile))
    sdat <- loadHaloDataFile(haloFile, filterExclusions = TRUE)

    ## for each cell, measure distance to closest point on tumor boundary
    ## and based on that distance, assign the cell to a boundary band/interval
    sbands <- addInterfaceBandInfo(sdat,
                                   haloAnnotations = boundaries,
                                   interfaceBins=bins, 
                                   numThreads = args$number_threads) 

    if(is.null(sbands) || nrow(sbands) == 0){
        log_debug("No band data to save. (No tumor boundaries?)")
        next
    }

    sbands <- sbands %>%
              select(UUID, CellDive_ID, FOV_number, X, Y, Distance, Band) %>%
              unique()

    ## calculate area for each band/interval in each FOV    
    sArea  <- calculateInterfaceArea(sbands, 
                                     haloAnnotations = boundaries, 
                                     maxG = args$max_g, 
                                     interfaceBins = bins,
                                     numThreads = args$number_threads)

    ## combine all cell level info
    bandDat <- sbands %>%
               select(CellDive_ID, FOV_number, UUID, X, Y, Distance, Band) %>%
               left_join(sArea, by=c("CellDive_ID", "FOV_number", "Band")) %>%
               rename(FOVBandArea = Area)
   
    bandAreaFile <- file.path(args$infiltration_area_dir, paste0(cdid, "_infiltration_band_areas.rda"))     
    log_debug(paste0("Saving to...", bandAreaFile))
    saveRDS(sArea, bandAreaFile) 

    bandFile     <- file.path(args$band_dir, paste0(cdid, "_infiltration_bands.rda"))
    log_debug(paste0("Saving to...", bandFile))
    saveRDS(bandDat, bandFile)  

}
