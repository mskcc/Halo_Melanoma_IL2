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
            --data_dir              directory containing all halo data RDA files for which area should
                                    be calculated; alternatively, provide a vector of specific files to
                                    --data_files argument
            --fov_area_dir          output directory for FOV area files (one per sample)
            --halo_boundaries_file  path to subdirectories of Halo XML files of boundary annotations
            --meta_dir              path to meta files in XLSX format, required IF meta_data_file is NULL           

          [OPTIONAL]
            --cell_dive_id                  CellDive_ID of sample to be processed; to be used only if a 
                                            single data file is provided
            --data_files                    halo data RDA file; must be specified if data_dir is not
            --manifest                      YAML file containing one or more parameter;
                                            NOTE: arguments on comman line override manifest arguments!!!
            --max_g                         default = 5 [ STILL NOT SURE WHAT THIS IS ]
            --meta_data_file                RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                                            meta_files are NULL
            --number_threads                number of threads to use for parallel processes
        \n"
    )
}

minReq <- list("fov_area_dir",
               c("data_files", "data_dir"),
               c("meta_dir","meta_files","meta_data_file"),
               "halo_boundaries_file")

defaults <- list(max_g = 5, number_threads = 6)

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- defaults
    args$data_files <- "preprocessing/02_rethresholded/mel_1_MegaTableV5b_Excl_SOX10reThres_Rule2__TnullFix.rda"
    args$cell_dive_id         <- "mel_1"
    args$halo_boundaries_file <- "halo_boundaries.rda"
    args$fov_area_dir         <- "processed/metrics/fovs/areas"
    args$meta_dir             <- "input/meta"
}

logParams(args, names(args)[!names(args) %in% c("no-restore", "slave", "args", "file")])

### load all previously created data
annot <- loadStudyAnnotations(metaFiles = getFiles(path = args$meta_dir, pattern = ".xlsx"), 
                              metaDataFile = args$meta_data_file)$flat

boundaries <- readRDS(args$halo_boundaries_file)
dataFiles  <- getFiles(path = args$data_dir, files = args$data_files) 
threads    <- min(args$number_threads, length(dataFiles))
maxG       <- args$max_g 
outD       <- args$fov_area_dir

###
### Calculate and save FOV areas
###
cl <- makeCluster(threads, type="FORK", outfile="")
clusterExport(cl, c("dataFiles","boundaries", "maxG", "outD"), envir=environment())

x <- parLapply(cl, unique(annot$CellDive_ID), function(cdid){

        haloFile <- dataFiles[grepl(paste0("^", cdid, "_"), basename(dataFiles))]
        if(length(haloFile) == 0){ return(NULL) }

        log_info(paste0("Reading file...",haloFile))
        fAreas <- getFOVAreas(cdid, haloFile, boundaries, maxG = maxG)

        outF <- file.path(outD, paste0(cdid, "_fov_areas.rda"))
        log_info(paste0("Saving FOV areas to...", outF))
        saveRDS(fAreas, outF) 
     })

stopCluster(cl)    

