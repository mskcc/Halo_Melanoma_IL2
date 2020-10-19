suppressMessages(library(logger))
suppressMessages(library(R.utils))
suppressMessages(library(funr))
suppressMessages(library(yaml))
suppressMessages(library(tidyverse))
suppressMessages(library(xlsx))
suppressMessages(library(parallel))

log_threshold(DEBUG)

sdir <- dirname(get_script_path())
log_debug(paste0("loading source files from: ",sdir))

source(file.path(sdir, "R/script_utils.R"))
source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/load_data.R"))
source(file.path(sdir, "R/exclusions.R"))

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript mark_exclusions.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --control_marker           name of marker used as control; all cells negative for 
                                       this marker are removed
            --drift_summary_dir        path to directory containing drift summary files
            --drift_threshold          maximum allowed percentage drift/loss to include in analyses
            --exclusion_data_dir       output directory; path to processed, exclusion-marked  
                                       RDA file of formatted halo object data
            --halo_boundaries_rda_dir  path to RDA files of halo boundaries (exclusions and interface) 
            --meta_dir                 path to meta files in XLSX format, required IF 
                                       meta_data_file is NULL
            --pad                      numeric; width in microns of pad around each FOV; cells falling 
                                       outside this border will be excluded
            --raw_data_dir             directory containing RDA files of formatted marker-level 
                                       Halo data; required when raw_data_files is NULL

          [OPTIONAL]
            --cell_dive_id      a single CellDive_ID to be marked for exclusion in the associated
                                data file in --raw_data_dir; when NULL (default) all files in that
                                directory will be processed
            --manifest          YAML file containing one or more parameter; NOTE: arguments on command 
                                line override manifest arguments!!!         
            --meta_files        comma-delimited list of meta files 
            --meta_data_file    RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                                meta_files are NULL
            --number_threads    number of threads to use for parallel processes
        \n"
    )
}

minReq   <- list("control_marker",
                 "data_dir",
                 c("meta_dir","meta_files","meta_data_file"),
                 "number_threads",
                 c("raw_data_files", "raw_data_dir"))

defaults <- list(number_threads = 4)

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, names(args)[!names(args) %in% c("no-restore", "slave", "args", "file")])

###
### get halo data & meta data 
###
metaFiles <- getFiles(path = args$meta_dir, files = args$meta_files, pattern = ".xlsx")
rawDataFiles <- getFiles(path = args$raw_data_dir, files = args$raw_data_files, pattern = ".rda")
if(!is.null(args$cell_dive_id)){
    rawDataFiles <- rawDataFiles[grep(paste0("^", args$cell_dive_id, "_"), basename(rawDataFiles))]
}

sampAnn <- read.xlsx(metaFiles[grepl("FOVs", metaFiles)], 1, check.names = F) %>% 
           as_tibble() %>%
           select(CellDive_ID, FOV_number, FOV_exclusion, Marker_exclusion) %>%
           filter(!is.na(FOV_exclusion) | !is.na(Marker_exclusion))

log_info("Marking for exclusion cells in file(s): ")
for(rdf in rawDataFiles){ log_info(paste0("  ", rdf)) }
    dat <- loadHaloDataFile(rdf, controlMarker = args$control_marker)
    

