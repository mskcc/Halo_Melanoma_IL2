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
source(file.path(sdir, "R/cell_annotation.R"))
source(file.path(sdir, "R/format_data.R"))

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript annotate_cells.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --annotated_cells_file  path to RDA file (output) that will contain a 
                                    single table where rows are cells and columns 
                                    are all data for that single cell
            --control_marker        name of marker used as control; all cells negative for 
                                    this marker are removed
            --data_dir              path to processed, exclusion-marked RDA files of 
                                    formatted halo object data; required if annotated_cells_file 
                                    does NOT exist and if data_files is NULL 
            --meta_dir              path to meta files in XLSX format, required IF meta_data_file 
                                    is NULL

          [OPTIONAL]
            --force_reannotation  when TRUE, will start with rethresholded halo object analysis files
                                  even if annotated_cells_file already exists; default: FALSE
            --data_files          full paths to each file to be included in analysis
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command 
                                  line override manifest arguments!!!         
            --meta_files          comma-delimited list of meta files 
            --meta_data_file      RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                                  meta_files are NULL
            --number_threads      number of threads to use for parallel processes
        \n"
    )
}

minReq   <- list("annotated_cells_file",
                 "control_marker",
                 c("data_dir","data_files"),
                 c("meta_dir","meta_files","meta_data_file"),
                 "number_threads")

defaults <- list(number_threads = 4)

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, names(args)[!names(args) %in% c("no-restore", "slave", "args", "file")])

###
### get meta data 
###
metaFiles <- getFiles(path = args$meta_dir, files = args$meta_files, pattern = ".xlsx")
dataFiles <- getFiles(path = args$data_dir, files = args$data_files, pattern = ".rda")

log_info("Annotating cells in file(s): ")
for(df in dataFiles){ log_info(paste0("  ", df)) }

annDat <- annotateCells(args$annotated_cells_file, 
                        dataFiles = dataFiles, 
                        metaFiles = metaFiles, 
                        metaDataFile = args$meta_data_file, 
                        numThreads = args$number_threads, 
                        filterExclusions = TRUE, 
                        controlMarker = args$control_marker,
                        forceReannotation = args$force_reannotation)

log_debug("Annotated:")
log_debug(paste("    UUIDs:\t", nrow(annDat)))
log_debug(paste("    FOV_IDs:\t", length(unique(annDat$FOV_ID))))
log_debug(paste("    Sample_IDs:\t", length(unique(annDat$Sample_ID))))
log_info("Done.")
