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
            --annotated_cells_file  path to RDA file (output) that will contain a single table where rows are
                                    cells and columns are all data for that single cell

          [OPTIONAL]
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command 
                                  line override manifest arguments!!!         
            --data_dir            path to processed, exclusion-marked RDA files of formatted halo object data 
            --data_files          full paths to each file to be included in analysis
            --meta_dir            path to meta files in XLSX format, required IF meta_data_file is NULL
            --meta_files          comma-delimited list of meta files 
            --meta_data_file      RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                                  meta_files are NULL
            --number_threads      number of threads to use for parallel processes
        \n"
    )
}

minReq <- list(c("meta_dir","meta_files","meta_data_file"),
               c("data_dir","data_files"),
               "annotated_cells_file",
               "number_threads")

used <- c("annotated_cells_file","data_dir","data_files","meta_dir","meta_files",
          "number_threads", "overwrite","meta_data_file")

defaults <- list('overwrite' = FALSE, 'number_threads' = 4)

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- processCMD(list(manifest="input/config/study_config.yaml"), defaults, minReq, usage) 
}

logParams(args, names(args)[!names(args) %in% c("no-restore", "slave", "args", "file")])

###
### get meta data 
###
metaFiles <- getFiles(path = args$meta_dir, files = args$meta_files, pattern = ".xlsx")
dataFiles <- getFiles(path = args$data_dir, files = args$data_files, pattern = ".rda")

annDat <- annotateCells(args$annotated_cells_file, dataFiles = dataFiles, 
                        metaFiles = metaFiles, metaDataFile = args$meta_data_file, 
                        numThreads = args$number_threads, filterExclusions = TRUE, 
                        controlMarker = "DAPI")

log_debug("Annotated:")
log_debug(paste("    UUIDs:\t", nrow(annDat)))
log_debug(paste("    FOV_IDs:\t", length(unique(annDat$FOV_ID))))
log_debug(paste("    Sample_IDs:\t", length(unique(annDat$Sample_ID))))
log_info("Done.")
