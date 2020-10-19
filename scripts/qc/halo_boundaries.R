suppressMessages(library(cowplot))
suppressMessages(library(funr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(logger))
suppressMessages(library(parallel))
suppressMessages(library(rgeos))
suppressMessages(library(R.utils))
suppressMessages(library(sp))
suppressMessages(library(tidyverse))
suppressMessages(library(xlsx))
suppressMessages(library(XML))
suppressMessages(library(yaml))

log_threshold(DEBUG)

sdir <- dirname(dirname(get_script_path()))
log_debug(paste0("loading source files from: ",sdir))

source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/load_data.R"))
source(file.path(sdir, "R/boundaries.R"))
source(file.path(sdir, "R/script_utils.R"))
source(file.path(sdir, "R/spatial_utils.R"))

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript mark_exclusions.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)]
            --raw_data_dir             directory containing raw halo object analysis data in RDA 
                                       format
            --halo_boundaries_rda_dir  output directory; one file per CellDive_ID will be written, 
                                       each containing a list named by FOV, each element containing
                                       another list named by boundary type code ('tumB', 'epiB', 
                                       'excB', 'glsB') where each of those elements is a tibble 
                                       containing coordinates of a single region
            --meta_dir                 path to meta files in XLSX format, required IF meta_data_file 
                                       is NULL
            --qc_dir                   QC directory, where plots will be saved

          [OPTIONAL]
            --cell_dive_id      a single CellDive_ID for which boundary files should be parsed and
                                formatted
            --manifest          YAML file containing one or more parameter; NOTE: arguments on command 
                                line override manifest arguments!!!         
            --meta_files        comma-delimited list of meta files 
            --meta_data_file    RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                                meta_files are NULL
            --number_threads    number of threads to use for parallel processes
        \n"
    )
}

minReq   <- list("halo_boundaries_rda_dir",
                 c("meta_dir","meta_files","meta_data_file"),
                 "qc_dir",
                 "raw_data_dir")

defaults <- list(number_threads = 2)
args     <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, names(args)[!names(args) %in% c("no-restore", "slave", "args", "file")])

fovFile <- getFiles(path = args$meta_dir, files = args$meta_files, pattern = "FOVs.xlsx")
fovAnn  <- read.xlsx(fovFile, 1, check.names = F) %>%
           as_tibble() 

ids <- args$cell_dive_id
threads <- 1
if(is.null(ids)){
    ids <- unique(fovAnn$CellDive_ID)
    threads <- args$number_threads
}

## get table of expected vs observed numbers of boundaries
discrep <- boundaryDiscrepancies(args$halo_boundaries_rda_dir, fovAnn)

## pseudo-colored plots of all FOVs for given ID(s) including
## manual boundary annotations
qcHaloBoundaries(ids, 
                 args$raw_data_dir, 
                 "SOX10", 
                 args$halo_boundaries_rda_dir, 
                 discrep, 
                 args$qc_dir, 
                 threads = threads) 


