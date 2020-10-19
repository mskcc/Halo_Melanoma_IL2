suppressMessages(library(funr))
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

sdir <- dirname(get_script_path())
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
            --halo_boundaries_xml_dir  directory containing subdirectories of all Halo boundary
                                       annotations in XML format; subdirectories include 'interface',
                                       'exclusions', 'glass', 'epidermis' and file names are in 
                                       format [CellDive_ID]_Spot[FOV_number].annotations
            --halo_boundaries_rda_dir  output directory; one file per CellDive_ID will be written, 
                                       each containing a list named by FOV, each element containing
                                       another list named by boundary type code ('tumB', 'epiB', 
                                       'excB', 'glsB') where each of those elements is a tibble 
                                       containing coordinates of a single region
            --meta_dir                 path to meta files in XLSX format, required IF meta_data_file 
                                       is NULL

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

minReq   <- list("halo_boundaries_xml_dir",
                 "halo_boundaries_rda_dir",
                 c("meta_dir","meta_files","meta_data_file"),
                 "number_threads")
defaults <- list(number_threads = 2)
args     <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, names(args)[!names(args) %in% c("no-restore", "slave", "args", "file")])

## NOTE: args$boundary_colors contains color map for Cohort 2. Color map for Cohort 1
##       is hard-coded in getAllHaloBoundaries()
if(is.null(args$boundary_colors)){
    args$boundary_colors <- list('65535' = 'Exc', '65280' = 'Epi', '255' = 'Gls', '16776960' = 'Tum')
}

ids <- args$cell_dive_id
threads <- 1
if(is.null(ids)){
    fovFile <- getFiles(path = args$meta_dir, files = args$meta_files, pattern = "FOVs.xlsx")
    ids <- read.xlsx(fovFile, 1, check.names = F) %>%
           as_tibble() %>%
           pull(CellDive_ID) %>% unique()
    threads <- args$number_threads
}

## map cohort 1 cell dive IDs to the corresponding IDs in annotation file names
idMap <- c(melanoma_untreated = "Untreated", melanoma_PR = "PR", melanoma_CR = "CR")

cl <- makeCluster(threads, type = "FORK", outfile = "")
clusterExport(cl, c('ids', 'args'), envir = environment())
tmp <- parLapply(cl, ids, function(id){
           print(id)
           bounds <- list()
           bounds <- getAllHaloBoundaries(ifelse(id %in% names(idMap), idMap[[id]], id), 
                                          getFiles(args$halo_boundaries_xml_dir), 
                                          boundaryColors = args$boundary_colors)
           outFile <- paste0(id, "_halo_boundaries.rda")
           saveRDS(bounds, file.path(args$halo_boundaries_rda_dir, outFile))
       })
stopCluster(cl)

