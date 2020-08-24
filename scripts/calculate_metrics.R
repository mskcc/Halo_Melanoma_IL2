suppressMessages(library(logger))
suppressMessages(library(R.utils))
suppressMessages(library(funr))
suppressMessages(library(yaml))
suppressMessages(library(tidyverse))
suppressMessages(library(xlsx))
suppressMessages(library(parallel))

log_threshold(DEBUG)

#####################################
#####        SOURCE CODE        #####
#####################################

sdir <- dirname(get_script_path())
log_debug(paste0("loading source files from: ",sdir))

source(file.path(sdir, "R/constants.R"))
source(file.path(sdir, "R/cell_annotation.R"))
source(file.path(sdir, "R/script_utils.R"))
source(file.path(sdir, "R/data_utils.R"))
source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/format_data.R"))
source(file.path(sdir, "R/metrics.R"))
source(file.path(sdir, "R/load_data.R"))
source(file.path(sdir, "R/filter_data.R"))

#####################################
#####       GET USER INPUT      #####
#####################################

usage <- function(){
    cat("\nUsage:  Rscript calculate_metrics.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)]
            --annotation_config_file      YAML file describing how to arrange and 
                                          index cell annotation
            --band_dir                    directory containing infiltration band files
            --cell_data_dir               path to RDA files each containing a single 
                                          table where rows are cells and columns are 
                                          all data for that single cell
            --fov_area_dir                directory containing fov area files 
            --meta_dir                    path to meta files in XLSX format
            --metrics_dir                 root directory for all area, fractions and 
                                          densities; subdirs will be created for each 
                                          cell region
            --neighborhood_counts_dir     directory of RDA files containing formatted macrophage 
                                          neighborhood counts 
            --neighborhood_dir            directory containing RDA files of all pairwise distances 
                                          between cells, at least those <= 30 microns. generally 
                                          each file will contain a table of all 'center' cells of
                                          a certain type, but this is not a requirement as all files
                                          will be loaded together
            --statistics_conditions_file  YAML file containing all 'conditions' or 
                                          cell states to analyze 
            --statistics_conditions_index XLSX file of pre-indexed conditions
            --statistics_questions_file   XLSX file describing questions to be answered (see docs)
            --tme_by_cell_dir             directory containing RDA files of tumor microenvironment 
                                          assignments; i.e., each cell (see docs for details); NOTE:
                                          not required if cell_dive_id specified and sample is NOT a tumor

          [OPTIONAL]
            --cell_dive_id            a single cell dive id for which metrics should be calculated
            --focus                   metrics level on which to focus ['fov'|'interface']
            --manifest                YAML file containing one or more parameter; NOTE: 
                                      arguments on command line override manifest arguments!!!         
            --number_threads          number of threads to use for parallel processes
        \n"
    )
}

## names of required args
minReq <- list("annotation_config_file", c("band_dir","fov_area_dir"), "cell_data_dir",
               c("meta_dir","meta_files"), "metrics_dir", "neighborhood_dir", 
               "neighborhood_counts_dir",
               "statistics_conditions_file", "statistics_conditions_index", 
               "statistics_questions_file",
               "tme_by_cell_dir") 

defaults <- list(cell_dive_id = "All", focus = NA, number_threads = 1)

if(!interactive()){
    suppressMessages(library(R.utils))
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- defaults
    args$annotation_config_file <- "input/config/annotation_config.yaml"
    args$band_dir <- "preprocessing/03_infiltration_band_data"
    args$cell_data_dir <- "preprocessing/03_annotated"
    args$fov_area_dir <- "processed/metrics/fovs/areas"
    args$meta_dir <- "input/meta"
    args$metrics_dir <- "processed/metrics"
    args$neighborhood_counts_dir <- "preprocessing/04_neighborhoods/macro_nbhd_counts"
    args$neighborhood_dir <- "preprocessing/04_neighborhoods/all_vs_all"
    args$statistics_conditions_file <- "input/config/stats_conditions.xlsx"
    args$statistics_conditions_index <- "results/statistics/conditions_index.xlsx"
    args$statistics_questions_file <- "input/config/stats_questions.xlsx"
    args$tme_by_cell_dir <- "preprocessing/05_microenvironments/cell_status"
    args$cell_dive_id <- "mel_1"
    args$focus <- "fov"
}

###############################################
###              CONFIGURATION              ###    
###############################################
aCfg  <- read_yaml(args$annotation_config_file)
cfg   <- resolveConfig(aCfg, args)

logParams(cfg, sort(names(cfg)[!names(cfg) %in% c("no-restore", "slave", "args", "file")]))

cellRegions <- c("fov", "interface", "interface inside", "interface outside", "neighborhood")

###############################################
###           INITIALIZE ALL DATA           ###
###############################################
if(tolower(cfg$focus) == "fov"){
    all <- list(cuFOV = list(CU = "FOV_ID", CR = c("fov")))

    loadGlobalStudyData(cfg, 
                        analyses = TRUE,
                        conditions = TRUE,
                        questions = TRUE,
                        neighborhoodCounts = TRUE,
                        cellsInTumorNeighborhood = FALSE,
                        tmeCellStatus = FALSE)

} else if(tolower(cfg$focus) == "interface"){

    all <- list(cuFOV = list(CU = c("FOV_ID"), CR = cellRegions[-1]),
                cuFOVBand = list(CU = c("FOV_ID","Band"), CR = cellRegions[-1]),
                cuSampleBand = list(CU = c("Sample_ID","Band"), CR = cellRegions[-1]))

    loadGlobalStudyData(cfg, all = T)

} else {

    all <- list(cuFOV = list(CU = c("FOV_ID"), CR = cellRegions),
                cuFOVBand = list(CU = c("FOV_ID","Band"), CR = cellRegions[-1]),
                cuSampleBand = list(CU = c("Sample_ID","Band"), CR = cellRegions[-1]))

    loadGlobalStudyData(cfg, all = T)

} 


###############################################
###          CALCULATE ALL METRICS          ###
###############################################

lapply(all, function(x){
    lapply(x$CR, function(cr, cuCols){
        precalculateMetrics(annCells, 
                            nbhdCounts, 
                            tumorNbhdCells, 
                            analysisList,
                            cfg$metrics_dir, 
                            cr, 
                            cuCols, 
                            numThreads = cfg$number_threads,
                            cellDiveID = cfg$cell_dive_id) 
    }, x$CU)
})


