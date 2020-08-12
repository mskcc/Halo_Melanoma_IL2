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
            --statistics_conditions_file  YAML file containing all 'conditions' or 
                                          cell states to analyze 
            --cell_data_dir               path to RDA files each containing a single 
                                          table where rows are cells and columns are 
                                          all data for that single cell
            --metrics_dir                 root directory for all area, fractions and 
                                          densities; subdirs will be created for each 
                                          cell region
            --meta_dir                    path to meta files in XLSX format
            --focus                       metrics level on which to focus ['fov'|'interface']
            --fov_area_dir                directory containing fov area files 
            --band_dir                    directory containing infiltration band files


          [OPTIONAL]
            --manifest                YAML file containing one or more parameter; NOTE: 
                                      arguments on command line override manifest arguments!!!         
            --number_threads          number of threads to use for parallel processes
            --cell_dive_id            a single cell dive id for which metrics should be calculated
        \n"
    )
}

## names of required args
minReq <- list("statistics_conditions_file",
               "annotation_config_file",
               c("meta_dir","meta_files"),
               "cell_data_dir",
               "fov_area_dir",
               "band_dir",
               "metrics_dir",
               "tme_by_cell_dir",
               "neighborhood_dir",
               "number_threads",
               "focus")

used <- c("metrics_dir","meta_dir","meta_files","cell_data_dir","fov_area_dir",
          "band_dir", "meta_data_file","number_threads","statistics_conditions_file",
          "tme_by_cell_dir", "neighborhood_dir",
          "cell_dive_id", "focus")

defaults <- list(focus = "fov")

if(!interactive()){
    suppressMessages(library(R.utils))
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- processCMD(list(manifest = "input/config/study_config.yaml",
                            cell_dive_id = "mel_1", 
                            focus = "fov"),
                       defaults, minReq, usage)
}

###############################################
###              CONFIGURATION              ###    
###############################################
aCfg  <- read_yaml(args$annotation_config_file)
cfg   <- resolveConfig(aCfg, args)
logParams(cfg, used)

cellRegions <- c("fov", "interface", "interface inside", "interface outside", "neighborhood")
cellDiveID  <- ifelse(is.null(cfg$cell_dive_id), "All", cfg$cell_dive_id)

###############################################
###           INITIALIZE ALL DATA           ###
###############################################
if(cfg$focus == "fov"){
    all <- list(cuFOV = list(CU = "FOV_ID", CR = c("fov")))

    loadGlobalStudyData(cfg, 
                        analyses = TRUE,
                        conditions = TRUE,
                        questions = TRUE,
                        neighborhoodCounts = TRUE,
                        cellsInTumorNeighborhood = FALSE,
                        tmeSampleStatus = FALSE,
                        tmeCellStatus = FALSE)
    tumorNbhdCells <<- NULL

} else if(cfg$focus == "interface"){

    all <- list(cuFOV = list(CU = c("FOV_ID"), CR = cellRegions[-1]),
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


