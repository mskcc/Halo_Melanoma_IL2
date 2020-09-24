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
source(file.path(sdir, "R/filter_data.R"))
source(file.path(sdir, "R/format_data.R"))
source(file.path(sdir, "R/metrics.R"))
source(file.path(sdir, "R/load_data.R"))
source(file.path(sdir, "R/statistics.R"))

#####################################
#####       GET USER INPUT      #####
#####################################

usage <- function(){

 cat("\nUsage:  Rscript report_statistics.R 
            
    [REQUIRED (may be defined on command line OR in manifest file)] 
      --annotation_config_file      YAML file describing how conditions are to be 
                                    arranged and indexed
      --band_dir                    path to RDA files, each containing a table of band 
                                    assignments and band areas for each cell
      --cell_data_dir               path to RDA files, each containing a single 
                                    table where rows are cells and columns are 
                                    all data for that single cell
      --fov_area_dir                path to RDA files, each containing table of FOVs and 
                                    total FOV area for all FOVs in a single sample 
      --meta_dir                    path to meta files in XLSX format
      --metrics_dir                 root directory for all area, fractions and 
                                    densities; subdirs will be created for each 
                                    cell region
      --neighborhood_dir            directory of cell to cell distances, specifically distances 
                                    between tumor cells and immune cells
      --neighborhood_counts_dir     directory of macro neighborhood counts files
      --statistics_conditions_file  XLSX file listing all cell states/conditions to 
                                    compare between two sample groups
      --statistics_config_file      YAML file containing stats config such as filters and
                                    calculations to run stats on (see docs for details)
      --statistics_conditions_index XLSX file with pre-indexed cell states/conditions
      --statistics_questions_file   XLSX file outlining all questions/comparisons for 
                                    which stats should be run
      --statistics_tables_dir       output directory where XLSX files of results should 
                                    be written

    [OPTIONAL]
      --question          a single QuestionNumber from statistics_questions_file to run stats on
      --manifest          YAML file containing one or more parameter; NOTE: arguments 
                          on command line override manifest arguments!!!        
      --number_threads    number of threads to use for parallel processes
      --tme_by_cell_dir   directory containing files of tumor microenvironment assignments by cell, 
                          required if one or more question to be run involves grouping or filtering
                          on these assignments
      --tme_by_sample_dir directory containing files of tumor microenvironment assignments by sample,
                          required if one or more question to be run involved such groupings or filtering
                          on these assignments 
  \n"
 )

}

## names of required args
minReq <- list("annotation_config_file", 
               "cell_data_dir", "band_dir", "fov_area_dir",
               c("meta_dir","meta_files","meta_data_file"),
               "metrics_dir",
               "neighborhood_counts_dir",
               "neighborhood_dir",
               "statistics_questions_file",
               "statistics_conditions_file", 
               "statistics_conditions_index",
               "statistics_config_file",
               "statistics_tables_dir")

defaults <- list(number_threads = 1)

used <- c(unlist(minReq), "question", "manifest", "number_threads", 
          "tme_by_cell_dir")

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

#####################################
###   CONFIGURATION & DATA INIT   ###
#####################################

cfg <- resolveConfig(args, 
                     read_yaml(args$annotation_config_file), 
                     read_yaml(args$statistics_config_file))
 
logParams(cfg, used) 
mkdir(cfg$statistics_tables_dir)

## set nbhd stuff to FALSE
stDat <- loadStudyData(cfg, 
                       analyses = T, 
                       questions = T, 
                       cellsInTumorNeighborhood = T, 
                       tmeSampleStatus = T) 

stDat$analysisList <- stDat$analysisList[c("fractions", "densities")]

######################################
###        ANSWER  QUESTIONS       ###
######################################
log_info(paste0("Running stats on FOV_ID level..."))

for(q in names(stDat$allQuestions)){
    log_info(paste0("QUESTION: ",q))

    allAnalyses <- stDat$analysisList

    ## for tumor neighborhood analyses, only run fractions analysis 
    if(any(tolower(unlist(stDat$allQuestions[[q]])) == "neighborhood")){ 
        allAnalyses <- allAnalyses[names(allAnalyses) == "fractions"] 
    }

    qRes <- tryCatch({
               compareSampleGroups(stDat$allQuestions[[q]], 
                                   stDat$annCells, 
                                   stDat$sampAnn, 
                                   allAnalyses,
                                   stDat$markers, 
                                   cfg$metrics_dir, 
                                   cfg$results_filters[[cfg$use_filter]], 
                                   filtName = cfg$use_filter,
                                   calcUnit = "FOV_ID",
                                   nbhdCounts = stDat$nbhdCounts, 
                                   tumorNbhdCells = stDat$tumorNbhdCells)
              }, error = function(e){
                  log_error(e)
                  log_error(paste0("Stats failed for question ",q))
             })

    if(is.null(qRes)){ next }

    ### Save results
    outFile <- file.path(cfg$statistics_tables_dir, paste0(q,".xlsx"))
    log_debug("writing XLSX file...")
    writeStatsQuestionXLSX(qRes, outFile)
    log_info("done.")
}

