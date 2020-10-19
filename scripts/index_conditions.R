suppressMessages(library(logger))
suppressMessages(library(R.utils))
suppressMessages(library(funr))
suppressMessages(library(yaml))
suppressMessages(library(tidyverse))
suppressMessages(library(xlsx))

log_threshold(DEBUG)

#####################################
#####        SOURCE CODE        #####
#####################################

sdir <- dirname(get_script_path())
log_debug(paste0("loading source files from: ",sdir))

source(file.path(sdir, "R/script_utils.R"))
source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/load_data.R"))
source(file.path(sdir, "R/format_data.R"))
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
            --statistics_conditions_index output XLSX file where indexed conditions
                                          will be saved

          [OPTIONAL]
            --manifest                YAML file containing one or more parameter; NOTE: 
                                      arguments on command line override manifest arguments!!!         
        \n"
    )
}

## names of required args
minReq <- list("statistics_conditions_file", "annotation_config_file", "statistics_conditions_index")
defaults <- list()

if(!interactive()){
    suppressMessages(library(R.utils))
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- processCMD(list(manifest = "input/config/study_config.yaml"),
                       defaults, minReq, usage)
}

###############################################
###              CONFIGURATION              ###    
###############################################
aCfg  <- read_yaml(args$annotation_config_file)
cfg   <- resolveConfig(aCfg, args)
logParams(cfg, names(cfg))

conds <- getConditionsIndex(cfg$statistics_conditions_index, 
                            cfg$statistics_conditions_file, 
                            cfg$arrange_annotation)

