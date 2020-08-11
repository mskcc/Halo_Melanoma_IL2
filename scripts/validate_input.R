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

source(file.path(sdir,"R/script_utils.R"))
source(file.path(sdir,"R/validation.R"))
source(file.path(sdir,"R/file_utils.R"))
source(file.path(sdir,"R/data_utils.R"))
source(file.path(sdir,"R/load_data.R"))
source(file.path(sdir,"R/cell_annotation.R"))
source(file.path(sdir,"R/format_data.R"))

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript annotate_cells.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --meta_dir     directory containing all required XLSX meta files
            --statistics_questions_file  XLSX file containing questions to be answered,
                                         including filtering criteria for each sample group
                                         in each question
          [OPTIONAL]
            --manifest     YAML file containing custom config
        \n"
    )
}

minReq   <- c("meta_dir")
used     <- c(minReq) 
defaults <- list()

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage, sourceFile = NULL)

if(args$debug == FALSE){ log_threshold(INFO) }

logParams(args, used)

cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n")
log_info("\nValidating all meta data...\n")
validateInput(args)


