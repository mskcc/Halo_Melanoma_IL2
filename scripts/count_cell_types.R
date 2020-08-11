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

source(file.path(sdir, "R/counts.R"))
source(file.path(sdir, "R/cell_annotation.R"))
source(file.path(sdir, "R/script_utils.R"))
source(file.path(sdir, "R/data_utils.R"))
source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/filter_data.R"))
source(file.path(sdir, "R/format_data.R"))
source(file.path(sdir, "R/load_data.R"))

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript annotate_cells.R 

          [REQUIRED (may be defined on command line OR in manifest file)]
            --cell_data_dir           path to annotated cells RDA files that each 
                                      contain single table where rows are cells
                                      and columns are all data for that single 
                                      cell (one file per sample)
            --cell_type_counts_file   path to output file - XLSX file to contain counts summaries

          [OPTIONAL]
            --manifest            YAML file containing one or more parameter; NOTE: arguments on command
                                  line override manifest arguments!!!
            --meta_dir            path to meta files in XLSX format, required IF meta_data_file is NULL
            --meta_files          comma-delimited list of meta files
            --number_threads      number of threads to use for parallel processes
        \n"
    )
}

## set up required args & defaults 
minReq <- list(c("meta_dir","meta_files"),
               "cell_data_dir",
               "cell_type_counts_file",
               "number_threads")

defaults <- list()

if(!interactive()){
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    ## for testing
    args <- list(manifest  = "input/config/study_config.yaml",
                 cell_data_dir = "preprocessing/03_annotated",
                 cell_type_counts_file = "results/counts/cell_type_counts.xlsx",
                 number_threads = 6,
                 meta_dir = "input/meta")
    args <- processCMD(args, defaults, minReq, usage) 
}

checkRequiredInput(args, minReq)

loadGlobalStudyData(args, conditions = FALSE, analyses = FALSE)

countSumms <- getCountSummaries(annCells, cellTypes, xlsxFile=args$cell_type_counts_file) 
openxlsx::write.xlsx(countSumms, args$cell_type_counts_file)

