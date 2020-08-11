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
source(file.path(sdir, "R/validation.R"))

usage <- function(){

    cat("\nUsage:  Rscript calculate_area.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --meta_dir            path to meta files in XLSX format, required if meta_data_file is NULL
            --meta_data_file      RDA or XLSX file to which newly flattened data should be saved

          [OPTIONAL]
            --meta_files          comma-delimited character string with each element containing path
                                  to one meta file; use this when there are multiples of one or more meta
                                  file in a directory
            --manifest            YAML file containing one or more required parameter; NOTE: arguments  
                                  on command line override manifest arguments!!!         
        \n"
    )
}

dflt <- list(meta_dir = NULL,
             meta_files = NULL,
             manifest = NULL,
             meta_data_file = NULL)

minReq <- list(c("meta_dir","meta_files"), "meta_data_file")

###
### process user input
###
if(!interactive()){
    ## if manifest of parameters is provided, any additional parameters given on command line
    ## will override those provided in manifest
    suppressMessages(library(R.utils))
    args <- processCMD(commandArgs(asValue=TRUE), dflt, minReq, usage)
} else {
    args <- list(manifest       = "input/config/study_config.yaml",
                 meta_data_file = "all_meta_data.rda")
    args <- resolveConfig(dflt, read_yaml(args$manifest), args)
}


###
### validate user input 
###
if(!any(is.null(args$meta_dir), is.null(args$meta_files))){
    err <- "Values were found for both 'meta_dir' and 'meta_files'. Please provide one or the other."
    log_error(err)
    stop(err)
}

if(!is.null(args$meta_files)){ 
    metaFiles <- args$meta_files %>% strsplit(',') %>% unlist
} else { 
    metaFiles <- getFiles(path = args$meta_dir, pattern = ".xlsx")
}
log_debug("Meta files:")
tmp <- sapply(metaFiles, function(x) log_debug(paste(" ", x))) 

oFile <- args$meta_data_file
oType <- validFileType(oFile)
log_debug("Compiled & flattened meta data saved in file:")
log_debug(paste(" ", oFile))


###
### process & save meta data
###
annot <- loadStudyAnnotations(metaFiles = metaFiles)
if(oType == 'rda'){
    saveRDS(annot, args$meta_data_file)
} else { 
    openxlsx::write.xlsx(annot, args$meta_data_file)
}

