#' Get list of files in a directory with full path
#' 
#' Get list of files in a directory with full path. Optionally,
#' provide a regex pattern to return only files matching the pattern.
#' If not provided, all files in directory will be returned
#'
#' @param path     directory to list
#' @param pattern  [optional]; character string containing a regex pattern
#'                 used to return only certain files
#'
#' @return vector of full file paths
ls_full_path <- function(path, pattern = NULL){

    fp <- file.path(path, dir(path))

    if(!is.null(pattern)){
        fp <- fp[grep(pattern, fp)]
    }

    fp
}

#' Create directory structure 
#' 
#' Create directory structure for input, temporary and final results
#' files
#' 
#' @param cfg  config in list form; character values for all key names 
#'             ending with "_dir" and NOT beginning with "source_" are assumed
#'             to be full paths to directories to be created
#'
#' @return nothing
dir_setup <- function(cfg){
    
    dirs <- names(cfg)[grepl("_dir$", names(cfg))]

    for(dr in dirs){
        if(!grepl("^/", dr) && !grepl("^source_", dr)){
            mkdir(cfg[[dr]])
        }
    }
}

#' Create symbolic links for all source data in local project directory
#' 
#' Link all source files in local directory
#' 
#' @param cfg  config in list form; character values for all key names
#'             starting with "source_" are assumed to be directories containing
#'             files to be linked
#' @param defaults  default config in list form containing keys corresponding to 
#'                  each 'source_*' key in cfg (e.g., 'source_meta_data' in cfg
#'                  should have 'meta_data' in defaults; this is the dir where
#'                  links will live
#'
#' @return nothing
link_input <- function(cfg, defaults){

    src_dirs <- names(cfg)[grep("source_", names(cfg))]

    for(sd in src_dirs){
        loc <- defaults[[gsub("source_", "", sd)]]
        rem <- ls_full_path(cfg[[sd]])
        for(rf in rem){
            suppressMessages(file.symlink(rf, loc))
        }
    } 
}
