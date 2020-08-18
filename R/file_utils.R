#' Extract valid file extension from path
#'
#' Separate file extension from file path basename. If vector of 
#' valid/expected file extensions is provided and extension is
#' not equal to one of them, an error is thrown.
#'
#' @param path      path containing basename from which file extension 
#'                  will be extracted
#' @param expected  optional character vector of expected file extensions
#'
#' @return if path contains an expected file type, the file extension is returned
validFileType <- function(path, expected = NULL){
    ext <- tools::file_ext(path)
    if(!is.null(expected) && !ext %in% expected){
        msg <- paste0("Unexpected file type: [", ext, "]")
        log_error(msg)
        stop(msg)
    }
    ext
}

#' If given PDF file is not NULL, open it
#' 
#' If given PDF file is not NULL, open it
#' 
#' @param pdfFile     path to file
#' @param pdfHeight   height of PDF in inches; default = 8.5
#' @param pdfWidth    width of PDF in inches; default = 11
#' @param singlePage  logical; when TRUE, a blank page will appear before anything is 
#'                    printed to file (?????)
#'
#' @return nothing
openPDF <- function(pdfFile, pdfHeight=8.5, pdfWidth=11, singlePage=FALSE){
    if(!is.null(pdfFile)){
        pdf(pdfFile, height=pdfHeight, width=pdfWidth, onefile = !singlePage)
    }
}

#' If given PDF is not null, close it
#' 
#' If given PDF is not null, close it
#' 
#' @param pdfFile  path to file to be closed
#'
#' @return nothing
closeOpenPDF <- function(pdfFile){
    if(!is.null(pdfFile)){
        dev.off()
    }
}

#' Get subset of files either from a directory or from a list of files
#' that have names matching a specific pattern
#' 
#' Given either a directory path or a vector of file names, return either
#' all files in the file list or, if a pattern is provided, the subset
#' of files matching that pattern
#'
#' @param path    full path to directory containing files of interest; if 'pattern' is 
#'                provided, return only the files within the directory that match; default: NULL
#' @param files   vector of files; if 'pattern' is provided, return only the files that
#'                match; default = NULL
#' @param pattern character pattern of interest; all files that match pattern will be returned
#'
#' @return  vector of the full set or a subset of given file list
getFiles <- function(path=NULL, files=NULL, pattern=NULL){

    if(is.null(path) && is.null(files)){
        stop("Need either a directory or a vector of files and a pattern to get files.")
    }
    fls <- c()
    if(!is.null(files)){
        fls <- files
    } else if(!is.null(path)){
        fls <- file.path(path, dir(path))
    } else {
        stop("Must provide either path or files argument.")
    }
    if(!is.null(pattern)){
        fls <- fls[grepl(pattern, fls)]
    }
    fls

}

#' Determine whether a file exists and is not empty
#'
#' Determine whether a file exists and is not empty
#' 
#' @param fileName  file name

#' @return logical
fileDone <- function(fileName){
    !is.null(fileName) && length(fileName) > 0 && file.exists(fileName) && file.size(fileName) > 0
}


#' Determine whether a file in a configuration list needs to/should be written
#' 
#' @param lst  list of study parameters
#' @param key  key of file to test

#' @return logical
needToWrite <- function(lst, key){
    is.null(lst[[key]]) || length(lst[[key]]) == 0 || !file.exists(lst[[key]])
}


#' Create a directory recursively without warnings
#' 
#' Create a directory recursively without warnings
#' 
#' @param path    full path of directory to create
#'
#' @return nothing
mkdir <- function(path){
    dir.create(path, recursive = T, showWarnings = F)
}
