suppressMessages(library(funr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(logger))
suppressMessages(library(parallel))
suppressMessages(library(R.utils))
suppressMessages(library(sp))
suppressMessages(library(tidyverse))
suppressMessages(library(xlsx))
suppressMessages(library(yaml))

log_threshold(DEBUG)

sdir <- dirname(get_script_path())
log_debug(paste0("loading source files from: ",sdir))

source(file.path(sdir, "R/constants.R"))
source(file.path(sdir, "R/boundaries.R"))
source(file.path(sdir, "R/exclusions.R"))
source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/load_data.R"))
source(file.path(sdir, "R/script_utils.R"))
source(file.path(sdir, "R/spatial_utils.R"))

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript mark_exclusions.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --control_marker           name of marker used as control; all cells negative for 
                                       this marker are removed
            --drift_summary_dir        path to directory containing drift summary files
            --drift_threshold          maximum allowed percentage drift/loss to include in analyses
            --exclusion_data_dir       output directory; path to processed, exclusion-marked  
                                       RDA file of formatted halo object data
            --halo_boundaries_rda_dir  path to RDA files of halo boundaries (exclusions and interface) 
            --meta_dir                 path to meta files in XLSX format, required IF 
                                       meta_data_file is NULL
            --pad                      numeric; width in microns of pad around each FOV; cells falling 
                                       outside this border will be excluded
            --raw_data_dir             directory containing RDA files of formatted marker-level 
                                       Halo data; required when raw_data_files is NULL

          [OPTIONAL]
            --cell_dive_id      a single CellDive_ID to be marked for exclusion in the associated
                                data file in --raw_data_dir; when NULL (default) all files in that
                                directory will be processed
            --manifest          YAML file containing one or more parameter; NOTE: arguments on command 
                                line override manifest arguments!!!         
            --meta_files        comma-delimited list of meta files 
            --meta_data_file    RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                                meta_files are NULL
            --number_threads    number of threads to use for parallel processes
        \n"
    )
}

minReq   <- list("control_marker",
                 "drift_summary_dir",
                 "drift_threshold",
                 "exclusion_data_dir",
                 "halo_boundaries_rda_dir",
                 c("meta_dir","meta_files","meta_data_file"),
                 "pad",
                 "raw_data_dir")

defaults <- list(number_threads = 1)

args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(args, names(args)[!names(args) %in% c("no-restore", "slave", "args", "file")])

###
### get halo data & meta data 
###
metaFiles    <- getFiles(path = args$meta_dir, files = args$meta_files, pattern = ".xlsx")
rawDataFiles <- getFiles(path = args$raw_data_dir)
driftFiles   <- getFiles(path = args$drift_summary_dir)
boundFiles   <- getFiles(path = args$halo_boundaries_rda_dir)

sampAnn <- read.xlsx(metaFiles[grepl("FOVs", metaFiles)], 1, check.names = F) %>% 
           as_tibble() %>%
           select(CellDive_ID, FOV_number, FOV_exclusion, Marker_exclusion) 

ids <- unique(sampAnn$CellDive_ID)
if(!is.null(args$cell_dive_id)){ ids <- args$cell_dive_id }

cl <- makeCluster(args$number_threads, type = "FORK", outfile = "")
clusterExport(cl, c("args", "rawDataFiles", "driftFiles", "boundFiles", "sampAnn"), envir = environment())
parLapply(cl, ids, function(id){ 
    log_debug(paste0("Marking exclusions in sample: ", id))
    ptrn <- paste0("^", id, "_") 
    rdf  <- rawDataFiles[grep(ptrn, basename(rawDataFiles))]
    df   <- driftFiles[grep(ptrn, basename(driftFiles))]
    bf   <- boundFiles[grep(ptrn, basename(boundFiles))]

    log_debug(paste0("Loading raw Halo data from file: ", rdf))
    dat    <- loadHaloDataFile(rdf, controlMarker = args$control_marker) %>%
              mutate(CellDive_ID = id)
    log_debug(paste0("Loading drift summary: ", df))
    drift  <- read.csv(df, sep = "\t") %>% as_tibble()
    bounds <- readRDS(bf)
    fovAnn <- sampAnn %>% filter(CellDive_ID == id) 
    
    log_debug(paste0("Adding EXCLUDE column to raw data"))
    exclDat <- markExclusions(dat, drift, fovAnn, bounds, 
                              driftThreshold = as.numeric(args$drift_threshold), 
                              borderPad = as.numeric(args$pad))

    outFile <- file.path(args$exclusion_data_dir, gsub(".rda", "_Excl.rda", basename(rdf)))
    saveRDS(exclDat, outFile) 

    #### TEMPORARY - QC
    log_debug("QCing...")
    old <- readRDS(file.path("preprocessing/03_annotated", paste0(id, "_annotated_cells.rda"))) 
    bbFOV <- getConstantBoundingBox(id)
    pad <- as.numeric(args$pad)
    bb <- list(X0 = bbFOV$X0 + pad/pixel2um,
               X1 = bbFOV$X1 - pad/pixel2um,
               Y0 = bbFOV$Y0 + pad/pixel2um,
               Y1 = bbFOV$Y1 - pad/pixel2um)
    clrs <- c(excB = "orange", glsB = "white", epiB = "lightblue", tumB = "black")

    pdf(file.path("qc", paste0(id, "_exclusion_qc.pdf")), height = 5, width = 10)
    for(fovID in unique(old$FOV_ID)){
        log_debug(paste0("fovID ", fovID))
        fov = gsub("^0", "", unlist(strsplit(fovID, "_"))[3])
        log_debug(paste0("fov ", fov))

        fovOld <- old %>% filter(FOV_ID == fovID) %>% mutate(X = (XMin + XMax)/2, Y = -(YMin + YMax)/2)
        fovB = bounds[[fov]]
        new <- exclDat %>%
               filter(FOV_number == as.numeric(fov), EXCLUDE == "") %>%
               select(FOV_number, UUID, Marker, X, Y, Value, EXCLUDE) %>%
               spread(Marker, Value)

        oldP <- ggplot(data = fovOld, aes(x = X, y = Y)) +
                geom_point(size = 0.2, color = "gray") +
                geom_rect(xmin = bb$X0, xmax = bb$X1, ymin = bb$Y0, ymax = bb$Y1,
                          fill = NA, color = "black", linetype = "dashed")
        newP <- ggplot(data = new, aes(x = X, y = Y)) +
                geom_point(size = 0.2, color = "gray") +
                geom_rect(xmin = bb$X0, xmax = bb$X1, ymin = bb$Y0, ymax = bb$Y1,
                          fill = NA, color = "black", linetype = "dashed")
        if(!is.na(fovB) && length(fovB) > 0){
            oldP <- addBoundaries(oldP, fovB, clrs)
            newP <- addBoundaries(newP, fovB, clrs)
        }
        grid.arrange(oldP, newP, ncol = 2)
    } 
    dev.off()
}) 



