options(dplyr.summarise.inform = FALSE)

srcRoot <- "/home/byrne/halo/dev/melanoma"

suppressMessages(library(assertthat))
suppressMessages(library(contoureR))
suppressMessages(library(cowplot))
suppressMessages(library(digest))
suppressMessages(library(effsize))
suppressMessages(library(egg))
suppressMessages(library(logger))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(kimisc))
suppressMessages(library(lemon))
suppressMessages(library(magrittr))
suppressMessages(library(parallel))
suppressMessages(library(pheatmap))
suppressMessages(library(plotrix))
suppressMessages(library(plyr))  ## load BEFORE tidyverse
suppressMessages(library(randtoolbox))
suppressMessages(library(raster))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reshape))
suppressMessages(library(rgeos))
suppressMessages(library(rJava))
suppressMessages(library(rlang)) ## load BEFORE tidyverse
suppressMessages(library(scales))
suppressMessages(library(SearchTrees))
suppressMessages(library(sp))
suppressMessages(library(tidyverse))
suppressMessages(library(testthat))
suppressMessages(library(tools))
suppressMessages(library(xlsx))
suppressMessages(library(xlsxjars))
suppressMessages(library(XML))
suppressMessages(library(yaml))

sourceDir <- file.path(srcRoot, "R") 
for(f in file.path(sourceDir, dir(sourceDir))){ 
    if(!is.dir(f)){ source(f) } 
} 

#dfltDir <- file.path(srcRoot, "defaults")
#studyDfltFile <- file.path(dfltDir, "study_params.yaml")
#defaults <- read_yaml(studyDfltFile)
