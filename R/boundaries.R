#' Recursively process an XML node for one FOV
#'
#' Visit each node in an object of class XMLNode and gather information
#' to print one line for each vertex which represents the location of a cell
#' on the slide [NOTE: DOUBLE CHECK THIS FOR ACCURACY] 
#' 
#' Output is a tab-delimited file containing: 
#'     1)  FOV name
#'     2)  Annotation number
#'     3)  Line color
#'     4)  Name
#'     5)  Visible
#'     6)  Region number
#'     7)  Region type
#'     8)  HasEndcaps
#'     9)  NegativeROA
#'     10) X coordinate
#'     11) Y coordinate ***(negative of what is output by Halo)***
#'     12) Region according to color (Tumor, Exclusion, Epidermis)
#' 
#' @param node           object of class XMLNode containing Halo annotations for one FOV
#' @param fovName        unique identifier of FOV
#' @param outFile        output file
#' @param LineColor      code for line color of the annotation ["655280"|"65535"|"255"]
#' @param Name           name of the annotation
#' @param Visible        logical indicating whether annotation is visible [????]
#' @param RegionType     region type *Note: any region with type other than Polygon is skipped
#' @param HasEndcaps     [???]
#' @param NegativeROA    [???]
#' @param boundaryColors a list including hex colors for Tum, Exc, and Epi boundaries; 
#'                       default=list("65280"="Tum", "65535"="Exc", "255"="Epi")
#' @return nothing
parseAnnotationXML <- function(node, fovName=NULL, outFile=NULL, LineColor="",
                               Name="", Visible="", RegionType="", HasEndcaps=0,
                               NegativeROA=0, boundaryColors=NULL){
  ## use global variables in order to keep track of annotation and region nums
  ## while traversing xml doc; remove them later
  annotationNum <<- 0
  regionNum <<- 0

  processNode <- function(node, fovName=NULL, outFile=NULL,
                          LineColor="", Name="", Visible="", RegionType="", HasEndcaps=0,
                          NegativeROA=0){

    if(is.null(node) | length(node) == 0){
        return()
    }

    nodeName <- xmlName(node)
    nodeAttrs <- xmlAttrs(node)

    if(nodeName == "Annotation"){
        annotationNum <<- annotationNum + 1
        #log_debug(paste0("AnnotNum: ",annotationNum))
        Name <- xmlGetAttr(node, "Name")
        Visible <- xmlGetAttr(node, "Visible")
        LineColor <- xmlGetAttr(node, "LineColor")
    } else if(nodeName == "Region"){
        regionNum <<- regionNum + 1
        #log_debug(paste0("RegionNum: ",regionNum))
        HasEndcaps <- xmlGetAttr(node, "HasEndcaps")
        NegativeROA <- xmlGetAttr(node, "NegativeROA")
        RegionType <- xmlGetAttr(node, "Type")
    } else if(nodeName == "V"){
        X = xmlGetAttr(node, "X")
        Y = -as.numeric(xmlGetAttr(node, "Y"))
        if(RegionType == "Polygon"){
            line <- paste(fovName, annotationNum, LineColor, Name, Visible,
                          regionNum, RegionType, HasEndcaps, NegativeROA, X, Y,
                          boundaryColors[[LineColor]], sep="\t")
            write(line, file=outFile, append=TRUE)
        }
        return()
    }

    if(xmlSize(node) > 0){
       for(i in 1:xmlSize(node)){
           processNode(node[[i]], fovName=fovName, LineColor=LineColor,
                       Name=Name, Visible=Visible, RegionType=RegionType,
                       HasEndcaps=HasEndcaps, NegativeROA=NegativeROA,
                       outFile=outFile)
       }
    }
  }

  res <- processNode(node, fovName=fovName, outFile=outFile,
           LineColor=LineColor, Name=Name, Visible=Visible, RegionType=RegionType,
           HasEndcaps=HasEndcaps, NegativeROA=NegativeROA)

  suppressWarnings(rm(annotationNum))
  suppressWarnings(rm(regionNum))
}


#' Parse boundary annotations of a single FOV in a single XML file
#' 
#' Parse XML file containing annotations for one FOV. Write tab-delimited file
#' of select information (see ?parseAnnotationXML for details)
#' 
#' @param annoteFile      XML file containing boundary annotations for a single FOV
#' @param boundaryColors  named list of integers, each representing the color used to draw
#'                        a single boundary type and named by the boundary type 
#'                        (Epi|Tum|Exc|Gls) 
#' @export
readHaloAnnotations <- function(annoteFile,boundaryColors=NULL,boundaryReassignmentFile=NULL){

    tfile <- tempfile(tmpdir=".")
    header = paste("File","AnnotationNum","LineColor","Name","Visible","RegionNum","Type",
                   "HasEndcaps","NegativeROA","X","Y","RegionCode",sep="\t")
    write(header, file=tfile, append=FALSE)

    rootNode <- xmlRoot(xmlParse(annoteFile))
    fovName <- gsub("\\.annotations","",basename(annoteFile))
    parseAnnotationXML(rootNode, fovName=fovName, outFile=tfile, boundaryColors=boundaryColors)

    aa <- read_tsv(tfile, col_types="ciccciciiiic")
    fixFile <- boundaryReassignmentFile
    if(!is.null(fixFile) && file.exists(fixFile)){
        log_debug("Reassign file=%s",fixFile)
        fix <- read_csv(fixFile, col_types="ccic")
        ff <- fix %>% filter(Spot==aa$File[1])
        if(nrow(ff)>0){
            for(kk in seq(nrow(ff))){
                aa$RegionCode[ aa$RegionNum==ff$RegionNum[kk] ] <- ff$RegionCode[kk]
            }
        }
    }

    file.remove(tfile)
    split(aa,aa$RegionNum)
}

#' Remove exclusion boundaries that are contained in another one
#' 
#' Remove data from boundaries table that represent exclusion boundaries
#' that are completely surrounded by another exclusion boundary
#' 
#' @param boundaries            table generated by readHaloAnnotations
#'                             
cleanBoundaries <- function(boundaries){

    regionTable <- boundaries %>% bind_rows %>% count(RegionNum,RegionCode)
    excB <- boundaries[regionTable$RegionCode=="Exc"]
    epiB <- boundaries[regionTable$RegionCode=="Epi"]
    glsB <- boundaries[regionTable$RegionCode=="Gls"]
    tumB <- boundaries[regionTable$RegionCode=="Tum"]

    excList <- list(excB=excB, epiB=epiB, glsB=glsB)
    for(x in names(excList)){
        ex <- excList[[x]]
        if(!is.null(ex) && length(ex) > 1){
            ## make list of spatial polygons, one element for each exclusion boundary
            spExcB <- seq(ex) %>% map(function(b){boundaryToSpatialPolygon(ex[[b]],b)})

            ## are any contained within another?
            containedBoundary <- rep(FALSE,length(ex))
            for(i in seq(length(ex))){
                containedBoundary[i] <- spExcB[-i] %>% map(gContains,spExcB[[i]]) %>% unlist %>% any
            }
            if(any(containedBoundary)){
                log_debug(paste0("Eliminating ", length(which(containedBoundary)), " contained boundaries."))
            }
            excList[[x]] <- ex[!containedBoundary]
        }
    }
    excList[["tumB"]] <- tumB
    return(excList)
}


getAllHaloBoundaries <- function(samp, boundaryDirs, boundaryColors=NULL){

    allFiles <- tibble()
    for(bd in boundaryDirs){
        fs <- getFiles(path = bd, pattern = paste0("^", samp, "_"))
        if(length(fs) == 0){ next } 
        for(f in fs){
            fov <- gsub(".annotations", "", gsub(".*Spot","",f))
            allFiles <- allFiles %>% 
                        bind_rows(tibble(File=f, Sample=samp, FOV=fov, BoundaryType=basename(bd)))
        }
    }
    if(nrow(allFiles) == 0){ 
        log_warn(paste0("No boundary files found for sample ", samp, " in directories."))
        return(NULL)
    }

    btCodes <- list(interface = 'tumB', exclusions = 'excB', epidermis = 'epiB', glass = 'glsB')

    allAnnotations <- list()
    for(fov in unique(allFiles$FOV)){
        clrs <- boundaryColors
        fovBoundaries <- list()
        aFiles <- allFiles %>% filter(FOV == fov) %>% select(File,BoundaryType)

        for(x in 1:nrow(aFiles)){
            clrs <- boundaryColors
            af <- aFiles$File[x]
            log_debug(paste0("getting boundaries from file ",af))

            ##### FIX FOR FOR COHORT 1 ONLY
            boundaryType <- aFiles$BoundaryType[x]
            bt <- btCodes[[boundaryType]]
            if(samp %in% c("Untreated", "PR", "CR")){
                clrs <- list('65280'='Epi','65535'='Exc','255'='Gls')
                if(boundaryType == "interface"){
                    if(samp == "Untreated"){
                        clrs <- list('65280'='Tum','65535'='IGNORE', '255'='IGNORE')
                    } else if(samp == "PR"){
                        clrs <- list('65535'='IGNORE', '65280'='Tum')
                        if(as.numeric(fov) %in% c(3,13,14,21)){
                            clrs <- list('65535'='Tum', '65280'='Exc', '255'='IGNORE')
                        }
                    }
                }
            }
            ##########################           

            boundaries <- readHaloAnnotations(af, boundaryColors=clrs) %>%
                          cleanBoundaries()
            fovBoundaries[[bt]] <- boundaries[[bt]][!is.na(names(boundaries[[bt]]))]
        }
        allAnnotations[[fov]] <- fovBoundaries
    }

    allAnnotations

}

summarizeHaloBoundaries <- function(boundary_rda_dir, cellDiveIDs){
    summ <- tibble()
    for(id in ids){
        bounds <- readRDS(file.path(args$halo_boundaries_rda_dir, paste0(id, "_halo_boundaries.rda")))
        for(fov in names(bounds)){
            fovSumm <- tibble(CellDive_ID = id,
                              FOV_number = as.integer(fov),
                              exclusions = length(bounds[[fov]][["excB"]]) +
                                           length(bounds[[fov]][["glsB"]]),
                              epidermis = length(bounds[[fov]][["epiB"]]),
                              interface = length(bounds[[fov]][["tumB"]]))
            summ <- summ %>% bind_rows(fovSumm)
        }
    }
    summ[summ == 0] <- NA
    summ
}

boundaryDiscrepancies <- function(boundaryDir, fovAnn){
    boundarySummary <- summarizeHaloBoundaries(boundaryDir, ids)

    fovAnn %>%
    select(CellDive_ID, FOV_number, dplyr::matches("Num_")) %>%
    full_join(boundarySummary, by = intersect(names(.), names(boundarySummary))) %>%
    select(CellDive_ID, FOV_number, Num_manual_exclusions, exclusions, dplyr::matches("epidermis"),
           dplyr::matches("interface")) %>%
    gather(3:ncol(.), key = "boundary", value = "Count") %>%
    mutate(BoundaryType = gsub("_.*", "", gsub("Num_", "", boundary))) %>%
    mutate(BoundaryType = ifelse(BoundaryType == "manual", "exclusions", BoundaryType),
           Cat = ifelse(grepl("Num", boundary), "Expected", "Observed")) %>%
    select(-boundary) %>%
    spread(Cat, Count) %>%
    filter(Expected != Observed)
}


addBoundaries <- function(p, bounds, clrs){

    if(!is.null(bounds$excB)){
        for(reg in names(bounds$excB)){
            p <- p + geom_point(data = bounds$excB[[reg]], aes(x = X, y = Y), 
                                color = clrs['excB'], size = 0.2)
        }
    }
    if(!is.null(bounds$tumB)){
        for(reg in names(bounds$tumB)){
            p <- p + geom_point(data = bounds$tumB[[reg]], aes(x = X, y = Y), 
                                color = clrs['tumB'], size = 0.2)
        }
    }
    if(!is.null(bounds$epiB)){
        for(reg in names(bounds$epiB)){
            p <- p + geom_point(data = bounds$epiB[[reg]], aes(x = X, y = Y), 
                                color = clrs['epiB'], size = 0.2)
        }
    }
    if(!is.null(bounds$glsB)){
        for(reg in names(bounds$glsB)){
            p <- p + geom_point(data = bounds$glsB[[reg]], aes(x = X, y = Y), 
                                color = clrs['glsB'], size = 0.2)
        }
    }
    return(p)
}


qcHaloBoundaries <- function(cellDiveIDs, rawDataDir, tumorMarker, boundaryDir, allDiscrep, qcDir, threads = 1){

    clrs <- c(excB = "orange", glsB = "white", epiB = "lightblue", tumB = "black")

    cl <- makeCluster(threads, type = "FORK", outfile = "")
    clusterExport(cl, c("rawDataDir", "tumorMarker", "boundaryDir",
                        "allDiscrep", "qcDir"), envir = environment())
    parLapply(cl, cellDiveIDs, function(id){

        print(id)
        rawFile <- getFiles(path = rawDataDir, pattern = paste0("^", id, "_"))
        cells <- readRDS(rawFile) %>%
                 mutate(Tumor = ifelse(Marker == tumorMarker & ValueType == "Positive" & Value == 1, "YES", NA),
                        X = (XMin + XMax)/2, Y = -(YMin + YMax)/2) %>%
                 select(FOV_number = SPOT, UUID, X, Y, Tumor) %>%
                 unique()
        fovs <- unique(cells$FOV_number)
        tumorCells <- cells %>% filter(Tumor == "YES")
        cells      <- cells %>% filter(is.na(Tumor))

        bounds <- readRDS(file.path(boundaryDir, paste0(id, "_halo_boundaries.rda")))

        pdf(file.path(qcDir, paste0(id, "_halo_boundaries_qc.pdf")), height = 5, width = 7)
        for(fov in sort(fovs)){
            fovDat <- cells %>% filter(FOV_number == fov)
            fovTum <- tumorCells %>% filter(FOV_number == fov)

            p <- ggplot(data = fovDat, aes(x = X, y = Y)) +
                 geom_point(color = "gray", size = 0.2) +
                 geom_point(data = fovTum, color = "red", size = 0.2) +
                 labs(title = paste(id, "FOV", fov)) +
                 scale_color_manual(name = "Boundary", values = clrs, labels = names(clrs))

            if(fov %in% names(bounds)){
                p <- addBoundaries(p, bounds[[as.character(fov)]], clrs)
            }

            fovDiscrep <- allDiscrep %>% filter(CellDive_ID == id, FOV_number == fov)
            if(nrow(fovDiscrep) > 0){
                tbl <- tableGrob(fovDiscrep, rows = NULL)
                grid.arrange(p, tbl, nrow = 2, heights = c(3.75, 0.75))
            } else {
                print(p)
            }
        }
        dev.off()
    })
}

