suppressMessages(library(logger))
suppressMessages(library(R.utils))
suppressMessages(library(cowplot))
suppressMessages(library(funr))
suppressMessages(library(grid))
suppressMessages(library(parallel))
suppressMessages(library(yaml))
suppressMessages(library(tidyverse))
suppressMessages(library(xlsx))

log_threshold(DEBUG)

#####################################
#####        SOURCE CODE        #####
#####################################

sdir <- dirname(get_script_path())
log_debug(paste0("loading source files from: ",sdir))

source(file.path(sdir, "R/cell_annotation.R"))
source(file.path(sdir, "R/data_utils.R"))
source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/filter_data.R"))
source(file.path(sdir, "R/format_data.R"))
source(file.path(sdir, "R/format_plot_data.R"))
source(file.path(sdir, "R/load_data.R"))
source(file.path(sdir, "R/metrics.R"))
source(file.path(sdir, "R/plot_utils.R"))
source(file.path(sdir, "R/script_utils.R"))
source(file.path(sdir, "R/statistics.R"))
source(file.path(sdir, "R/stats_detail_plot.R"))

#####################################
#####       GET USER INPUT      #####
#####################################

usage <- function(){
    cat("\nUsage:  Rscript calculate_metrics.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)]
            --band_dir                       path to RDA files, each containing a table of band 
                                             assignments and band areas for each cell
            --cell_data_dir                  path to RDA files, each containing a single 
                                             table where rows are cells and columns are 
                                             all data for that single cell
            --figure_config_file             YAML file of figure configuration including
                                             facet order & biological filters
            --fov_area_dir                   path to RDA files, each containing table of FOVs and 
                                             total FOV area for all FOVs in a single sample 
            --meta_dir                       path to meta files in XLSX format
            --metrics_dir                    root directory of all precalculated metrics
            --plot_color_file                YAML file containing all color assignments 
            --statistics_conditions_file     XLSX file listing all cell states/conditions to 
                                             compare between two sample groups
            --statistics_conditions_index    XLSX file with pre-indexed cell states/conditions
            --statistics_questions_file      XLSX file outlining all questions/comparisons for 
                                             which stats should be run
            --statistics_tables_dir          output directory where XLSX files of results should 
                                             be written

          [OPTIONAL]
            --detail_figure_conditions_file  XLSX file containing a column for each question, 
                                             each column containing any number of Cell State IDs
                                             that should be included in the figure
            --neighborhood_counts_dir        directory of RDA files containing formatted macrophage 
                                             neighborhood counts 
            --neighborhood_dir               directory containing RDA files of all pairwise distances 
                                             between cells, at least those <= 30 microns. generally 
                                             each file will contain a table of all 'center' cells of
                                             a certain type, but this is not a requirement as all files
                                             will be loaded together; required if question data is restricted
                                             to tumor or macrophage neighborhoods
            --question_number                question number of comparison to be plotted
            --plot_calculation               calculation to be plotted (densities or fractions)
            --manifest                       YAML file containing one or more parameter; NOTE: 
                                             arguments on command line override manifest arguments!!!         
            --tme_by_cell_dir                directory containing RDA files of tumor microenvironment 
                                             assignments; i.e., each cell (see docs for details); NOTE:
                                             not required if cell_dive_id specified and sample is NOT a tumor
        \n"
    )
}

## names of required args
minReq <- c("band_dir", "cell_data_dir", "figure_config_file",
            "fov_area_dir", "meta_dir", "metrics_dir", "plot_color_file",
            "statistics_conditions_file", "statistics_conditions_index",
            "statistics_questions_file", "statistics_tables_dir") 

used <- c(names(minReq), "detail_figure_conditions_file", "neighborhood_counts_dir",
          "neighborhood_dir", "question_number", "plot_calculation", "manifest",
          "tme_by_cell_dir")
defaults <- list()

if(!interactive()){
    suppressMessages(library(R.utils))
    args <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)
} else {
    args <- processCMD(list(manifest = "input/config/study_config.yaml",
                       figure_config_file = "input/config/stats_detail_figure_config.yaml",
                       plot_color_file = "input/config/global_plot_colors.yaml"),
                       defaults, minReq, usage)
}

###################################
##   CONFIGURE & INITIALIZE DATA ##
###################################
cfg <- resolveConfig(args, read_yaml(args$plot_color_file), read_yaml(args$figure_config_file))

logParams(cfg, used)

loadGlobalStudyData(cfg, all = T)

if(!is.null(cfg$question_number)){ 
    allQuestions <- allQuestions[names(allQuestions) == cfg$question_number] 
}
if(!is.null(cfg$plot_calculation)){
    cfg$calculations <- cfg$calculations[names(cfg$calculations) == cfg$plot_calculation]
}

calcUnit      <- "FOV_ID"
calcInfo <- cfg$calculations
statsFiles <- file.path(cfg$statistics_tables_dir, dir(cfg$statistics_tables_dir))

rel_widths <- list("fractions" = c(1.5, 0.7, 0.55),
                   "densities" = c(1, 0.7, 0.55))

###################################
##              PLOT             ##
###################################
for(calc in names(calcInfo)){

    conds2plot <- list()
    if(!is.null(cfg$detail_figure_conditions_file) && file.exists(cfg$detail_figure_conditions)){
        conds2plot <- read.xlsx(cfg$detail_figure_conditions_file, sheetName = calc, check.names=F)
        conds2plot <- conds2plot[,which(colnames(conds2plot) != "NA")] %>% as_tibble()
    }

    for(q in names(allQuestions)){
        log_info(paste0(q, ": Formatting data for plotting ", toupper(calc), "..."))
        ids <- NULL
        if(q %in% names(conds2plot)){
            ids <- conds2plot[[q]][!is.na(conds2plot[[q]])] %>% as.numeric
            log_debug(paste0("Including ", length(ids), " conditions specified in file ",
                             cfg$detail_figure_conditions_file))
        }

        plotDat <- getStatsDetailData(statsFiles[grepl(paste0(q, ".xlsx"), basename(statsFiles))], 
                                      sampAnn, 
                                      annCells, 
                                      allQuestions[[q]], 
                                      conds, 
                                      cellTypes, 
                                      analysisList[names(analysisList) == calc], 
                                      markers, 
                                      cfg$metrics_dir,
                                      sheetName      = "all_fractions_and_densities",
                                      calcUnit       = "FOV_ID", 
                                      calculation    = calc, 
                                      calcColumn     = calcInfo[[calc]]$calc_column, 
                                      cellStateIDs   = ids,
                                      statsFilter    = cfg$bio_filter, 
                                      facetOrder     = cfg$facet_order, 
                                      popOrder       = facetOrder$Tag, 
                                      orderBy        = calcInfo[[calc]]$effect_col, 
                                      facets         = calcInfo[[calc]]$y_facets, 
                                      nbhdCounts     = nbhdCounts, 
                                      tumorNbhdCells = tumorNbhdCells)

        openxlsx::write.xlsx(plotDat, file.path(cfg$statistics_detail_dir, paste0(q, calc, "detail.xlsx", sep="_"))) 

        
        clrKeys     <- gsub(" \\(.*", "", levels(plotDat$fovVals$GroupLabel))
        clrs        <- cfg$clinical_colors[clrKeys]
        names(clrs) <- levels(plotDat$fovVals$GroupLabel)

        pdfFile   <- file.path(cfg$statistics_detail_dir, paste0(q,"_",calc,"_detail.pdf"))
        pdfHeight <- (0.10 * nrow(plotDat$labels)) + 2.5
        pdfWidth  <- ifelse(is.null(cfg$pdfWidth), 10, as.numeric(cfg$pdfWidth))

        log_info("  Generating figure...")
        tryCatch({

            p <- plotQuestionResults(plotDat$fovVals, 
                                     plotDat$stats, 
                                     plotDat$labels, 
                                     clrs, 
                                     "GroupLabel",
                                     calcType    = calc,
                                     cellTypes   = cellTypes,
                                     xVar        = calcInfo[[calc]]$calc_column,
                                     yVar        = "y",
                                     yNudge      = 0.1,
                                     facetY      = calcInfo[[calc]]$y_facets,
                                     fontsize    = 12,
                                     spacerColor = "#e0e0e0",
                                     rel_widths  = rel_widths[[calc]])

            log_info(paste0("  Saving figure to file: ", pdfFile))
            ggsave(p, height = pdfHeight, width = pdfWidth, units = "in",
                      device = cairo_pdf, filename = pdfFile)

          }, error = function(e){
                print(e)
        })
    }
}
