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
            --figure_config_file             YAML file of figure configuration including
                                             facet order & biological filters
            --plot_color_file                YAML file containing all color assignments 
            --metrics_dir                    root directory of all precalculated metrics
            --statistics_detail_dir          output directory where figure(s) will be saved
            --statistics_conditions_file     YAML file containing all 'conditions' or 
                                             cell states to analyze 
            --statistics_conditions_index    output XLSX file where indexed conditions
                                             will be saved
            --detail_figure_conditions_file  XLSX file containing a column for each question, 
                                             each column containing any number of Cell State IDs
                                             that should be included in the figure
            --meta_dir                       directory containing meta data XLSX files
            --statistics_questions_file      XLSX file describing questions by indicating how to
                                             filter data to form two sample groups for each question

          [OPTIONAL]
            --question_number         question number of comparison to be plotted
            --plot_calculation        calculation to be plotted (densities or fractions)
            --manifest                YAML file containing one or more parameter; NOTE: 
                                      arguments on command line override manifest arguments!!!         
        \n"
    )
}

## names of required args
minReq <- c("statistics_conditions_file",
            "figure_config_file",
            "plot_color_file",
            "metrics_dir",
            "meta_dir",
            "statistics_conditions_index",
            "statistics_questions_file")

used <- names(minReq)
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
