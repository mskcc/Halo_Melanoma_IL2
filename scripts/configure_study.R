suppressMessages(library(funr))
suppressMessages(library(R.utils))
suppressMessages(library(yaml))
suppressMessages(library(logger))

sdir <- dirname(get_script_path())
log_debug(paste0("loading source files from: ",sdir))

source(file.path(sdir, "R/config.R"))
source(file.path(sdir, "R/file_utils.R"))
source(file.path(sdir, "R/validation.R"))
source(file.path(sdir, "R/script_utils.R"))

log_threshold(DEBUG)

#####################################
#####        SET UP INPUT       #####
#####################################
usage <- function(){

    cat("\nUsage:  Rscript configure_study.R 
            
          [REQUIRED (may be defined on command line OR in manifest file)] 
            --source_halo_csv_dir              path to CSV files output by Halo
            --source_halo_boundaries_xml_dir   path to subdirectories of Halo XML files of boundary annotations
            --source_halo_image_dir            path to images written by Halo
            --source_meta_dir                  path to XLSX files of all study meta data

          [OPTIONAL]
            --study_config                     YAML file containing custom config
            --plot_config_file                 YAML file containing plot colors to be applied to all analyses
            --statistics_config_file           YAML file containing configuration of statistical analyses & visualizations
            --statistics_questions_file        XLSX file of questions to ask, including filtering criteria
            --statistics_conditions_file       XLSX file of all conditions to be analyzed      
 
        \n"
    )
}

### SET DEFAULTS:
defaults <- list(log_dir = "logs",

                 ### input and data preprocessing
                 config_dir = "input/config",
                 meta_dir = "input/meta",
                 halo_csv_dir = "input/halo_csv",
                 halo_image_dir = "input/halo_images",
                 halo_boundaries_xml_dir ="input/halo_boundaries",

                 raw_data_dir = "preprocessing/00_raw",
                 halo_boundaries_rda_dir = "preprocessing/00_halo_boundaries",
                 drift_summary_dir = "preprocessing/00_drift_summaries",
                 exclusion_data_dir = "preprocessing/01_exclusions",
                 data_dir = "preprocessing/02_rethresholded",
                 cell_data_dir = "preprocessing/03_annotated",
                 band_dir = "preprocessing/03_infiltration_band_data",
                 neighborhood_dir = "preprocessing/04_neighborhoods/all_vs_all",
                 neighborhood_data_dir = "preprocessing/04_neighborhoods/formatted",
                 microenvironment_dir = "prepocessing/05_microenvironments",
                 tme_by_sample_dir = "prepocessing/05_microenvironments/sample_status", 
                 tme_by_cell_dir = "prepocessing/05_microenvironments/cell_status",
                 qc_dir = "qc",

                 statistics_conditions_index = "processed/statistics/conditions_index.xlsx",
                 statistics_conditions_file = "input/config/stats_conditions.xlsx",
                 statistics_questions_file = "input/config/stats_questions.xlsx",

                 ### compiled intermediate data
                 annotated_cells_file = "annotated_cells.rda",
                 halo_boundaries_file = "halo_boundaries.rda",
                 meta_data_file = "all_meta_data.rda",

                 ### calculations
                 metrics_dir = "processed/metrics",
                 fov_metrics_dir = "processed/metrics/fovs",
                 fov_area_dir = "processed/metrics/fovs/areas",
                 infiltration_metrics_dir = "processed/metrics/infiltration",
                 band_area_dir = "processed/metrics/infiltration/areas",
                 neighborhood_metrics_dir = "processed/metrics/neighborhoods",

                 ### deliverables
                 counts_dir = "results/counts",
                 cell_type_marker_counts_file = "results/counts/cell_type_marker_combo_counts.xlsx",
                 statistics_dir = "results/statistics",
                 statistics_tables_dir = "results/statistics/tables",
                 statistics_overview_dir = "results/statistics/figures/overviews",
                 statistics_detail_dir = "results/statistics/figures/detail",
                 statistics_misc_dir = "results/statistics/figures/misc",

                 ### other config files
                 annotation_config_file = "~",
                 statistics_config_file = "~",
                 plot_config_file = "~",
                 detail_figure_config_file = "~",
                 detail_figure_conditions_file = "~",
                 condition_detail_config = "~",
                 infiltration_density_config = "~",
                 
                 ### parameters
                 pad = 20,
                 drift_threshold = 0.1,
                 max_g = 5.0,
                 band_width = 10,
                 max_distance_from_interface = 360,
                 number_threads = 6.0,
                 debug = "yes"
            )

minReq <- c("source_halo_csv_dir",
            "source_halo_boundaries_xml_dir",
            "source_halo_image_dir",
            "source_meta_dir")

used <- c(minReq, names(defaults))
cfg  <- processCMD(commandArgs(asValue=TRUE), defaults, minReq, usage)

logParams(cfg, names(cfg)[!names(cfg) %in% c("no-restore", "slave", "file", "args")])
cat("Directory structure will be created according to the configuration above.")
cat("\n Type 'stop' to write config to file for editing. Press return to proceed: ") 

proceed <- readLines("stdin", n=1)
if(tolower(proceed) == "stop"){
    cat("\nWriting template config to file
        
           TEMPLATE_CONFIG.yaml

           Please modify as necessary, rename and rerun configuration with '--study_config [config_file]'\n\n")
    write_yaml(cfg, "TEMPLATE_CONFIG.yaml")
    q()
}

log_info("Creating directories...")
dir_setup(cfg)

log_info("Linking remote input files to local directories...")
link_input(cfg, defaults)

## clean up
cfg$slave <- cfg$`no-restore` <- cfg$file <- cfg$args <- NULL

## save config
cfg_file <- file.path(cfg$config_dir, "study_config.yaml")
log_info(paste0("Saving all study configuration to: ", cfg_file))
write_yaml(cfg, cfg_file) 

## validate
log_info("Validating study config...")
scValid <- validateStudyConfig(cfg)

## copy pipeline snakefile to current directory to be run after configuration
log_info("Saving pipeline script: Snakefile")
tmp <- file.copy(file.path(sdir,"../Snakefile"), getwd())
log_info("Configuration successful.")
