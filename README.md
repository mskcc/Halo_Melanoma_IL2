# Analysis pipeline for Melanoma IL2 study

This is the first pass at a pipeline to analyze Melanoma IL2 sample data from Halo.

## Curate input data

Input includes paths to data output directly from Halo (CSV, images, XML boundaries) and manually curated study meta data.

#### Meta data files to be created are:

    - [StudyName]_Patients.xlsx
    - [StudyName]_Samples.xlsx
    - [StudyName]_FOVs.xlsx
    - [StudyName]_CellTypes.xlsx
    - [StudyName]_Markers.xlsx
    - [StudyName]_Files.xlsx

See [study_meta_data.md](docs/study_meta_data.md) for instructions on how to create each of these files.

## Configuration

Run script to create project directory structure and generate YAML file of configuration to be used throughout the pipeline.

Create project directory and run configuration
<b>TO DO: update this (and config script) to parse [StudyName]_Files.xlsx instead of taking in each path </b>
```
mkdir [project_dir]
cd [project_dir]
Rscript [src_dir]/scripts/configure_study.R \
    --source_halo_csv_dir             /path/to/halo/csv/files \
    --source_halo_image_dir           /path/to/halo/image/files \
    --source_halo_boundaries_xml_dir  /path/to/halo/boundary/xml/files \
    --source_meta_dir                 /path/to/meta/data
```

In addition to creating the directory structure for your study, the configuration will link all source files to the newly created project directory and set the configuration to point to the local directories. Study configuration will be validated and written to:

    **input/config/study_config.yaml**

**NOTE: One of the most important checks done during configuration is the list of samples to be analyzed. PLEASE review warnings carefully, fix errors and check the 'SAMPLES' parameter in this config file before proceeding!!**

## Notes on R

Currently using local R installation:

    /home/byrne/R/R-3.6.1

Check docs for dependencies

## Input data validation

Check all manually-created meta files for errors. You may either pass to the script the entire configuration or just the directory of meta data files.

```
Rscript [src_dir]/scripts/validate_input.R --manifest input/config/study_config.yaml
```
or if configuration is not yet set up/validated:
```
Rscript [src_dir]/scripts/validate_input.R --meta_dir /path/to/meta/files
```

## Run pipeline

Copy Snakefile from source code directory to project directory.

```
snakemake -S Snakefile all
```
