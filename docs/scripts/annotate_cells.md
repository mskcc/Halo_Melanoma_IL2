## annotate_cells.R

Usage:
```
Rscript annotate_cells.R 
            
    [REQUIRED (may be defined on command line OR in manifest file)] 
        --annotated_cells_file  path to RDA file (output) that will contain a 
                                single table where rows are cells and columns 
                                are all data for that single cell
        --control_marker        name of marker used as control; all cells negative for 
                                this marker are removed
        --data_dir              path to processed, exclusion-marked RDA files of 
                                formatted halo object data; required if annotated_cells_file 
                                does NOT exist and if data_files is NULL 
        --meta_dir              path to meta files in XLSX format, required IF meta_data_file 
                                is NULL

    [OPTIONAL]
        --force_reannotation  when TRUE, will start with rethresholded halo object analysis files
                              even if annotated_cells_file already exists; default: FALSE
        --data_files          full paths to each file to be included in analysis
        --manifest            YAML file containing one or more parameter; NOTE: arguments on command 
                              line override manifest arguments!!!         
        --meta_files          comma-delimited list of meta files 
        --meta_data_file      RDA file with pre-compiled/flattened meta data, required IF meta_dir &
                              meta_files are NULL
        --number_threads      number of threads to use for parallel processes
```

