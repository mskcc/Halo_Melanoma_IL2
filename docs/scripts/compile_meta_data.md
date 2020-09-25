## compile_meta_data.R

Usage:  
```
Rscript compile_meta_data.R 
            
    [REQUIRED (may be defined on command line OR in manifest file)] 
      --meta_dir        path to meta files in XLSX format, required if 
                        meta_data_file is NULL
      --meta_data_file  RDA or XLSX file to which newly flattened data 
                        should be saved

      [OPTIONAL]
        --meta_files    comma-delimited character string with each element 
                        containing path to one meta file; use this when there 
                        are multiples of one or more meta file in a directory
        --manifest      YAML file containing one or more required parameter; 
                        NOTE: arguments on command line override manifest 
                        arguments!!!     
```
## Output:

RDA file containing a list with the following elements:

**Patients:** tibble created from [StudyName]_Patients.xlsx
**Samples:** tibble created from [StudyName]_Samples.xlsx
- any samples for which ALL FOVs are marked for exclusion in [StudyName]_FOVs.xlsx are excluded
- Sample_ID is generated in the form [Patient_ID]_[Sample_number]
**FOVs:** tibble created from [StudyName]_FOVs.xlsx
- any FOVs marked for exclusion are excluded
- FOV_ID is generated in the form [Patient_ID]_[Sample_number]_[FOV_number]
**flat:** tibble that is a full join of Patients, Samples and FOVs by CellDive_ID 
**IDs:** tibble that is a map of all FOV identifiers 
