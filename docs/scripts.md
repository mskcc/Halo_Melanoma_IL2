##Individual scripts run by pipeline
To print usage of any of these scripts, run:

```
Rscript [script.R] -h
```

### Validate meta data
#### &nbsp; &nbsp; [../scripts/validate_meta.R](../scripts/validate_meta.R)
Validate all manually created meta data XLSX files (files described [here](study_meta_data.md))

### Configure study
#### &nbsp; &nbsp;  [../scripts/configure_study.R](../scripts/configure_study.R)
Run within study directory, subfolders will be created and input files linked in a standard format unless config has been customized

### Compile all meta data
####  &nbsp; &nbsp; [../scripts/data/compile_meta.R](../scripts/data/compile_meta.R)
Parse and join all meta data. Store in a list each individual table in addition to a single flattened table and one of just a map of all IDs. List is saved in a RDA file.

### Convert boundaries from XML to RDA
####  &nbsp; &nbsp; [../scripts/data/parse_halo_boundaries.R](../scripts/data/parse_halo_boundaries.R)
Convert boundary data from XML to flat tables & store in RDA file.

### Mark cells for exclusion
####  &nbsp; &nbsp; [../scripts/data/mark_exclusions.R](../scripts/data/mark_exclusions.R)
Add column to Halo data indicating whether a marker, cell or FOV should be excluded from all analyses. Any field in this column that is not NA or the empty string will be excluded.

### Annotate cells
####  &nbsp; &nbsp; [../scripts/data/annotate_cells.R](../scripts/data/annotate_cells.R)
Filter out exclusions and consolidate marker-level data from Halo into cell-level data for downstream analysis. Assign classifiers to each cell according to the definitions in meta file [StudyName]_CellTypes.xlsx (description [here](study_meta_data.md))

### Summarize cell counts
####   &nbsp; &nbsp; [../scripts/data/count_cell_types.R](../scripts/data/count_cell_types.R)
Count occurrences of each possible positive marker combination and summarize by counting number of cells that fall into each classifier both in each FOV and in each sample. Output is a single XLSX file with a separate sheet for each summary

### Calculate area of cell regions
#### [../scripts/data/calculate_area.R](../scripts/calculate_area.R)
Calculate total cell area in each FOV and also in each 10micron band around the tumor interface, from -360:360 microns.

### Pre-calculate all metrics for downstream analyses
#### [../scripts/data/calculate_metrics.R](../scripts/calculate_metrics.R)
Precalculate all metrics to be analyzed including fractions, densities, neighborhood averages and neighborhood fractions as described in [StudyName]_Conditions.xlsx. Default is run calculation per FOV, but can also be run per Sample or per interface band. Calculations are saved to RDA files.

### Run statistics and generate report for study questions
#### [../scripts/data/report_statistics.R](../scripts/data/report_statistics.R)
For each question outlined in [StudyName]_Questions.xlsx, run Wilcox test to compare two sets of samples. Generate a single XLSX file containing tabs for both complete and filtered results for each calculation type (fraction, density, etc.)
<br>
