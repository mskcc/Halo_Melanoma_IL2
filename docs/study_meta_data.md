# Study meta data

#### General file rules

* all file names should start with prefix [StudyName]
* all files must contain at minimum the columns described below but may contain additional columns as well
* a column must exist for any clinical variables to be analyzed, including variables to be compared and variables on which data will be filtered
* it is strongly recommended that column names NOT contain any spaces or special characters
* column names may NOT start with a number
* repeated values within a column must match EXACTLY including case and spaces
* no values should contain leading or trailing whitespace
* column names should not be used in multiple XLSX files except for variables to used for joining tables (CellDive_ID, Patient_ID)

<br>

#### [StudyName]_Patients.xlsx
All meta data corresponding to individual patients

| Column | Description |
| ----------- | ----------- |
| Patient_ID | integer representing one patient |
Other possible columns include: Patient_response

<br>

#### [StudyName]_Samples.xlsx
All meta data corresponding to individual samples
| Column | Description |
| ----------- | ----------- |
| CellDive_ID | Unique sample identifier contained in names of Halo files (images, CSV and XML boundary files)
| Patient_ID| Patient_ID from [StudyName]_Patients.xlsx corresponding to the patient from which the sample came from |
| Sample_number | integer representing a single sample within a single patient (i.e., this number is unique within a patient but may be repeated in multiple patients) |
Other possible columns include: Lesion_response, Treatment, Treated, etc.

<br>

#### [StudyName]_FOVs.xlsx
All meta data corresponding to individual FOVs

| Column | Description |
| ----------- | ----------- |
| CellDive_ID | Unique sample identifier of sample that FOV came from |
| FOV_number | integer representing a single FOV within a single sample (i.e., this number is unique within a sample but may be repeated across samples) |
| FOV_type | e.g., 'Center', 'Interface', 'Regressed' |
| FOV_used_for_thresholding | an 'X' in this column indicates the FOV was used for manual thresholding of the sample from which it came |
| FOV_exclusion | an 'X' here indicates the FOV should NOT be included in ANY analyses |
| Marker_exclusion | comma-separated string of marker names to be excluded from an FOV |
| Num_manual_exclusions | number of manually drawn regions of cells to be excluded including epidermis, exclusion and glass regions (for QA & debugging marked/filtered exclusions) |
| Num_epidermis_exclusions | number of manually drawn epidermis regions (in case it is decided that epidermis should NOT be excluded)|
| Num_interface_exclusions | number of manually drawn regions of tumor cells |

<br>

#### [StudyName]_Markers.xlsx
names and descriptions of all markers used in study
| Column | Description |
| ----------- | ----------- |
| Marker_name | name of marker exactly as it appears in  |
| Description | what kind of marker is this? e.g., Identity or Functional |
| Cell_type | an 'X' here means the marker is used in one or more cell type definitions; cells that do not contain any positive cell_type markers are classified as 'super negative'|
| Threshold_compartment | |
| In_cell_types | for functional markers, value is a comma-delimited list of cell types in which marker may be positive |

<br>

#### [StudyName]_Files.xlsx
manifest of all halo files per FOV
| Column | Description |
| ----------- | ----------- |
| CellDive_ID | unique sample identifier to which files correspond |
| FOV_number | integer representing a single FOV within a single sample (i.e., this number is unique within a sample but may be repeated across samples) |
| Halo_data_file | full path to Halo CSV file |
| Halo_image_file | full path to Halo image file |
| Halo_boundary_file | full path to XML file containing all exclusion and tumor boundaries for a single FOV |
