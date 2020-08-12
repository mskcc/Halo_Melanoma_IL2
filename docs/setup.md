## Prerequisites & setup instructions

1) Install [Python 3](https://www.python.org/downloads/) if necessary.

3) Install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

4) Install [R verion >= 3.6](https://cloud.r-project.org/) if necessary.

4) Make sure the correct R version is on your PATH.

5) Check for the following R libraries: 
     - cowplot
     - funr
     - ggplot2
     - gridExtra
     - logger
     - parallel
     - randtoolbox
     - rgeos
     - R.utils
     - SearchTrees
     - sp
     - tidyverse
     - xlsx
     - yaml 

6) Modify [run.sh](../pipeline/run.sh) and/or [serial_pipeline.sh](../pipeline/serial_pipeline.sh) to set correct PYTHONPATH and activate virtual environment if necessary. Alternatively, these can be set outside of the script and removed or commented out. 
