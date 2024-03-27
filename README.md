#  _In vitro_ biofilm formation of antibiotic-resistant _E. coli_ is not correlated with colonisation of lamb’s lettuce leaves

This repository contains the scripts used to analyse the data published in: xxxxxx 

These scripts were developed using R version 4.3.3.

Raw data is storaged in [Zenodo](https://zenodo.org/records/10871054):

File | Data type | Content
:---: | :---: | :---:
`biofilm.csv` | CSV | _In vitro_ screening of 174 _E. coli_ strains
`cfu.csv` | CSV | CFU data _E. coli_ strains in _V. locusta_
`strain.csv` | CSV | _E. coli_ strain IDs
`strain_abr.csv` | CSV | _E. coli_ antibiotic resistance
`bacimg.tar.gz` | Images | Processed live/dead stains of _E. coli_ strains in the _V. locusta_ phyllosphere

##  Scripts
Script | Information
:---: | :---:
`download_data.sh` | Bash script to download the raw data
`analysis.sh` | Bash script to run the _in vitro_ and _in planta_ analyses
`01data_analysis_invitro.R` | Data cleaning, analysis and figures related to _in vitro_ screening
`02data_analysis_inplanta.R` | Data cleaning, analysis and figures related to _in planta_ experiments
`03mixmodel.R` | Code for mixture model analysis
`dependencies.R` | Activate the R environment with required packages
`download_zenodo.R` | Code for installing and downloading the raw data from Zenodo
`theme_rs.R` | Plot template 
`macro_imgproc.ijm` | FIJI/ImageJ Macro script for image processing

##  Setting up 
First, you need to clone this repository.
```
git clone https://github.com/roschlec/schlechter_biofilm_2024.git
```

To download dependencies and packages used for data analysis, run the following command from the project directory.
```
Rscript src/dependencies.R
```

##  Download datasets
To download the datasets associated to the manuscript, run the following command from the project directory. This script will create a data directory and download the raw data and images files from Zenodo.
```
bash src/download_data.sh
```

##  Run analysis
To run the full _in vitro_ and _in planta_ analyses, you must run the `analysis.sh` script from the project home directory.
```
bash src/analysis.sh
```
This will generate the draft plots used in the manuscript.


