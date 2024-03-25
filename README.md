# _In vitro_ biofilm formation of antibiotic-resistant _E. coli_ is not correlated with colonisation of lamb’s lettuce leaves

This repository contains the scripts used to analyse the data published in: xxxxxx 

Raw data is storaged in [Zenodo](add link here):

File | Data type | Content
:---: | :---: | :---:
`biofilm.csv` | CSV | In vitro screening of 174 _E. coli_ strains
`cfu.csv` | CSV | CFU data _E. coli_ strains in _V. locusta_
`strain.csv` | CSV | _E. coli_ strain IDs
`strain_abr.csv` | CSV | _E. coli_ antibiotic resistance
`bacimg.tar.gz` | Images | Processed live/dead stains of _E. coli_ strains in the _V. locusta_ phyllosphere

## Scripts
Script | Information
:---: | :---:
`theme_rs.R` | Plot template
`01data_analysis_invitro.R` | Data cleaning, analysis and figures related to _in vitro_ screening
`02data_analysis_inplanta.R` | Data cleaning, analysis and figures related to _in planta_ experiments
`03mixmodel.R` | Code for mixture model analysis
`macro_imgproc.ijm` | FIJI/ImageJ Macro script for image processing
