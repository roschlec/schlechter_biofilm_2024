#!/bin/bash

#   Create a data folder
mkdir -p data

#   Download data from zenodo
Rscript src/00download_zenodo.R
echo "Data was downloaded successfully!" 

#   Create output folder
mkdir -p output

#   Run R scripts
Rscript src/01data_analysis_invitro.R
echo "Finished with the in vitro analysis!"

Rscript src/02data_analysis_inplanta.R
echo "Finished with the in planta analysis!"

Rscript src/03mixmodel.R
echo "Finished with the mixture model analysis!"