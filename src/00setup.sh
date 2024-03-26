#!/bin/bash

#   Create a data folder
mkdir data

#   Download data from zenodo
Rscript src/00download_zenodo.R

#   Create output folder
mkdir output

#   Run R scripts
Rscript src/01data_analysis_invitro.R
echo "Finished with the in vitro analysis!"

Rscript src/02data_analysis_inplanta.R
echo "Finished with the in planta analysis!"

Rscript src/03mixmodel.R
echo "Finished with the mixture model analysis!"