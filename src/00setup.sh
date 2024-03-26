#!/bin/bash/env bash

#   Install R dependencies
Rscript src/dependencies.R

#   Create a data folder
mkdir -p data

#   Download data from zenodo
if Rscript src/download_zenodo.R ; then
    echo "Data was downloaded successfully!" 
else
    echo "Fail!"
fi

#   Create output folder
mkdir -p output

#   Run R scripts
if Rscript src/01data_analysis_invitro.R ; then   
    echo "Finished with the in vitro analysis!"  
else
    echo "Fail!"
fi

if Rscript src/02data_analysis_inplanta.R ; then   
    echo "Finished with the in planta analysis!" 
else
    echo "Fail!"
fi

if Rscript src/03mixmodel.R ; then   
    echo "Finished with the mixture model analysis!" 
else
    echo "Fail!"
fi