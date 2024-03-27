#!/bin/bash/env bash

#   Create a data folder
mkdir -p data

#   Download data from zenodo
if Rscript src/download_zenodo.R ; then
    echo "Data was downloaded successfully!" 
else
    echo "Fail!"
fi