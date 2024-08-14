#!/usr/bin/env Rscript

if(!require(renv)){
    install.packages("renv")
    library(renv)
}

source("renv/activate.R")

renv::activate()

# Install all packages expect inborutils, which is downloaded 
# in download_zenodo.R
renv::install(exclude = "inborutils")
