#!/usr/bin/env Rscript

if(!require(renv)){
    install.packages("renv")
    library(renv)
}

source("renv/activate.R")