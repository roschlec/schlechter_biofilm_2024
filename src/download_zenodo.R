#!/usr/bin/env Rscript

# Install the package inborutis and load into library
if(!require(devtools)){install.packages("devtools")}
if(!require(inborutils)){devtools::install_github("inbo/inborutils")}

library(inborutils)

doi <- "10.5281/zenodo.10871054"
local_path <- "data"
inborutils::download_zenodo(doi, local_path, quiet = TRUE)
list.files(local_path)