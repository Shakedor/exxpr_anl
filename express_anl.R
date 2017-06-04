#!/usr/bin/env Rscript

source("express_anl_aux.R")


#get arguments
c <- commandArgs(trailingOnly = TRUE)
print(c)
analysis_main(c)


