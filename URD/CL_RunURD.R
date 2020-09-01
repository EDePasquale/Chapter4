#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#pull in all Monocle function files
source('R/RunURD.R', local = TRUE)


urd_object=suppressWarnings(RunURD(args[1], #expression
                             args[2], #metadata
                             args[3], #wd
                             args[4], #filename
                             args[5], #rootcluster
                             args[6], #endcluster
                             args[7], #logs.value
                             as.numeric(args[8]), #knn
                             ))


save.image(paste0(args[4],"_results.RData"))