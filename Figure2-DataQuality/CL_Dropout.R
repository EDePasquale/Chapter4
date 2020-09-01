#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#pull in all Monocle function files
source('R/Dataset_Gen_Code_Fig3.R', local = TRUE)

setwd(args[4])
Downsample_Genes(args[1], args[2], as.numeric(args[3]))