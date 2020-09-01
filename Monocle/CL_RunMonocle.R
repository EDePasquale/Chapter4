#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#pull in all Monocle function files
source('R/monocle_DR_functions.R', local = TRUE)
source('R/KT_RunMonocle2_v2.R', local = TRUE)


cds=suppressWarnings(RunMonocle(args[1], #set directory
                              args[2], #expr file
                              args[3], #final groups file
                              args[4], #expressionFamily
                              args[5], #fullModelFormulaStr for differentialGeneTest
                              args[6], #initial_method for RGE technique
                              args[7], #reduction_method (RGE technique)
                              args[8],#filename
                              args[9], #logs.value
                              args[10]))  #root

save.image(paste0(args[8], ".RData"))

