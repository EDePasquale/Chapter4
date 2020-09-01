#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#pull in all Monocle function files
source('R/TI_tool_functions.R', local = TRUE)
source('R/t_col.R', local = TRUE)


final_object=suppressWarnings(RunSlingshot(args[1], #labels
                                           args[2], #dr_coordinates
                                           args[3], #wd
                                           args[4], #filename
                                           args[5], #start_cluster
                                           if(args[5]=="NULL"){NULL}
                                           args[6], #end_cluster
                                           if(args[6]=="NULL"){NULL}
                                           args[7], #labels_name
                                           if(args[7]=="NULL"){NULL}
                                           args[8], #cols vector
                                           if(args[8]=="NULL"){NULL}
                                           as.numeric(args[9]), #pch
                                           as.numeric(args[10]), #cex
                                           as.numeric(args[11]), #percent of transparency in DR plot
))

save.image(paste0(args[4],"_results.RData"))