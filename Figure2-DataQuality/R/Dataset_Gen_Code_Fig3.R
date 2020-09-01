###########################
#                         #
# DATASET GENERATION CODE #
#         Figure 3        #
#                         #
###########################


####################
# Downsample Cells #
#    Figure 3A     #
####################

# EXP1 = filename of tab-delimited expression file log-scaled containing no meta data
# EXP2 = filename of tab-delimited expression file counts containing no meta data
# VALUE = proportion of cells that will remain (1 -> 0)
# GROUPING = labels in the format of SPRING
Downsample_Cells <- function(EXP1, EXP2, VALUE=1, GROUPING){
  
  
  #read in expression file
  FILENAME1=sub(".txt", "", EXP1)
  FILENAME2=sub(".txt", "", EXP2)
  EXP1=read.table(EXP1, sep="\t", header=T, row.names=1, stringsAsFactors = F)
  EXP2=read.table(EXP2, sep="\t", header=T, row.names=1, stringsAsFactors = F)
  
  #read in grouping file
  GROUPING=read.csv(GROUPING, sep=",", header=F, row.names=1)
  colnames(GROUPING) <- colnames(EXP1) #same for both EXP1 and EXP2
  
  #sample cells without replacement
  SAMPLED_CELLS=sample(colnames(EXP1), size=ncol(EXP1)*VALUE, replace=FALSE)
  NEW_EXP1 <- EXP1[,SAMPLED_CELLS]
  NEW_EXP2 <- EXP2[,SAMPLED_CELLS]
  
  #reduce groupings file
  GROUPING2 <- GROUPING[,SAMPLED_CELLS]
  
  #write new files
  write.table(NEW_EXP1, paste0(FILENAME1, "_Sampled_", VALUE, ".txt"), sep="\t", row.names = TRUE, col.names=NA, quote=FALSE)
  write.table(NEW_EXP1, paste0(FILENAME1, "_Sampled_", VALUE, ".csv"), sep=",", row.names = FALSE, col.names=FALSE, quote=FALSE)
  write.table(NEW_EXP2, paste0(FILENAME2, "_Sampled_", VALUE, ".txt"), sep="\t", row.names = TRUE, col.names=NA, quote=FALSE)
  write.table(NEW_EXP2, paste0(FILENAME2, "_Sampled_", VALUE, ".csv"), sep=",", row.names = FALSE, col.names=FALSE, quote=FALSE)
  write.table(GROUPING2, paste0(FILENAME1, "_Grouping_", VALUE, ".csv"), sep=",", col.names = F, row.names = TRUE, quote=FALSE)
}

####################
# Downsample Genes #
#    Figure 3B     #
####################

# EXP1 = filename of tab-delimited expression file log-scaled containing no meta data
# EXP2 = filename of tab-delimited expression file counts containing no meta data
# VALUE = proportion of genes that will be removed (0 -> 1)
Downsample_Genes <- function(EXP1, EXP2, VALUE=0){
  
  #read in expression file
  FILENAME1=sub(".txt", "", EXP1)
  FILENAME2=sub(".txt", "", EXP2)
  NEW_EXP1=read.table(EXP1, sep="\t", header=T, row.names=1, stringsAsFactors = F)
  NEW_EXP2=read.table(EXP2, sep="\t", header=T, row.names=1, stringsAsFactors = F)
  
  #exponentiate CPTT input
  NEW_EXP1[]=lapply(NEW_EXP1, exp)
  NEW_EXP1=NEW_EXP1-1
  
  #create dropout by randomly sampling from an exponential distribution then subtracting from each cellxgene
  num_elements=nrow(NEW_EXP1)*ncol(NEW_EXP1) #same for both EXP1 and EXP2
  SCALE_FACTOR=NEW_EXP2/NEW_EXP1
  DOWNSAMPLE_VAL=rexp(num_elements, rate=unlist(list(-log(VALUE)/NEW_EXP1)))
  NEW_EXP1[]=NEW_EXP1-DOWNSAMPLE_VAL
  NEW_EXP2[]=NEW_EXP2-DOWNSAMPLE_VAL*SCALE_FACTOR
  #NEW_EXP[]=NEW_EXP-rexp(num_elements, rate=unlist(list(-log(VALUE)/NEW_EXP)))

  #fix
  NEW_EXP1[NEW_EXP1<0]<-0
  NEW_EXP1[is.na(NEW_EXP1)]<-0
  NEW_EXP2[NEW_EXP2<0]<-0
  NEW_EXP2[is.na(NEW_EXP2)]<-0
  
  #rescale data
  NEW_EXP1=NEW_EXP1+1
  NEW_EXP1[]=lapply(NEW_EXP1, log)

  #ceiling the counts data
  NEW_EXP2=ceiling(NEW_EXP2)
  
  #write new file
  write.table(NEW_EXP1, paste0(FILENAME1, "_Dropout_", VALUE, ".txt"), sep="\t", row.names = FALSE)
  write.table(NEW_EXP1, paste0(FILENAME1, "_Dropout_", VALUE, ".csv"), sep=",", row.names = FALSE, col.names=FALSE)
  write.table(NEW_EXP2, paste0(FILENAME2, "_Dropout_", VALUE, ".txt"), sep="\t", row.names = FALSE)
  write.table(NEW_EXP2, paste0(FILENAME2, "_Dropout_", VALUE, ".csv"), sep=",", row.names = FALSE, col.names=FALSE)
  
}
