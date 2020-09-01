library(slingshot)
library(RColorBrewer)

RunSlingshot <- function(labels,plot,wd,filename,startcluster=NULL,endcluster=NULL, labels.name=NULL,cols=NULL,pch=16,cex=0.5,perc=70){
  
  #Set directory
  setwd(wd)
  
  #Load and read labels and plot files 
  
  labels_data=read.table(labels, sep="\t", header=T)
  plot_data=read.table(plot, sep="\t", header=T)
  
  data = as.matrix(plot_data)
  tokeep <- which(sapply(labels_data,is.numeric))
  clusterLabels = labels_data[, tokeep]
  
  if (is.character(endcluster)){
    endcluster = read.table(endcluster,sep="\t",header = FALSE)
    endcluster = as.character(as.matrix(endcluster))
  }
  
  #Run Slingshot
  
  if (is.null(startcluster) & is.null(endcluster)){
    results = slingshot(data, clusterLabels, allow.breaks=FALSE)
  }
  
  if (is.null(startcluster) & is.character(endcluster)){
    results = slingshot(data, clusterLabels, allow.breaks=FALSE,end.clus = endcluster)
  }
  
  if (is.null(endcluster) & is.character(startcluster)){
    results = slingshot(data, clusterLabels, allow.breaks=FALSE,start.clus=startcluster)
  }
  
  if (is.character(endcluster) & is.character(startcluster)){
    results = slingshot(data, clusterLabels, allow.breaks=FALSE,start.clus=startcluster,end.clus=endcluster)
  }
 
  #if labels.name file provided, make vector of label names
  if(!is.null(labels.name)){
    labels_name = read.table(labels.name, sep="\t", header=T)
    labels_name = as.matrix(labels_name[,1])
  }
  

  
  #Plot the data and results and save them
  nb.cols <- length(unique(clusterLabels))
  clusternumber_cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(nb.cols)
  clusternumber_vec <- unique(clusterLabels)
  names(clusternumber_cols) <- clusternumber_vec
  
  #generate colors if not provided
  if (is.character(cols)){
    cols = read.table(cols, sep="\t",header=TRUE)
    cols = as.matrix(cols[,1])
  }
  
  if (is.null(cols)){
    brewercols = vector(length = length(clusterLabels))
    
    for(i in 1:length(clusternumber_vec)){
      idx = which(clusterLabels==clusternumber_vec[i])
      brewercols[idx] = clusternumber_cols[i]
    }
    
    cols = brewercols
  }
  

  
  #generate transparent version of colors  
  transparent_cols = vector()
  for(i in 1:length(clusterLabels)){
    transparent_cols[i] <- t_col(cols[i],perc=perc,name="transp")
  }
  

  #Generate colors for lineages 
  n_lineage=length(results@lineages)
  lineage_cols=rainbow(n_lineage,s=1,v=1,start = 0, end = max(1, n_lineage - 1)/n_lineage,alpha=1)
  
  #Generate colors of curves
  n_curve=length(results@curves)
  curve_cols=rainbow(n_curve,s=1,v=1,start = 0, end = max(1, n_curve - 1)/n_curve,alpha=1)
  
  #Create new folder to save results
  dir.create(filename)
  setwd(paste0(wd, "/", filename))
  
  
  #DR plot with no trajectories 
  pdf(paste0(filename,"_DR_plot.pdf"), width = 10, height = 8)
    DR_plot= plot(data, col=cols, pch=pch, cex=cex)
    if(is.null(labels.name)){
      legend("topright",legend=as.matrix(unique(clusterLabels)), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    if (!is.null(labels.name)){
      legend("topright",legend=unique(labels_name), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    print(DR_plot)
  dev.off()
  
  #DR plot with lineages (in black)
  
  pdf(paste0(filename,"_DR_plot_lineage.pdf"), width = 10, height = 8)
    DR_plot_lineage= plot(data, col=brewercols, pch=pch, cex=cex)
    if(is.null(labels.name)){
      legend("topright",legend=as.matrix(unique(clusterLabels)), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    if (!is.null(labels.name)){
      legend("topright",legend=unique(labels_name), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    lines(SlingshotDataSet(results), lwd=2,type='lineage',col='black',show.constraints=TRUE)
    print(DR_plot_lineage)
  dev.off()
  
  #DR plot with curves (in black)
  
  pdf(paste0(filename,"_DR_plot_curves.pdf"), width = 10, height = 8)
    DR_plot_curves= plot(data, col=cols, pch=pch, cex=cex)
    if(is.null(labels.name)){
      legend("topright",legend=as.matrix(unique(clusterLabels)), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    if (!is.null(labels.name)){
      legend("topright",legend=unique(labels_name), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    lines(SlingshotDataSet(results), lwd=2,type='curve',col='black',show.constraints=TRUE)
    print(DR_plot_curves)
  dev.off()
  
  #DR plot with lineage (in color)
  
  pdf(paste0(filename,"_DR_plot_lineage_color.pdf"), width = 10, height = 8)
    DR_plot_lineage_color= plot(data, col=transparent_cols, pch=pch, cex=cex)
    
    for(i in 1:n_lineage) {
      lines(results, linInd = i, type='l',cex=1.5,col=lineage_cols[i],lwd = 2,show.constraints=TRUE)
    }
    
    if(is.null(labels.name)){
      legend("topright",legend=as.matrix(unique(clusterLabels)), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    if (!is.null(labels.name)){
      legend("topright",legend=unique(labels_name), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    print(DR_plot_lineage_color)
  dev.off()
  
  #DR plot with curves (in color)
  
  pdf(paste0(filename,"_DR_plot_curve_color.pdf"), width = 10, height = 8)
    DR_plot_curve_color= plot(data, col=transparent_cols, pch=pch, cex=cex)
    for(i in 1:n_curve) {
      lines(results, linInd = i, type='c',cex=1.5,col=curve_cols[i],lwd = 2,show.constraints=TRUE)
    }
    
    if(is.null(labels.name)){
      legend("topright",legend=as.matrix(unique(clusterLabels)), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    if (!is.null(labels.name)){
      legend("topright",legend=unique(labels_name), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    print(DR_plot_curve_color)
  dev.off()
  
  #DR plot with lineage and curve 
  
  pdf(paste0(filename,"_DR_plot_lineage_curve.pdf"), width = 10, height = 8)
    DR_plot_lineage_curve= plot(data, col=transparent_cols, pch=pch, cex=cex)
    for(i in 1:n_lineage) {
      lines(results, linInd = i, type='l',cex=1.5,col=lineage_cols[i],lwd = 2,show.constraints=TRUE)
    }
    
    lines(SlingshotDataSet(results), lwd=1,type='curve',col='black',show.constraints=TRUE)
    
    if(is.null(labels.name)){
      legend("topright",legend=as.matrix(unique(clusterLabels)), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    
    if (!is.null(labels.name)){
      legend("topright",legend=unique(labels_name), col=unique(cols), pch=pch, cex=cex,pt.cex = 0.8)
    }
    print(DR_plot_lineage_curve)
  dev.off()
  
  #Text file output with lineage info 
  
  for(i in 1:length(results@lineages)){
    cat(paste("Lineage",as.character(i)),file='lineage_info.txt',append = TRUE)
    write.table(results@lineages[[i]],file='lineage_info.txt',append=TRUE, row.names=FALSE,col.names =FALSE,quote = FALSE)
  }
  
  
  save.image(paste0(filename, "_env.Rdata"))
  saveRDS(results, file = paste0(filename, "_slingshot_object.rds"))
  
  final_object = list(labels_data = labels_data,plot_data = plot_data,data_matrix = data,slingshot_results=results)
  
  return(final_object) 
}