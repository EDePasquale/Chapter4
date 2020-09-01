# ERICA MONOCLE 2 FUNCTION


# Load Monocle and other required libraries
library(monocle)
require(reshape2)
library(Seurat)
library(Biobase)
library(Seurat)
library(plotly)
library(dplyr)
library(cowplot)
library(htmlwidgets)

RunMonocle <- function(wd,path_to_rna,path_to_groups,expressionFamily, fullModelFormulaStr,initial_method,reduction_method, filename, log.values=FALSE, ROOT="auto", BRANCH=1, HMCLUS=20){
  
  fullModelFormulaStr = paste0("~",fullModelFormulaStr)
  
  # Set file location
  setwd(wd)
  exp<- read.table(file=path_to_rna, sep="\t", stringsAsFactors=F, header=T, check.names=F)
  exp <- exp[order(rowSums(exp[,2:length(colnames(exp))]), decreasing=T),]
  exp <- exp[which(!duplicated(exp[,1])),]
  rownames(exp) <- exp[,1]
  exp <- exp[,c(-1)]
  #exp <- (2^exp)-1
  exp = ceiling(exp)
  idx_r = which(rowSums(exp)==0)
  if (length(idx_r)==0){exp=exp}else{exp = exp[-idx_r,]}
  #exp = exp[-(which(rowSums(exp)==0)),]
  idx_c = which(colSums(exp)<1)
  if (length(idx_c)==0){exp=exp}else{exp = exp[,-idx_c]}
  
  
  
  # Read in the groups
  groups <- read.table(file=path_to_groups, sep="\t", stringsAsFactors=F, header=T, 
                       check.names=F)
  
  #colnames(groups) <- c("Cluster_Name","Groups")
  #rownames(groups) <- groups$UID
  if (length(idx_c)==0){groups=groups}else{groups = groups[-idx_c,]}
  #groups  = groups[-idx,]
  rownames(groups) <- colnames(exp)
  #groups = cbind(row.names(groups),groups)
  colnames(groups) <- c("UID","Groups","Cluster_Name")
  
  cell_intersect <- intersect(colnames(exp), groups$UID)
  groups <- groups[cell_intersect,]
  exp <- exp[,cell_intersect]
  
  # Data is currently log normalized TPM, so we'll exponentiate it
  if (log.values=="TRUE"){log.values=TRUE}
  
  if(log.values=="FALSE"){log.values=FALSE}
  
  if(log.values==TRUE){exp <- (2^exp)-1}
  
  # Create cell annotation df (pd) and gene annotation df (fd) for CellDataSet
  pd <- new("AnnotatedDataFrame", data=groups)
  gene_groups <- data.frame(gene_short_name=rownames(exp),
                            row.names=rownames(exp),
                            stringsAsFactors=F)
  fd <- new("AnnotatedDataFrame", data=gene_groups)
  
  # Make CDS Object and estimate size factors and dispersions if neg binomial is used 
  if(expressionFamily == "tobit"){
    myCDS <- newCellDataSet(as(data.matrix(exp), "sparseMatrix"),phenoData=pd,featureData=fd,expressionFamily=tobit())
  }
  
  if(expressionFamily == "negbinomial.size"){
    myCDS <- newCellDataSet(as(data.matrix(exp), "sparseMatrix"),phenoData=pd,featureData=fd,expressionFamily=negbinomial.size())
    myCDS <- estimateSizeFactors(myCDS)
    myCDS <- estimateDispersions(myCDS)
  }
  
  
  ##############
  # CLEAN DATA #
  ##############
  
  # Estimate size factors and dispersions
  # myCDS <- estimateSizeFactors(myCDS)
  # myCDS <- estimateDispersions(myCDS)
  
  #Make a directory for saving results
  dir.create(filename)
  setwd(paste0(wd, "/", filename))
  
  # Filter out low quality cells
  myCDS <- detectGenes(myCDS, min_expr = 0.1) #genes expressed per cell
  expressed_genes <- row.names(subset(fData(myCDS), num_cells_expressed >= 4))
  valid_cells <- row.names(subset(pData(myCDS))) #filtering on genes expressed in 4 or more cells
  myCDS <- myCDS[,valid_cells]
  pData(myCDS)$Total_mRNAs <- Matrix::colSums(exprs(myCDS)) #mRNA distribution across all cells
  myCDS <- myCDS[,pData(myCDS)$Total_mRNAs < 1e6]
  upper_bound <- 10^(mean(log10(pData(myCDS)$Total_mRNAs)) + 2*sd(log10(pData(myCDS)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(myCDS)$Total_mRNAs)) - 2*sd(log10(pData(myCDS)$Total_mRNAs)))
  
  #Set colors for seperated plot by group
  group_vec <- unique(pData(myCDS)$Groups)
  NUMCLUS=length(unique(pData(myCDS)$Groups))
  group_cols <- colorRampPalette(RColorBrewer::brewer.pal(name="Set1", n = 8))(NUMCLUS)
  group_cols[6] <- "#6A3D9A"
  names(group_cols) <- group_vec
  pdf(paste0("mRNA-Counts-plot-", filename, ".pdf"), width = 10, height = 8) #create and save plot
  plot<-qplot(Total_mRNAs, data = pData(myCDS),color = Groups, geom = "density") + geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound) + scale_color_manual(values = group_cols, name = "Groups")
  print(plot)  
  dev.off()
  myCDS <- myCDS[,pData(myCDS)$Total_mRNAs > lower_bound & pData(myCDS)$Total_mRNAs < upper_bound] #remove cells with very low and very high mRNA
  myCDS <- detectGenes(myCDS, min_expr = 0.1)
  L <- log(exprs(myCDS[expressed_genes,])+1) #Log-transform each value in the expression matrix.
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L)))) #Standardize each gene, so that they are all on the same scale, then melt the data with plyr so we can plot it easily
  qplot(value, geom = "density", data = melted_dens_df)+ stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(FPKM)") + ylab("Density") #Plot the distribution of the standardized gene expression values.
  
  
  ###################################
  # CLUSTERING WITHOUT MARKER GENES #
  ###################################
  
  # Filter by high expression level and unusual variability
  disp_table <- dispersionTable(myCDS)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  myCDS <- setOrderingFilter(myCDS, unsup_clustering_genes$gene_id)
  pdf(paste0("plot-ordering-genes-", filename, ".pdf"), width = 10, height = 8) #create and save plot
  plot<-plot_ordering_genes(myCDS)
  print(plot)
  dev.off()
  pdf(paste0("variance-", filename, ".pdf"), width = 10, height = 8) #create and save plot
  plot<-plot_pc_variance_explained(myCDS, norm_method = "none",return_all = F) 
  print(plot)
  dev.off()
  
  # Try clustering cells
  myCDS <- reduceDimension(myCDS, max_components = 2, num_dim = 6, reduction_method = 'tSNE', verbose = T,check_duplicates=F)
  myCDS <- clusterCells(myCDS, num_clusters = 26) #set the number of clusters at or higher than what you think you need
  pdf(paste0("Cluster-cells-", filename, ".pdf"), width = 10, height = 8) #create and save plot
  plot<-plot_cell_clusters(myCDS, 1, 2, color='as.factor(Groups)') + scale_color_manual(values = group_cols, name = "Groups")
  print(plot)
  dev.off()
  
  
  #########################################
  # CONSTRUCTING SINGLE CELL TRAJECTORIES #
  #########################################
  
  # Choose genes that define a cell's progress
  diff_test_res <- differentialGeneTest(myCDS[expressed_genes,], fullModelFormulaStr = fullModelFormulaStr) #group-based DE genes
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) #genes used for pseudotemporal ordering
  myCDS <- setOrderingFilter(myCDS, ordering_genes)
  pdf(paste0("Select-cells-best-for-clustering-", filename, ".pdf"), width = 10, height = 8) #create and save plot
  plot<-plot_ordering_genes(myCDS)
  print(plot)
  dev.off()
  
  # Dimensionality reduction
  #myCDS <- reduceDimension(myCDS, method=reduction_method, initial_method = match.fun(initial_method), max_components = 2)
  if (length(valid_cells) < 200) {myCDS <- reduceDimension(myCDS, method=reduction_method, max_components = 2,auto_param_selection = F)}else{
    myCDS <- reduceDimension(myCDS, method=reduction_method, max_components = 2)
  }
  
  
  # Order cells along the trajectory
  myCDS <- orderCells(myCDS)
  pdf(paste0("order-cells-trajectory-Groups-", filename, ".pdf"), width = 10, height = 8) #create and save plot
  plot=plot_cell_trajectory(myCDS, color_by = 'as.factor(Groups)') + scale_color_manual(values = group_cols, name = "Groups")
  print(plot)
  dev.off()
  # ggsave(paste0("order-cells-trajectory-Groups1-", filename, ".png"), width=15, height=15) #create and save plot (better aspect ratios)
  # htmlwidgets::saveWidget(as_widget(ggplotly()), 'order-cells-trajectory-Groups1.html')
  # plot<-plot_cell_trajectory(myCDS, color_by = 'as.factor(Groups)') + scale_color_manual(values = group_cols, name = "Groups")
  # print(plot)
  # dev.off()
  
  # Save progress
  saveRDS(myCDS, file = paste0("myCDS_", filename, ".rds"))
  
  # Selecting start state
  ##  NOTE: Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning", so we often have to call 
  ##        orderCells again using the root_state argument to specify the beginning. First, we plot the trajectory, this time 
  ##        coloring the cells by "State"
  pdf(paste0("order-cells-trajectory-state-",filename,".pdf"), width = 10, height = 8)
  plot<-plot_cell_trajectory(myCDS, color_by = "State")
  print(plot)
  dev.off()
  
  # Function for identifying the State which contains most of the cells from time zero
  GM_state <- function(myCDS){
    if (length(unique(pData(myCDS)$State)) > 1){
      T0_counts <- table(pData(myCDS)$State, pData(myCDS)$Groups)[,4]
      return(as.numeric(names(T0_counts)[which
                                         (T0_counts == max(T0_counts))]))
    } else {
      return (1)
    }
  }
  if(ROOT=="auto"){
    root_state=GM_state(myCDS)
  }else{
    root_state=ROOT
  }
  myCDS <- orderCells(myCDS, root_state)
  pdf(paste0("Order-Pseudotime-", filename, ".pdf"), width = 10, height = 8) #create and save plot
  plot<-plot_cell_trajectory(myCDS, color_by = "Pseudotime")
  print(plot)
  dev.off()
  
  ##  NOTE: If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
  ##        Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located
  pdf(paste0("Order-Pseudotime-states-", filename, ".pdf"), width = 10, height = 8) #create and save plot
  plot<-plot_cell_trajectory(myCDS, color_by = "State") + facet_wrap(~State, nrow = 1)
  print(plot)
  dev.off()
  
  # Black and white state v group graph
  state_cluster_stat <- table(pData(myCDS)[, c('State', 'Groups')])
  
  state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
  state_cluster_stat_ordered <- t(state_cluster_stat)
  
  options(repr.plot.width=3, repr.plot.height=3)
  pdf(paste0("State-Group-Heatmap-", filename, ".pdf"), width = 8, height = 10) #create and save plot
  plot<-pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(25))
  print(plot)
  dev.off()
  
  ##########################
  # CREATE COMPLICATE TREE #
  ##########################
  
  # Set colors for full plot
  type_vec <- unique(pData(myCDS)$Groups)
  type_cols <- colorRampPalette(RColorBrewer::brewer.pal(name="Set1", n = 8))(NUMCLUS)
  type_cols[6] <- "#6A3D9A"
  names(type_cols) <- type_vec
  
  # Set plot options
  options(repr.plot.width=6, repr.plot.height=4)
  
  #Groups
  pdf(paste0("Complicate-Tree-", filename,".pdf"), width = 10, height = 8)
  plot=plot_complex_cell_trajectory(myCDS[, ], color_by = "Groups", show_branch_points = T, cell_size = 0.5, cell_link_size = 0.5, root_states=root_state) + 
    scale_color_manual(values = group_cols, name = "Groups")
  print(plot)
  dev.off()
  
  #State
  state_vec <- unique(pData(myCDS)$State)
  state_cols <- RColorBrewer::brewer.pal(9, name = 'Set1')
  state_cols[6] <- "#6A3D9A"
  names(state_cols) <- state_vec
  pdf(paste0("Complicate-Tree-States-", filename, ".pdf"), width = 10, height = 8)
  plot<-plot_complex_cell_trajectory(myCDS[, ], color_by = "State", show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3, root_states=root_state)
  print(plot)
  dev.off()
  
  #Set colors for seperated plot by group
  group_vec <- unique(pData(myCDS)$Groups)
  group_cols <- colorRampPalette(RColorBrewer::brewer.pal(name="Set1", n = 8))(NUMCLUS)
  group_cols[6] <- "#6A3D9A"
  names(group_cols) <- group_vec
  
  pdf(paste0("Complicate-Tree-Groups-", filename, ".pdf"), width = 15, height = 8)
  plot<-plot_complex_cell_trajectory(myCDS[, ], color_by = "Groups", show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3, root_states=root_state) + facet_wrap(~Groups, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = group_cols, name = "Groups")
  print(plot)
  dev.off()
  
  ####################
  # SAVE ENVIRONMENT #
  ####################
  
  save.image(paste0("myCDS_", filename, "_env.Rdata"))
  
  
  #################
  # BEAM ANALYSIS #
  #################
  
  BEAM_res <- BEAM(myCDS, branch_point = BRANCH, cores = 1) # branch 2 is the first, then 1 in this case
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  # pdf(paste0("BEAM_",BRANCH,"_", filename, ".pdf"), width = 8, height = 15)
  # plot<-plot_genes_branched_heatmap(myCDS[row.names(subset(BEAM_res,
  #                                                          qval < 1e-20)),],
  #                                   branch_point = BRANCH,
  #                                   #num_clusters = HMCLUS,
  #                                   cores = 1,
  #                                   use_gene_short_name = T,
  #                                   show_rownames = T)
  # print(plot)
  # dev.off()
  # # 
  # BEAM_hm=plot_genes_branched_heatmap(myCDS[row.names(subset(BEAM_res,
  #                                                            qval < 1e-20)),],
  #                                     branch_point = BRANCH,
  #                                     #num_clusters = HMCLUS,
  #                                     cores = 1,
  #                                     use_gene_short_name = T,
  #                                     show_rownames = T,
  #                                     return_heatmap = T)
  # # 
  # 
  # #########################
  # # RECREATE BEAM HEATMAP #
  # #########################
  # 
  # library(pheatmap)
  # library(pals)
  # 
  # if(!is.null(MYCOLORS)){
  #   my_colour=MYCOLORS
  # }else{
  #   my_colour = list(
  #     `Cell Type` = BEAM_hm[["annotation_colors"]][["Cell Type"]],
  #     Cluster = c(`1`="#0000FF", `2`="#FF0000" ,`3`="#00FF00", `4`="#000033", `5`="#FF00B6", `6`="#005300", `7`="#FFD300",
  #                 `8`="#009FFF", `9`="#9A4D42", `10`="#00FFBE", `11`="#783FC1", `12`="#1F9698", `13`="#FFACFD", `14`="#B1CC71",
  #                 `15`="#F1085C", `16`="#FE8F42", `17`="#DD00FF", `18`="#02AD24", `19`="#720055", `20`="#766C95")
  #   )
  # }
  # 
  # pdf(paste0("BEAM_manual-",BRANCH,"_", filename, ".pdf"), width = 8, height = 15)
  # plot<-ph_res <- pheatmap(BEAM_hm$heatmap_matrix, #ph$tree_row$order
  #                          useRaster = T,
  #                          cluster_cols=FALSE, 
  #                          cluster_rows=TRUE, 
  #                          show_rownames=FALSE, 
  #                          show_colnames=F, 
  #                          #scale="row",
  #                          clustering_distance_rows=BEAM_hm$row_dist, #row_dist
  #                          clustering_method = "ward.D2", #ward.D2
  #                          cutree_rows=HMCLUS,
  #                          # cutree_cols = 2,
  #                          annotation_row=BEAM_hm$annotation_row,
  #                          annotation_col=BEAM_hm$annotation_col,
  #                          #annotation_colors=BEAM_hm_first_state$annotation_colors,
  #                          #annotation_colors=colorRampPalette(RColorBrewer::brewer.pal(name="Set1", n = 8))(20),
  #                          annotation_colors=my_colour,
  #                          gaps_col = BEAM_hm$col_gap_ind,
  #                          treeheight_row = 20, 
  #                          breaks=seq(-3.1,3.1, length.out = length(BEAM_hm$hmcols)),
  #                          fontsize = 6,
  #                          color=BEAM_hm$hmcols, 
  #                          border_color = NA,
  #                          silent=FALSE)
  # print(plot)
  # dev.off()
  
  ####################################
  # GO BIOLOGICAL PROCESS (prefered) #
  ####################################
  #
  # x.clust <- cutree(BEAM_hm[["ph_res"]][["tree_row"]], k = HMCLUS)
  # x.uniq=unique(x.clust)
  # #BiocManager::install("GOFunction")
  # #BiocManager::install("org.Hs.eg.db")
  # library("GOFunction")
  # library("org.Hs.eg.db")
  # library("mygene")
  # refGenes=queryMany(as.character(Genes[,1]), species="human", fields="entrezgene", scopes=c("symbol","alias","name"), return.as="DataFrame") #convert to entrezid
  # GO_results=NA
  # for(rrow in x.uniq){
  #   print(rrow)
  #   sink("/dev/null") #hides mygene output
  #   geneEquiv=queryMany(names(which(x.clust==rrow)), species="human", fields="entrezgene", scopes=c("symbol","alias","name"), return.as="DataFrame") #convert to entrezid
  #   sink()
  #   entrezIDs=geneEquiv@listData$entrezgene #grab entrezids
  #   sigTerm <- GOFunction(entrezIDs, refGenes$`_id`, organism="org.Hs.eg.db",
  #                         ontology="BP", fdrmethod="BY", fdrth=0.25, ppth=0.25, pcth=0.25,
  #                         poth=0.25, peth=0.25, bmpSize=2000, filename=paste0("sigTerm_", filename, "_", rrow))
  #   print(sigTerm$name[1:min(10,length(sigTerm$name))])
  #   GO_results=rbind(GO_results, cbind(paste0("Cluster ", rrow, " BP"),"Pvalue"), cbind(sigTerm$name, sigTerm$pvalue), "")
  # }
  # write.table(GO_results[-1,], paste0("GO_BP_summary-",BRANCH,"_", filename, ".txt"), sep="\t", row.names=F, col.names = F)
  # 
  save.image(paste0("myCDS_",BRANCH, "_",filename, "_env.Rdata"))
  
}
