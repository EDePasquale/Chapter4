RunURD <- function(expression,metadata,wd, filename, rootcluster, endcluster,log.values=TRUE, knn=100){
  
  library(URD)
  library(ggplot2)
  library(Matrix)
  library(tools)
  
  #Set directory
  setwd(wd)
  
  #Read the files
  if (file_ext(expression) == "csv"){
    count.object <- as.matrix(read.table(expression,sep=",",header = FALSE))
  }
  
  
  
  #write.table[count.object, file="VanGalen_CHHi_strip_26_Sampled_0.25.txt",quote=F,sep="\t",row.names=F,colnames=F]
  
  if (file_ext(expression) == "txt"){
    count.object <- read.table(expression,sep="\t",header = TRUE,check.names = F,stringsAsFactors = FALSE)
  }
  
  #count.object<- read.table(file=expression, sep="\t", stringsAsFactors=F, header=T, check.names=F)
  count.object <- count.object[order(rowSums(count.object[,2:length(colnames(count.object))]), decreasing=T),]
  count.object <- count.object[which(!duplicated(count.object[,1])),]
  rownames(count.object) <- count.object[,1]
  count.object <- count.object[,c(-1)]
  count.object = as.matrix(count.object)
  #count.object = as(count.object,"dgCMatrix")

  
  if (log.values=="TRUE"){log.values = TRUE}
  if (log.values == "FALSE"){log.values=FALSE}
  
  if (log.values==TRUE){count.object = (2^count.object)}else{count.object=count.object}
  
  #metadata_try = "downsample_urd_metadata.txt"
  #meta.object_try <- read.table(metadata_try,sep="\t", stringsAsFactors=F, header=T, check.names=F)
  
  meta.object <- read.table(metadata,sep="\t", stringsAsFactors=F, header=T, check.names=F)
  #genes <- read.table(genes,header=T,sep=",")
  endcluster = as.matrix(read.table(endcluster,sep="\t",header = FALSE))
  endcluster = endcluster[,1]
  
  if (nrow(meta.object) != ncol(count.object)){
    print("ncol of counts and nrow of metadata don't match.")
  }
  
  #Ensure the row names of meta data and colnames of counts matrix match (set to barcodes)
  row.names(meta.object) = meta.object[,1]
  #colnames(count.object) = rownames(meta.object)
  
  #Set the rownames of the count matrix with genes 
  #rownames(count.object) = genes[,1] 
  
  print('Creating URD object...')
  object <- createURD(count.data = count.object, meta = meta.object, min.cells=3, min.counts=3)
  
  
  
  setwd(wd)
  dir.create(filename)
  setwd(paste0(wd, "/", filename))
  
  
  # Copy stage from @meta to @group.ids 
  object@group.ids$group <- as.character(object@meta[["group"]])
  
  
  # Get variable genes for each group
  print('Calculating variable genes...')
  uniqueGroups <- unique(object@group.ids$group)
  
  # new_list_clusters = c()
  # for (i in 1:length(uniqueGroups)){
  #   if (length(cellsInCluster(object,'group',uniqueGroups[i] > 5))){
  #     new_list_clusters = c(new_list_clusters, i)
  #   }
  # }
  
  # myfunction  <- function(n,object,uniqueGroups)
  #   if(length(cellsInCluster(object,clustering = 'group',cluster=uniqueGroups[n])) == 1){
  #     uniqueGroups_filter = uniqueGroups[-n]
  #     return(uniqueGroups_filter)
  #   }
  
  uniqueGroups_original = uniqueGroups
  
  for (i in 1:length(uniqueGroups)){
    if(length(cellsInCluster(object,clustering = 'group',cluster=uniqueGroups[i])) == 1){
      uniqueGroups = uniqueGroups[-i]}
  }
  
  
  
  var.by.group <- lapply(1:length(uniqueGroups),function(n){
    pdf(file=paste0(filename,'var_genes_',uniqueGroups[n],'.pdf'))
    x <- findVariableGenes(object,cells.fit=cellsInCluster(object, 'group', uniqueGroups[n]), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=0.005, mean.max=100, main.use=paste0('Group ', uniqueGroups[n]), do.plot=T)
    dev.off()
    return(x)
  })
  
  var.genes <- sort(unique(unlist(var.by.group)))
  object@var.genes <- var.genes
  
  
  print('Calculating PCA plot...')
  object <- calcPCA(object, mp.factor = 2,do.print=TRUE)
  
  pca_score_matrix_sig = as.matrix(object@pca.scores[, which(object@pca.sig)])
  ncols_pca_sig_matrix = ncol(pca_score_matrix_sig)

  remove_cells_list = c()

  for (i in 1:ncols_pca_sig_matrix){
    remove_cells_list = c(remove_cells_list, which(duplicated(pca_score_matrix_sig[,i])))
  }

  remove_cells_list = unique(remove_cells_list)

  if (length(remove_cells_list) == 0){
    object = object
    print("No duplicates of PCA scores. 0 cells removed.")
    }else{
    object_original = object
    all_barcodes = object@count.data@Dimnames[[2]]
    perc_cell_removed = (length(remove_cells_list)/length(all_barcodes))*100
    print(paste("Removing ",perc_cell_removed," percent of cells because of duplicate PCA scores"))
    object = urdSubset(object,cells.keep = all_barcodes[-remove_cells_list])
  }
  
  
  
  print('Calculating tsne plot...')
  set.seed(19)
  object <- calcTsne(object = object)
  
  print('Calculating diffusion map...')
  object <- calcDM(object, knn = knn, sigma='local')
  print('Finished Diffusion maps...')
  
  pdf(paste0(filename,"_DiffusionMaps_Grid.pdf"), width = 10, height = 8) #create and save plot
  DM_array = plotDimArray(object, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map (Sigma Local, 100 NNs): Groups", label="group", plot.title="", legend=F)
  print(DM_array)
  dev.off()
  
  pdf(paste0(filename,"_DiffusionMaps_Transitions.pdf"), width = 10, height = 8) #create and save plot
  transition_DM = plotDim(object, "group", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")
  print(transition_DM)
  dev.off()
  
  print('Setting root cells...')
  # Here we use all cells from the first stage as the root
  root.cells <- cellsInCluster(object, "group", cluster = rootcluster)
  
  print('Running floodPseudotime simulation...')
  # Then we run 'flood' simulations
  object.floods <- floodPseudotime(object, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)
  indexing <- which(!(is.na(object.floods[1,])))
  object.floods <- object.floods[,indexing]
  
  print('Processing flood simulation into pseudotime...')
  # The we process the simulations into a pseudotime
  object <- floodPseudotimeProcess(object, object.floods, floods.name="pseudotime")
  
  pdf(paste0(filename,"_PseudotimeStability.pdf"), width = 10, height = 8) #create and save plot
  PT_Stability = pseudotimePlotStabilityOverall(object)
  print(PT_Stability)
  dev.off()
  
  pdf(paste0(filename,"_Pseudotime_tSNE.pdf"), width = 10, height = 8) #create and save plot
  PT_tSNE = plotDim(object, "pseudotime")
  print(PT_tSNE)
  dev.off()
  
  pdf(paste0(filename,"_Pseudotime_Dists.pdf"), width = 10, height = 8) #create and save plot
  PT_Dists = plotDists(object, "pseudotime", "group", plot.title="Pseudotime by Groups")
  print(PT_Dists)
  dev.off()
  
  print('Setting endclusters...')
  endcluster_original = endcluster 
  
  for (i in 1:length(endcluster)){
    if(length(cellsInCluster(object,clustering = 'group',cluster=endcluster[i])) == 1){
      endcluster = endcluster[-i]}
  }
  
  
  # Create a subsetted object of just those cells from the final stage
  object.terminal <- urdSubset(object, cells.keep=cellsInCluster(object, "group", endcluster))
  
  # Use the variable genes that were calculated only on the final group of stages (which
  # contain the last stage).
  print('Variable genes for terminal cluster...')
  
  terminal_var_genes = c()
  for (i in 1:length(endcluster)){
    if (length(which(uniqueGroups==endcluster[i])) == 0){endcluster=endcluster[-i]} else {
      terminal_var_genes = c(terminal_var_genes,var.by.group[[which(uniqueGroups==endcluster[i])]])
    }
  }
  
  object.terminal@var.genes <- terminal_var_genes
  
  print('Calculating PCA Tsne of terminal clusters...')
  # Calculate PCA and tSNE
  object.terminal <- calcPCA(object.terminal, mp.factor = 1.5)
  
  pca_score_matrix_sig_terminal = as.matrix(object.terminal@pca.scores[, which(object.terminal@pca.sig)])
  ncols_pca_sig_matrix_terminal = ncol(pca_score_matrix_sig_terminal)
  
  remove_cells_list_terminal = c()
  
  for (i in 1:ncols_pca_sig_matrix_terminal){
    remove_cells_list_terminal = c(remove_cells_list_terminal, which(duplicated(pca_score_matrix_sig_terminal[,i])))
  }
  
  remove_cells_list_terminal = unique(remove_cells_list_terminal)
  
  if (length(remove_cells_list_terminal) == 0){
    object.terminal = object.terminal
    print("No duplicates of PCA scores in terminal URD subset. 0 cells removed.")
  }else{
    object_original_terminal = object.terminal
    all_barcodes_terminal = object.terminal@count.data@Dimnames[[2]]
    perc_cell_removed_terminal = (length(remove_cells_list_terminal)/length(all_barcodes_terminal))*100
    print(paste("Removing ",perc_cell_removed_terminal," percent of cells from end state clusters because of duplicate PCA scores"))
    object.terminal = urdSubset(object.terminal,cells.keep = all_barcodes_terminal[-remove_cells_list_terminal])
  }
  
  
  #object@group.ids[rownames(object.terminal@group.ids), "tip.clusters"] <- object.terminal@meta[rownames(object.terminal@group.ids),'cluster_number']
  
  # pdf(paste0(filename,"PCA_elbow_terminal.pdf"), width = 10, height = 8) #create and save plot
  #   PCA_elbow_terminal = pcSDPlot(object.terminal)
  #   print(PCA_elbow_terminal)
  # dev.off()
  
  set.seed(20)
  object.terminal <- calcTsne(object.terminal)
  
  # Calculate graph clustering of these cells
  object.terminal <- graphClustering(object.terminal, num.nn = 50, do.jaccard=T, method="Louvain")
  
  # pdf("Terminal_clusters_louvain.pdf", width = 10, height = 8) #create and save plot
  #   terminal_clusters_louvain = plotDim(object.terminal, "Louvain-50", plot.title = "Louvain (50 NN) graph clustering", point.size=3)
  #   print(terminal_clusters_louvain)
  # dev.off()
  
  object.terminal@group.ids$group_number = character(length = length(object.terminal@group.ids[["group"]]))
  
  for (i in 1:length(object.terminal@group.ids[["group"]]))
  {
    for (g in 1:length(endcluster)){
      if (object.terminal@group.ids$group[i] == endcluster[g]){
        object.terminal@group.ids$group_number[i] = which(uniqueGroups==endcluster[g])
      }
    }
  }
  #   if (object.terminal@group.ids$group[i] == "ERP4"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="ERP4")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "immNeu"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="immNeu")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "MP"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="MP")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "cMoP"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="cMoP")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "Mast cell"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="Mast cell")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "Basophil"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="Basophil")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "Eosinophils"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="Eosinophils")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "CD19+"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="CD19+")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "MKP"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="MKP")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "NK"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="NK")
  #   }
  #   
  #   if (object.terminal@group.ids$group[i] == "DC"){
  #     object.terminal@group.ids$group_number[i] = which(groups=="DC")
  #   }
  #   
  # }
  
  # Copy cluster identities from object.terminal object to a new clustering ("tip.clusters") in the full object object.
  object@group.ids[rownames(object.terminal@group.ids), "tip.clusters"] <- object.terminal@group.ids$group_number
  
  # Determine the parameters of the logistic used to bias the transition probabilities. The procedure
  # is relatively robust to this parameter, but the cell numbers may need to be modified for larger
  # or smaller data sets.
  object.ptlogistic <- pseudotimeDetermineLogistic(object, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)
  
  
  # Bias the transition matrix acording to pseudotime
  object.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(object, "pseudotime", logistic.params=object.ptlogistic))
  
  # Noticed that the biased transition matrix doesn't include all cells, so I have to segment
  # the original URD object to just the cells included
  object.seg <- urdSubset(object,cells.keep=colnames(object.biased.tm))
  
  #urd_object = list(full_urd = object,seg_urd = object.seg,tree = object.tree,terminal_urd = object.terminal,ptlogistic = object.ptlogistic)
  
  save.image(paste0(filename,"_results_before_Walks.RData"))
  
  # Simulate the biased random walks from each tip
  object.walks <- simulateRandomWalksFromTips(object.seg, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = object.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = T)
  
  # Process the biased random walks into visitation frequencies
  object.seg <- processRandomWalksFromTips(object.seg, walks=object.walks, verbose = F)
  
  object.seg@tree$pseudotime.breakpoint.details <- 0
  
  pdf(paste0(filename,"Tip_clusters_tSNE.pdf"), width = 10, height = 8) #create and save plot
  tip_tsne = plotDim(object, "tip.clusters", plot.title="Cells in each tip")
  print(tip_tsne)
  dev.off()
  
  # pdf("Transition_visitfreq_1.pdf", width = 10, height = 8) #create and save plot
  #   transition_visitfreq1 = plotDim(object, "visitfreq.log.1", plot.title="Visitation frequency from tip 1 (log10)", transitions.plot=10000)
  #   print(transition_visitfreq1)
  # dev.off()
  # 
  # pdf("Transition_visitfreq_2.pdf", width = 10, height = 8) #create and save plot
  #   transition_visitfreq2 = plotDim(object, "visitfreq.log.2", plot.title="Visitation frequency from tip 2 (log10)", transitions.plot=10000)
  #   print(transition_visitfreq2)
  # dev.off()
  # 
  # pdf("Transition_visitfreq_3.pdf", width = 10, height = 8) #create and save plot
  #   transition_visitfreq3 = plotDim(object, "visitfreq.log.3", plot.title="Visitation frequency from tip 3 (log10)", transitions.plot=10000)
  #   print(transition_visitfreq3)
  # dev.off()
  
  # Load the cells used for each tip into the URD object
  object.tree <- loadTipCells(object.seg, "tip.clusters")
  
  tips_use = c()
  for (i in 1:length(endcluster)){
    ec_idx = which(uniqueGroups==endcluster[i])
    tips_use = c(tips_use,ec_idx)
  }
  
  
  # Build the tree
  object.tree <- buildTree(object.tree, pseudotime = "pseudotime", tips.use=tips_use, divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)
  
  # Name the segments based on our previous determination of the identity of tips 1 and 2.
  object.tree <- nameSegments(object.tree, segments=tips_use, segment.names = endcluster, short.names = endcluster)
  
  pdf(paste0(filename,"_URDTree_Groups.pdf"), width = 10, height = 8) #create and save plot
  URD_tree = plotTree(object.tree, "group", title="Developmental Stage")
  print(URD_tree)
  dev.off()
  
  pdf(paste0(filename,"_URDTree_Pseudotime.pdf"), width = 10, height = 8) #create and save plot
  URD_tree_pt = plotTree(object.tree, "pseudotime", title="Developmental Stage")
  print(URD_tree_pt)
  dev.off()
  
  pdf(paste0(filename,"_Tree_segment_tSNE.pdf"), width = 10, height = 8) #create and save plot
  URD_tree_segment_tsne = plotDim(object.tree, "segment", plot.title="URD tree segment")
  print(URD_tree_segment_tsne)
  dev.off()
  
  pdf(paste0(filename,"_URDTree_segment.pdf"), width = 10, height = 8) #create and save plot
  URD_tree_segment = plotTree(object.tree, "segment", title="URD Tree Segment")
  print(URD_tree_segment)
  dev.off()
  
  visitation = data.frame(cell= rownames(object.tree@diff.data),seg=object.tree@diff.data$segment,stringsAsFactors = F,row.names = rownames(object.tree@diff.data))
  
  visitation$visit= log10(apply(visitation,1,function(cr) object.tree@diff.data[as.character(cr["cell"]),paste0("visitfreq.raw.",as.character(cr["seg"]))])+1)
  
  robustly.visited.cells = visitation[visitation$visit>=0.5,"cell"]
  
  final.tips = segTerminal(object.tree)
  
  object.tree = treeForceDirectedLayout(object.tree,num.nn=120,method="fr",cells.to.do = robustly.visited.cells,tips=final.tips,cut.unconnected.segments = 2,min.final.neighbors = 4,verbose=F)
  
  numclus = length(unique(uniqueGroups_original))
  color_vec <- unique(uniqueGroups_original)
  cluster_cols <- colorRampPalette(RColorBrewer::brewer.pal(name="Set1", n = 8))(numclus)
  
  pdf(paste0(filename,"_URDTreeForce_2D.pdf"), width = 10, height = 8)
  plot = plotTreeForce2D(object.tree, "group", title = "Groups",point.alpha = 0.8) + scale_color_manual(values = cluster_cols)+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(colour = guide_legend(override.aes = list(size=7)))
  print(plot)
  dev.off()
  
  urd_object = list(full_urd = object,seg_urd = object.seg,tree = object.tree,terminal_urd = object.terminal,ptlogistic = object.ptlogistic,biased_tm = object.biased.tm,walks = object.walks)

  return(urd_object)
}
