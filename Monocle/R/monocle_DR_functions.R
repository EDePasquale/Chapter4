#following are functions for different dimension reduction methods:
PCA <- function(data, ...) {
  res <- prcomp(t(data), center = F, scale = F)
  res$x
}
ICA <- function(data, max_components = 3, ...) {
  data <- data[, !(duplicated(t(data)))]
  monocle:::ica_helper(t(data), n.comp = max_components,  use_irlba = F)$S #dimension need to be less than variable
}

LLE2 <- function(data, m = 2, k = 5){
  lle_reduced_data <- lle(t(data), m = m, k = k)
  res <- t(lle_reduced_data$Y)
  colnames(res) <- colnames(data)
  return(t(res))
}

LLE <- function (data, max_components = 3, num_neigh = NULL, reg = 2, ss = FALSE,
                 id = TRUE, v = 0.9, iLLE = FALSE) {
  if(is.null(num_neigh)) {
    # num_neigh_list <- calc_k(t(data), m = max_components, kmin = 1,
    #                          kmax = 20, plotres = TRUE, parallel = T, cpus = detectCores(), iLLE = iLLE)
    # num_neigh <- num_neigh_list$k[which(num_neigh_list$rho ==
    #                                       min(num_neigh_list$rho))]
    
    #use slicer's approach: 
    k = SLICER::select_k(t(data), kmin=5)
    message('k is ', k)
  }
  lle_reduced_data <- lle(t(data), m = max_components, k = k,
                          reg = reg, ss = ss, id = id, v = v, iLLE = iLLE)
  
  res <- t(lle_reduced_data$Y)
  colnames(res) <- colnames(data)
  return(t(res))
}
ISOMAP <- function(data, max_components = 3){
  tmp <- isomap(dist(t(data)), ndim = max_components, k = 3, fragmentedOK = T)
  res <- tmp$points
  row.names(res) <- colnames(X)
}

destiny_diffusionMaps <- function(data, max_components = 3, sigma = NULL, k = find_dm_k(nrow(data) - 1L),
                                  n.eigs = min(20L, nrow(data) - 2L), density.norm = TRUE,
                                  ..., distance = c("euclidean"), censor.val = NULL,
                                  censor.range = NULL, missing.range = NULL, vars = NULL, verbose = !is.null(censor.range)){
  data <- t(data)
  tmp <- destiny::DiffusionMap(data, sigma = sigma, k = k,
                               n.eigs = n.eigs, density.norm = density.norm, distance = distance, censor.val = censor.val,
                               censor.range = censor.range, missing.range = missing.range, vars = vars, verbose = verbose)
  res <- tmp@eigenvectors
  row.names(res) <- row.names(data)
  
  return(res[, 1:max_components])
}
