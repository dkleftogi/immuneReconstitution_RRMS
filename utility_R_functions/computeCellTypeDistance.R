#@Dimitrios Kleftogiannis NeuroSys-Med, UiB 2023
#computes the cosine similarity between the centroids of cell populations

#' @param healthy.AB.data is the data matrix
#' @param col.names are the channels/markers that will be used for distance estimation
#' @param filtering_T The threshold used to filter edges in the graph. If this value is \code{< 1} the edges are filtered based on the
#'   cosine similarity value. If this is an integer \code{>= 1} they are filtered based on rank for each node (i.e. for each node only
#'   the edges with rank less than the threshold are retained)


computeCellTypeDistance <- function(healthy.AB.data.tab,col.names,filtering_T){
  
  att <- as.matrix(healthy.AB.data.tab[,col.names])
  row.names(att) <- healthy.AB.data.tab$cellType
  #compute the cosine similarity between the populations found in the dataset
  dd <- cosine_similarity_matrix(att)
  diag(dd) <- 0 #here we need to decide if we need to keep the self loops
  dd[is.na(dd)] <- 0 #This can happen if one of the attractors has all 0's for the markers of interest
  
  if(filtering_T >= 1){
    filter_matrix_by_rank(dd, filtering_T)
  }else{
    filter_matrix(dd, filtering_T)
  }
  return(dd)
}
#load the filtering functions from vite package too
filter_matrix <- function(m, threshold) {
  invisible(.Call('_vite_filter_matrix', PACKAGE = 'vite', m, threshold))
}

filter_matrix_by_rank <- function(m, threshold) {
  invisible(.Call('_vite_filter_matrix_by_rank', PACKAGE = 'vite', m, threshold))
}