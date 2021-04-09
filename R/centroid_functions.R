# Calculate cluster centroids

# Takes a Seurat object with cluster assignments at each step in meta data
# Returns list of matrices. Each matrix representing all the centroids for a 
# given step. 
centroids <- function(object, reduction = NULL, dims = NULL) {
  
  if (is.null(reduction)) {
    space <- t(as.matrix(GetAssayData(object, slot = "data")[VariableFeatures(object),]))
  } else if (reduction == "pca") {
    space <- object@reductions[[reduction]]@cell.embeddings[,dims]
  } else {
    space <- object@reductions[[reduction]]@cell.embeddings
  }
  
  total_steps <- length(names(object@meta.data)[grep(names(object@meta.data), pattern = "step_")])
  all_clusters <- data.frame(row.names = colnames(object))
  output <- list()
  for (step in 1:total_steps) {
    all_clusters[,step] <- object@meta.data[paste0("step_", step)][,1]
    mat <- matrix(NA, nrow = length(levels(all_clusters[,step])), ncol = dim(space)[2])
    n <- 1
    for (i in levels(all_clusters[,step])) {
      # Space of a specific cluster
      cluster_features <- space[all_clusters[,step] == as.numeric(i),]
      mat[n,] <- apply(cluster_features, MARGIN = 2, FUN = mean)
      n <- n + 1
    }
    colnames(mat) <- colnames(space)
    rownames(mat) <- as.numeric(levels(all_clusters[,step]))
    output[[paste0("step_", step)]] <- mat
  } 
  return(output)
}

# Calculates euclidean distance between all points (cells) in an matrix
# and all centroids in a matrix. 
# points1 = matrix with single cells
# points2 = matrix with centroids
# cells/centroids should be rows and 
# principle components/expression should be columns
# Returns a matrix where rows are single cells, columns are centroids
myEuclid <- function(points1, points2) {
  distanceMatrix <- matrix(NA, nrow=dim(points1)[1], ncol=dim(points2)[1])
  colnames(distanceMatrix) <- rownames(points2)
  rownames(distanceMatrix) <- rownames(points1)
  for (i in 1:nrow(points2)) {
    distanceMatrix[,i] <- sqrt(rowSums(t(t(points1)-points2[i,])^2))
  }
  return(distanceMatrix)
}

CorDist <- function(points1, points2, method = "pearson") {
  output <- 1-(cor(t(points1), t(points2), method = method))
  return(output)
}


# Calculates all distances from every cell to each cluster centroid at every step. 
# Returns a list, with each element representing each step as a matrix of distances
# with single cells as the rows and cluster centroids as columns.
sc_dist <- function(object, distance.method = "euclidean", corr.method = NULL, 
                    reduction = NULL, dims = NULL) {
  output <- list()
  
  print("calculating centroids")
  centers <- centroids(object, reduction = reduction, dims = dims)
  
  if (is.null(reduction)) {
    space <- t(as.matrix(GetAssayData(object, slot = "data")[VariableFeatures(object),]))
  } else if (reduction == "pca") {
    space <- object@reductions[[reduction]]@cell.embeddings[,dims]
  } else {
    space <- object@reductions[[reduction]]@cell.embeddings
  }
  
  if (distance.method == "euclidean") {
    print("calculating Euclidean distances")
    for (step in names(centers)) output[[step]] <- myEuclid(space, centers[[step]])
  }
  if (distance.method == "correlation") {
    print("calculating correlation distances")
    for (step in names(centers)) output[[step]] <- CorDist(space, centers[[step]],
                                                           method = corr.method)
  }
  
  return(output)
}


# output
foo <- function(object, distance.method = "euclidean", corr.method = NULL, 
                reduction = NULL, dims = NULL) {
  distances <- sc_dist(object, reduction = reduction, dims = dims, 
                       distance.method = distance.method, corr.method = corr.method)
  output <- list()
  for (step in names(distances)) {
    output[[step]] <- list()
    
    nearestc_vector <- c()
    cutoff_vector <- c()
    bool_mat <- data.frame(row.names = colnames(object))
    
    for (cell in rownames(distances[[step]])) {
      # 2 shortest distances for each cell
      if (step == "step_1") {
        nearest_clusters <- colnames(distances[[step]])
      } else {
        nearest_clusters <- names(sort(distances[[step]][cell,]))[1:2] 
      } 
      nearestc_vector <- append(nearestc_vector, values = nearest_clusters, 
                                after = length(nearestc_vector))
    }
    for (cluster in colnames(distances[[step]])) {
      seurat_cluster_ids <- object@meta.data[step] == as.numeric(cluster)
      # grab the 2nd quantile of dist values for each cluster
      cutoff <- quantile(distances[[step]][seurat_cluster_ids,cluster])[3] 
      bool_mat[,cluster] <- distances[[step]][,cluster] < cutoff
    }
    
    if (step == "step_1") {
      top_clusters_matrix <- t(matrix(nearestc_vector, nrow = 1,
                                      ncol = length(colnames(object))))
    } else {
      top_clusters_matrix <- t(matrix(nearestc_vector, nrow = 2, 
                                      ncol = length(colnames(object))))
    }

    rownames(top_clusters_matrix) <- colnames(object)
    output[[step]][["dist"]] <- distances[[step]]
    output[[step]][["top_clusters"]] <- top_clusters_matrix
    output[[step]][["bool_mat"]] <- bool_mat
  }
  return(output)
}
