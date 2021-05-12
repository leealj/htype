################################################################
################################################################

# Hierarchical clustering of high granularity clustering results
# 
# March 2021:
# - Improve speed of scan_dendrogram() by assigning cells to cluster # instead of
#   iterating through each step (which results in re-calculating the clusters several times).
# - Enabled subsetting with cell ids which will also subset the appropriate branches of 
#   the phylo tree using treeio::tree_subset().
# - Implemented a streamlined multi-gene/multi-gene set analysis with sig_cells() and 
#   scan_dendrogram(). Also updated feature_phylo() to provide a list of plots
#   to accomodate for the multi-gene set analyses that are produced. 
#
# April 2021:
# - Added parameter in recluster_all() to choose between euclidean or correlation 
#   distance when performing BuildClusterTree.
# - Modified BuildClusterTree from Seurat to be able to use correlation distance
#   to calculate distance between centroids.
# - Separated some functions to util.R and vis.R
#
# May 2021:
# - Modified recluster_all and BuildClusterTree to take distance method and linkage 
#   method parameters; default = correlation distance with average linkage
# - Added option to choose best combination of distance and linkage methods based on 
#   cophenetic correlation.
################################################################
################################################################

# Function to add all clustering information across all steps to metadata
# 1. determines all possible clusterings at increments of k + 1 with cutree.
# 2. generates matrix with all cluster assignments based on original cluster labels
# 3. for each step (increments of k + 1), add cluster assignments to meta data. 
# 4. returns Seurat object with phylo object in slot = "BuildClusterTree" and 
#    with all clustering assignments stored in meta data. 
recluster_all <- function(object, 
                          ident = "seurat_clusters", 
                          dims = NULL, 
                          dist.method = "correlation", 
                          linkage = "average") {
  time.start <- Sys.time()
  object.phylo <- Tool(object, slot = "BuildClusterTree")
  
  
  orig_clusters <- object[[ident]][,1]
  if (ident != "seurat_clusters") {
    orig_clusters <- as.factor(orig_clusters-1)
  } # temporary fix on cluster number indexing for SC3 results
  
  object <- AddMetaData(object, metadata = orig_clusters, col.name = "ident")
  
  if (length(VariableFeatures(object)) == 0) object <- FindVariableFeatures(object)
  if (is.null(object.phylo)) {
    object <- BuildClusterTree(object, dims = dims, 
                               features = VariableFeatures(object), 
                               dist.method = dist.method, linkage = linkage)
    object.phylo <- Tool(object, slot = "BuildClusterTree")
  }
  
  total_n_clusters <- max(as.numeric(levels(orig_clusters)))+1
  new_clusterings <- data.frame(row.names = levels(orig_clusters))

  print("slicing tree")
  for (n_clusters in 1:total_n_clusters) {
    new_clusterings[,n_clusters] <- cutree_nodes(object.phylo, k = n_clusters)
  }

  colnames(new_clusterings) <- 1:total_n_clusters
  new_clusterings <- new_clusterings-1  # adjust cluster numbering to start from 0

  print("setting meta data")
  for (step_n in 1:dim(new_clusterings)[2]) {
    # print(paste("setting meta data for step", step_n))
    new_clusters <- new_clusterings[,step_n]
    names(new_clusters) <- levels(orig_clusters)
    new_clusters <- as.factor(new_clusters)
    levels(new_clusters) <- sort(unique(new_clusters))

    for (cluster in levels(orig_clusters)) {
      cells.use <- colnames(object)[orig_clusters == cluster]
      object <- SetIdent(object, cells = cells.use, value = new_clusters[cluster])
    }
    object <- AddMetaData(object, metadata = object@active.ident,
                          col.name = paste0("step_", step_n))
  }

  all_clusters <- sort(unique(concat_steps(object)$cluster))
  internal_nodes_labels <- all_clusters[!all_clusters %in% object.phylo$tip.label]
  object@tools$BuildClusterTree$node.label <- internal_nodes_labels

  time.end <- Sys.time()
  print(time.end-time.start)
  return(object)
}
################################################################
################################################################

# Automatically determines all cells with significant expression of a marker gene
# sig_threshold is only set if test = "wilcoxon".
# marker_genes parameter should be a character vector of genes or a list of gene sets. 
# Returns a character vector of all cell barcodes determined to have significant 
# expression of designated marker gene.
sig_cells <- function(object, scan_output=NULL, test = "wilcoxon",
                      marker_genes, sig_threshold=1e-3) {
  bool_mat <- cluster_bool_mat(object)
  sig_cells_bool <- data.frame(row.names = colnames(object))
  sig_clusters <- list()
  
  # Check if scan_output has been provided
  if (is.null(scan_output)) {
    scan_out <- scan_dendrogram(object = object, test = test, marker_genes = marker_genes)
    marker_genes <- scan_out[5:length(scan_out)]
  } else {
    scan_out <- scan_output
    marker_genes <- scan_out[5:length(scan_out)]
  }
  
  if (is.null(sig_threshold) & (unique(scan_out$test) == "wilcoxon")) {
    sig_threshold <- quantile(scan_out$score)[2]
  }
  
  # Check if recluster_all() has already been run or not
  if (length(grep(names(object@meta.data), pattern = "step_")) == 0) {
    object <- recluster_all(object)
  }
  
  for (set in names(marker_genes)) {
    print(paste("Identifying clusters with", set, "expression"))
    scan_copy <- scan_out
    positive_clusters <- c()
    
    iteration <- 1
    while (min(scan_copy$scores[,set]) < sig_threshold) {
      print(paste("iteration", iteration))
      significant_cluster <- most_sig_cluster(scan_copy, set)
      positive_clusters <- append(positive_clusters, significant_cluster, 
                                  after = length(positive_clusters))
      cells_to_remove <- rownames(bool_mat)[bool_mat[,significant_cluster]]
      clusters_to_remove <- colnames(bool_mat)[apply(bool_mat[cells_to_remove,], 
                                                     MARGIN = 2, FUN = any)]
      scan_copy$scores <- scan_copy$scores[!(rownames(scan_copy$scores) %in% clusters_to_remove),]
      iteration <- iteration + 1
    }
    
    sig_cluster_mat <- bool_mat[,colnames(bool_mat) %in% positive_clusters]
    if (length(positive_clusters) > 1) {
      significant_cells <- apply(sig_cluster_mat, MARGIN = 1, FUN = any)
    } else {significant_cells <- sig_cluster_mat}
    significant_cells <- rownames(bool_mat)[significant_cells]
    sig_cells_bool[,set] <- colnames(object) %in% significant_cells
    
    sig_clusters[[set]] <- positive_clusters
  }
  
  # Generate output
  output <- list()
  output[["cells"]] <- sig_cells_bool
  output[["clusters"]] <- sig_clusters
  output[["scan_output"]] <- scan_out
  output[["phylo"]] <- Tool(object, slot = "BuildClusterTree")
  output[["cell_cluster_matrix"]] <- bool_mat
  output[["project.name"]] <- object@project.name
  
  return(output)
}

# Calculate p-value of a particular marker gene across all steps of dendrogram
# recluster_all() needs to be run on the Seurat object first before running.
# marker_genes parameter should be a character vector of genes or a list of gene sets. 
# test = "wilcoxon" or "ssgsea"
# Called scan_dendrogram() is to be called within sig_cells()
# Returns a data.frame 
scan_dendrogram <- function(object, marker_genes, 
                            test = "wilcoxon", background = NULL) {
  time.start <- Sys.time()
  expr_mat <- as.matrix(GetAssayData(object, slot = "data"))
  
  # Identify root node clusters
  bool_mat <- cluster_bool_mat(object)
  root_clusters <- names(which(apply(bool_mat, MARGIN = 2, FUN = all)))
  exclude_roots <- root_clusters[root_clusters != tail(root_clusters, n = 1)]
  
  # Identify all cluster names
  cluster_cells <- concat_steps(object) 
  all_clusters <- sort(unique(cluster_cells$cluster))
  all_clusters <- all_clusters[!(all_clusters %in% exclude_roots)]
  
  if (class(marker_genes) != "list") {
    marker_genes <- list(marker_genes)
    names(marker_genes) <- paste(unlist(marker_genes), collapse = "_")
  } else if (is.null(names(marker_genes))) {
    names(marker_genes) <- paste0("set_", 1:length(marker_genes))
  }
  
  all_scores <- data.frame(row.names = all_clusters)
  all_scores[,"cluster"] <- all_clusters
  all_avg_expr <- data.frame(row.names = all_clusters)
  all_avg_expr[,"cluster"] <- all_clusters
  
  for (gene_set in names(marker_genes)) {
    print(paste("calculating gene set", gene_set))
    test_set <- marker_genes[[gene_set]]
    
    score <- c()
    avg_expr <- c()
    
    # Use only genes present in data set
    keep_genes <- test_set[test_set %in% rownames(object)]
    
    if (test == "wilcoxon") {
      i <- 0
      for (cluster in all_clusters) {
        # Print every 5 clusters
        if ((i/5)%%1==0) {
          print(paste("calculating cluster", cluster)) 
        }
        cells <- unique(cluster_cells$cell_id[cluster_cells$cluster == cluster])
        coi_expr <- expr_mat[keep_genes, (colnames(object) %in% cells)]
        rest_expr <- expr_mat[keep_genes, !(colnames(object) %in% cells)]
        
        if (cluster %in% root_clusters) { 
          # set p-value to 1 on step 1 because there is only one cluster.
          score <- append(score, values = 1, after = length(score))
          avg_expr <- append(avg_expr, values = mean(coi_expr), after = length(avg_expr))
        }
        if (!(cluster %in% root_clusters)) {
          if (length(keep_genes) > 1) {
            coi_expr <- as.numeric(apply(coi_expr, MARGIN = 2, FUN = mean))
            rest_expr <- as.numeric(apply(rest_expr, MARGIN = 2, FUN = mean))
          }
          
          cluster_pval <- wilcox.test(coi_expr, rest_expr, alternative = "greater")[["p.value"]]
          score <- append(score, values = cluster_pval, after = length(score))
          avg_expr <- append(avg_expr, values = mean(coi_expr), after = length(avg_expr))
        }
        i <- i + 1
      }
      score[score < 1e-250] <- 1e-250
      all_scores[,gene_set] <- score
      all_avg_expr[,gene_set] <- avg_expr
    }
  }
  output <- list()
  output[["scores"]] <- all_scores
  output[["avg_expr"]] <- all_avg_expr
  output[["test"]] <- test
  output[["root_clusters"]] <- root_clusters
  for (set in names(marker_genes)) {
    output[[set]] <- marker_genes[[set]]
  }
  time.end <- Sys.time()
  print(time.end-time.start)
  return(output)
}
#####
