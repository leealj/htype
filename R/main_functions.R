################################################################
################################################################

# Hierarchical clustering of high granularity clustering results
# 
# Development version 4:
# - Improve speed of scan_dendrogram() by assigning cells to cluster # instead of
#   iterating through each step (which results in re-calculating the clusters several times).
# - Enabled subsetting with cell ids which will also subset the appropriate branches of 
#   the phylo tree using treeio::tree_subset().
# - Implemented a streamlined multi-gene/multi-gene set analysis with sig_cells() and 
#   scan_dendrogram(). Also updated feature_phylo() to provide a list of plots
#   to accomodate for the multi-gene set analyses that are produced. 

################################################################
################################################################

# Function to add all clustering information across all steps to metadata
# 1. determines all possible clusterings at increments of k + 1 with cutree.
# 2. generates matrix with all cluster assignments based on original cluster labels
# 3. for each step (increments of k + 1), add cluster assignments to meta data. 
# 4. returns Seurat object with phylo object in slot = "BuildClusterTree" and 
#    with all clustering assignments stored in meta data. 
recluster_all <- function(object, ident = "seurat_clusters", dims = NULL) {
  time.start <- Sys.time()
  object.phylo <- Tool(object, slot = "BuildClusterTree")
  
  
  orig_clusters <- object[[ident]][,1]
  if (ident != "seurat_clusters") {
    orig_clusters <- as.factor(orig_clusters-1)
  } # temporary fix on cluster number indexing for SC3 results
  
  object <- AddMetaData(object, metadata = orig_clusters, col.name = "ident")
  
  if (length(VariableFeatures(object)) == 0) object <- FindVariableFeatures(object)
  if (is.null(object.phylo)) {
    object <- BuildClusterTree(object, dims = dims, features = VariableFeatures(object), )
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

# same as stats::cutree(), but returns node number instead. 
cutree_nodes <- function(phylo_obj, k=NULL) {
  cut_output <- cutree(as.hclust(phylo_obj), k=k)
  tmp <- list()
  for (group in unique(cut_output)) {
    orig_identity <- names(cut_output)[cut_output == group] 
    all_orig_parents <- list()
    for (id in orig_identity) {
      all_orig_parents[[id]] <- GetParents(phylo_obj, node=as.numeric(id)+1)
    }
    tmp[[group]] <- Reduce(intersect, all_orig_parents)
  }
  unique_identities <- unique(unlist(tmp))
  tally <- c()
  for (id in unique_identities) {
    tally[as.character(id)] <- sum(unlist(tmp) == id)
  }
  node_ids <- names(tally)[tally==1]
  output <- cut_output
  for (orig_group in unique(cut_output)) {
    output[cut_output==orig_group] <- as.numeric(node_ids[orig_group])
  }
  return(output)
}

# Use (cluster/node number + 1) as a numeric to determine all parent clusters
# minus one from all of the outputs to get the proper cluster/node numbers.
GetParents <- function(tree, node) {
  output <- c(node)
  while (sum(tree$edge[,2]==node)>0) {
    node <- tree$edge[,1][tree$edge[,2] == node]
    output <- c(output, node)
  }
  return(rev(output))
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

# Subtracts mean expression per gene (given log-normalized expression)
NormalizeGenes <- function(object) {
  sub.vectors <- function(x, y) {
    nonzero <- which(x != 0)
    x[nonzero] <- x[nonzero] - y[nonzero] # where x and y are vectors of equal length
    return(x)
  }
  
  object <- NormalizeData(object, normalization.method = "LogNormalize", 
                          scale.factor = 10000, verbose = FALSE)
  expr_mat <- GetAssayData(object, slot = "data")
  means <- apply(expr_mat, MARGIN = 1, FUN = mean)
  expr_mat <- as.sparse(apply(expr_mat, MARGIN = 2, FUN = sub.vectors, y = means))
  
  object <- SetAssayData(object, slot = "data", new.data = expr_mat)
  return(object)
}

# Function to identify largest cluster/node (at lowest step) with lowest p-value.
# Returns a single cluster name. 
most_sig_cluster <- function(scan_output, gene_set) {
  score_matrix <- scan_output$scores
  output <- rownames(score_matrix)[score_matrix[,gene_set] == min(score_matrix[,gene_set])]
  n_tips <- min(as.numeric(scan_output$root_clusters))-1
  sig_clusters <- as.numeric(output)
  if (any(sig_clusters > n_tips)) {
    output <- min(sig_clusters[sig_clusters > n_tips])
  } else {
    # if all results are tips, choose largest cluster
    output <- min(sig_clusters)
  }
  return(as.character(output))
}

################################################################
################################################################
 
# Explore DEG between clusters determined from sig_cells()

# Generates a matrix of TRUE/FALSE with rownames set to cell ids and colnames
# set to clsuter/node number. Boolean indicates whether a cell is a member
# of a particular cluster/node or not.
# Returns a matrix.
cluster_bool_mat <- function(object) {
  concat_out <- concat_steps(object)
  cluster_ids <- sort(unique(concat_out$cluster))
  
  total_cells <- dim(object)[2]
  cell_cluster_mat <- matrix(data = NA, nrow = total_cells, ncol = length(cluster_ids))
  rownames(cell_cluster_mat) <- colnames(object)
  colnames(cell_cluster_mat) <- as.character(cluster_ids)
  
  for (cluster in cluster_ids) {
    cells_in_clust <- unique(concat_out$cell_id[concat_out$cluster == cluster])
    cell_cluster_mat[,as.character(cluster)] <- rownames(cell_cluster_mat) %in% cells_in_clust
  }
  return(cell_cluster_mat)
}

# Sets active.ident of Seurat object based on output from sig_cells
# parameter gene_set must be a column from sig_cells_output$cells
# cluster.labels must be a character vector in the same order as clusters to label 
# from parameters node or gene_set, and must therefore be the same length.
# Returns a Seurat object
SetSigIdent <- function(object, sig_cells_output, gene_set = NULL, 
                        nodes = NULL, label = "other", cluster.labels = NULL) {
  if (!is.null(gene_set)) sig_clusters <- sig_cells_output$clusters[[gene_set]]
  if (!is.null(nodes)) {
    if (exists("sig_clusters")) {
      sig_clusters <- append(sig_clusters, as.character(nodes), length(sig_clusters))
    } else {
      sig_clusters <- as.character(nodes)
    }
  }
  
  bool_mat <- cluster_bool_mat(object)
  new.ident <- matrix(data = NA, nrow = length(colnames(object)))
  for (cluster in 1:length(sig_clusters)) {
    cluster.number <- sig_clusters[cluster]
    if (is.null(cluster.labels)) {
      new.ident[bool_mat[,as.character(cluster.number)], ] <- cluster.number
    } else {
      new.ident[bool_mat[,as.character(cluster.number)], ] <- cluster.labels[cluster]
    }
  }
  
  new.ident[is.na(new.ident)] <- label
  rownames(new.ident) <- colnames(object)
  object <- SetIdent(object, cells = colnames(object), value = new.ident)
  return(object)
}

# Subset a Seurat object based on output from sig_cells
# parameter gene_set must be a column from sig_cells_output$cells
# Set gene_set OR node. 
# Calls subset_sigout to perform selection and phylo object subsetting. 
subset_object <- function(object, sig_cells_output, gene_set=NULL, node=NULL) {
  selection <- subset_sigout(sig_cells_output, node = node, gene_set = gene_set)
  object <- object[,selection$cells[,gene_set]]
  object@tools$BuildClusterTree <- selection$phylo
  return(object)
}

# Subset a sig_cells_output object based on a column in "cells" OR by any desired
# node from the tree. Setting gene_set will take the largest significant node
# from the phylo object, ignoring smaller nodes/clusters deemed to be significant.
subset_sigout <- function(sig_cells_output, node=NULL, gene_set=NULL) {
  node <- as.numeric(node)+1
  phylo <- sig_cells_output$phylo
  
  if (!is.null(gene_set)) {
    clusters <- sig_cells_output$clusters[[gene_set]]
    d <- lapply(clusters, FUN = phytools::getDescendants, tree = phylo)
    d <- unlist(lapply(d, FUN = length))
    names(d) <- clusters
    major.branch <- names(d)[d == max(d)]
    node <- major.branch
  } else {
    gene_set <- as.character(node-1)
  }
  
  tips <- ape::Ntip(phylo)
  if (node < tips) is.tip <- TRUE else is.tip <- FALSE
  if (is.tip) {
    phylo <- treeio::tree_subset(phylo, node = node, levels_back = 3)
  } else {
    phylo <- treeio::tree_subset(phylo, node = node, levels_back = 0)
  }
  
  select_cl <- c(phylo$node.label, phylo$tip.label)
  
  bool_mat <- sig_cells_output$cell_cluster_matrix[,select_cl]
  bool_mat <- bool_mat[apply(bool_mat, MARGIN = 1, FUN = any),]
  cells <- data.frame(row.names = rownames(sig_cells_output$cells))
  cells[,gene_set] <- rownames(sig_cells_output$cells) %in% rownames(bool_mat)
  scores <- sig_cells_output$scan_output$scores
  scores <- scores[scores$cluster %in% as.numeric(select_cl),]
  avg_expr <- sig_cells_output$scan_output$avg_expr
  avg_expr <- avg_expr[avg_expr$cluster %in% as.numeric(select_cl),]
  
  sig_cells_output$phylo <- phylo
  sig_cells_output$cells <- cells
  sig_cells_output$cell_cluster_matrix <- bool_mat
  sig_cells_output$scan_output$scores <- scores
  sig_cells_output$scan_output$avg_expr <- avg_expr

  return(sig_cells_output)
}

################################################################
################################################################

# Plotting functions

# Generates plots with UMAP projections of all clusters at every step
# Returns a list of ggplot objects; each element for each step.
cut_dimplots <- function(object, label.size = 6) {
  require(ggplot2)
  output <- list()
  for (i in 1:length(names(object@meta.data)[grep(names(object@meta.data), pattern = "step_")])) {
    output[[paste0("step_", i)]] <- DimPlot(object,
                                            group.by = paste0("step_", i),
                                            label = TRUE, 
                                            label.size = label.size,
                                            pt.size = 0.1) + 
      NoLegend() +
      ggtitle(paste("n_clusters =", i))
  }
  return(output)
}

# Given a Seurat object and the sig_cells() output, plot either the p-values
# or the cluster average expression values on a dendrogram.
# Edges are sized proportionally to cluster sizes. 
# Output is a list of ggplot objects. 
feature_phylo <- function(sig_cells_output, values = "mean", 
                          cap = 1e-50, title = NULL, subtitle = NULL) {
  phylo <- sig_cells_output$phylo
  scan_output <- sig_cells_output$scan_output
  marker_genes <- scan_output[5:length(scan_output)]
  output <- list()
  
  for (set in names(marker_genes)) {
    node_info <- cluster_stats(sig_cells_output)[[set]]
    edge_thickness <- (node_info$cluster_size)*5/(max(node_info$cluster_size))
    if (is.null(subtitle)) {target.name <- set} else {target.name <- subtitle}
    
    if (values == "pval") {
      if (!is.null(cap)) {node_info$score[node_info$score < cap] <- cap}
      base_tree <- ggtree::ggtree(phylo, aes(color = -log10(score)), 
                                  size = edge_thickness)
      base_tree <- ggtree::"%<+%"(base_tree, node_info)
      subtitle_final <- paste("Cluster -log10(pval) for", target.name, 
                              "\n p-value capped at", cap, "for plotting")
    }
    
    if (values == "mean") {
      base_tree <- ggtree::ggtree(phylo, aes(color = avg_expr), size = edge_thickness)
      base_tree <- ggtree::"%<+%"(base_tree, node_info)
      subtitle_final <- paste("Average", target.name, "expression per cluster")
    }
    
    if (is.null(title)) {title <- sig_cells_output$project.name} else {title <- title}
    
    final_tree <- base_tree +
      ggtree::geom_nodelab(geom = "label") +
      ggtree::geom_tiplab(geom = "text") +
      ggplot2::scale_color_gradient(low="grey", high="blue") +
      ggplot2::ggtitle(label = title, subtitle = subtitle_final)
    
    output[[set]] <- final_tree
  }
  return(output)
}

# Combines all of the "step_" meta.data slots from a Seurat object after recluster_all()
# has been run. the column "cell_id" will have as many repeats of every cell in a 
# data set as there are steps. 
# This makes it easier to identify which cells are in a particular cluster/node
concat_steps <- function(object) {
  metadata <- object@meta.data
  all_steps <- names(metadata)[grep(x = names(metadata), pattern = "step_")]
  cells <- c()
  clusters <- c()
  for (step in all_steps) {
    cells <- append(cells, values = colnames(object), after = length(cells))
    clusters <- append(clusters, values = as.numeric(as.character(metadata[[step]])), 
                       after = length(clusters))
  }
  output <- data.frame(row.names = 1:length(cells))
  output[,"cell_id"] <- cells
  output[,"cluster"] <- clusters
  return(output)
}

# Given the scan_dendrogram() output, construct a data.frame with node/cluster 
# number in first column, p-values in the second column, and mean expression 
# values in the third column. 
cluster_stats <- function(sig_cells_output) {
  output <- list()
  
  phylo <- sig_cells_output[["phylo"]]
  scan_output <- sig_cells_output[["scan_output"]]
  scores <- scan_output$scores
  expr <- scan_output$avg_expr
  
  all_clusters <- c(phylo$tip.label, phylo$node.label)
  bool_mat <- sig_cells_output[["cell_cluster_matrix"]]
  bool_mat <- bool_mat[,colnames(bool_mat) %in% all_clusters]
  all_sizes <- apply(bool_mat, MARGIN = 2, FUN = sum)
  
  marker_genes <- scan_output[5:length(scan_output)]
  
  for (set in names(marker_genes)) {
    cluster_score <- c()
    cluster_expr <- c()
    cluster_size <- c()
    
    for (cluster in all_clusters) {
      cluster_score <- append(cluster_score, scores[,set][rownames(scores) == cluster],
                              after = length(cluster_score))
      cluster_expr <- append(cluster_expr, expr[,set][rownames(expr) == cluster],
                             after = length(cluster_expr))
      cluster_size <- append(cluster_size, all_sizes[cluster],
                             after = length(cluster_size))
    }
    
    set_data <- data.frame(row.names = all_clusters)
    set_data[,"cluster"] <- all_clusters
    set_data[,"score"] <- cluster_score
    set_data[,"avg_expr"] <- cluster_expr
    set_data[,"cluster_size"] <- cluster_size
    
    output[[set]] <- set_data
  }
  return(output)
}

# Create an animated ggplot object (via gganimate::transition_states())
# animate.by should be a character vector of column names from object metadata.
DimPlot_animation <- function(object, reduction = "umap", color.by = NULL, 
                              animate.by, data.type = "continuous") {
  mat <- object@meta.data
  embeddings <- object@reductions[[reduction]]@cell.embeddings
  m <- mat[,animate.by]
  if (data.type == "continuous") {m <- apply(m, MARGIN = 2, FUN = as.numeric)}
  if (data.type == "discrete") {m <- apply(m, MARGIN = 2, FUN = as.character)}
  rownames(m) <- rownames(mat)
  
  df.all <- do.call(rbind, lapply(1:length(colnames(m)), function(i) {
    assignments <- m[rownames(embeddings), colnames(m)[i]]
    df <- data.frame(embeddings, assignments, metadata.colname=colnames(m)[i], order=i)
  }))

  x <- df.all[,1]
  y <- df.all[,2]

  plot <- ggplot(data = df.all, aes(x = x, y = y, color = assignments)) +
    geom_point() +
    theme_bw() +
    gganimate::transition_states(order, transition_length = 1, state_length = 3)

  return(plot)
}

