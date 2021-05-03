# Functions to help analyze biological differences between nodes in a Seurat object
# Explore DEG between clusters determined from sig_cells()

# Sets active.ident of Seurat object based on output from sig_cells
# parameter gene_set must be a column from sig_cells_output$cells
# cluster.labels must be a character vector in the same order as clusters to label 
# from parameters node or gene_set, and must therefore be the same length.
# Returns a Seurat object
SetSigIdent <- function(object, sig_cells_output = NULL, gene_set = NULL, 
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

# Perform GSEA between different nodes from a Seurat object
# Can input a vector of node numbers into ident.2 
gsea_nodes <- function(object, pathways, ident.1, ident.2) {
  object <- SetSigIdent(object, nodes = c(ident.1, ident.2))
  
  if (length(ident.2) > 1) {
    ident = paste(ident.2, collapse = "_")
    labels = levels(object@active.ident)
    labels[labels %in% ident.2] <- ident
  } else {
    ident <- ident.2
  }
  
  object <- AddMetaData(object, metadata = object@active.ident, col.name = "selection")
  stats = presto::wilcoxauc(object, "selection", groups_use = c(ident.1, ident))
  
  stats <- stats %>% 
    dplyr::filter(group==ident.1) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  ranked = tibble::deframe(stats)
  
  output = fgsea::fgsea(pathways = pathways, stats = ranked)
  return(output)
}

