# Visualization functions to create ggtree dendrogram or UMAP plots along the
# dendrogram of a Seurat object


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
  require(ggtree)
  
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
      base_tree <- ggtree(phylo, aes(color = -log10(score)), 
                                  size = edge_thickness)
      base_tree <- "%<+%"(base_tree, node_info)
      subtitle_final <- paste("Cluster -log10(pval) for", target.name, 
                              "\n p-value capped at", cap, "for plotting")
    }
    
    if (values == "mean") {
      base_tree <- ggtree(phylo, aes(color = avg_expr), size = edge_thickness)
      base_tree <- "%<+%"(base_tree, node_info)
      subtitle_final <- paste("Average", target.name, "expression per cluster")
    }
    
    if (is.null(title)) {title <- sig_cells_output$project.name} else {title <- title}
    
    final_tree <- base_tree +
      geom_nodelab(geom = "label") +
      geom_tiplab(geom = "text") +
      ggplot2::scale_color_gradient(low="grey", high="blue") +
      ggplot2::ggtitle(label = title, subtitle = subtitle_final)
    
    output[[set]] <- final_tree
  }
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

