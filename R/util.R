# Utility functions to help the main_functions.R run the hierarchical analysis

# Modified BuildClusterTree function from Seurat
# Added "dist.method" parameter to indicate euclidean or correlation distance
# to calculate distance between centroids.
BuildClusterTree <- function(
  object,
  assay = NULL,
  features = NULL,
  dims = NULL,
  graph = NULL,
  slot = 'data',
  reorder = FALSE,
  reorder.numeric = FALSE,
  verbose = TRUE,
  dist.method = "correlation",
  linkage = "average"
) {
  require(stats)
  require(utils)
  PackageCheck <- function(..., error = TRUE) {
    pkgs <- unlist(x = c(...), use.names = FALSE)
    package.installed <- vapply(
      X = pkgs,
      FUN = requireNamespace,
      FUN.VALUE = logical(length = 1L),
      quietly = TRUE
    )
    if (error && any(!package.installed)) {
      stop(
        "Cannot find the following packages: ",
        paste(pkgs[!package.installed], collapse = ', '),
        ". Please install"
      )
    }
    invisible(x = package.installed)
  }
  if (!PackageCheck('ape', error = FALSE)) {
    stop(cluster.ape, call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = graph)) { 
    idents <- levels(x = object)
    nclusters <- length(x = idents)
    data.dist <- matrix(
      data = numeric(length = 1L),
      nrow = nclusters,
      ncol = nclusters,
      dimnames = list(idents, idents)
    )
    graph <- object[[graph]]
    cxi <- CellsByIdentities(object = object)
    cpairs <- na.omit(object = unique(x = t(x = apply(
      X = expand.grid(1:nclusters, 1:nclusters)[, c(2, 1)],
      MARGIN = 1,
      FUN = function(x) {
        if (length(x = x) == length(x = unique(x = x))) {
          return(sort(x = x))
        }
        return(c(NA, NA))
      }
    ))))
    if (verbose) {
      pb <- txtProgressBar(style = 3, file = stderr())
    }
    for (i in 1:nrow(x = cpairs)) {
      i1 <- cpairs[i, ][1]
      i2 <- cpairs[i, ][2]
      graph.sub <- graph[cxi[[idents[i1]]], cxi[[idents[i2]]]]
      d <- mean(x = graph.sub)
      if (is.na(x = d)) {
        d <- 0
      }
      data.dist[i1, i2] <- d
      if (verbose) {
        setTxtProgressBar(pb = pb, value = i / nrow(x = cpairs))
      }
    }
    if (verbose) {
      close(con = pb)
    }
    diag(x = data.dist) <- 1
    data.dist <- dist(x = data.dist)
  } else if (!is.null(x = dims)) { 
    my.lapply <- ifelse(test = verbose, yes = pbapply::pblapply, no = lapply)
    embeddings <- Embeddings(object = object, reduction = 'pca')[, dims]
    data.dims <- my.lapply(
      X = levels(x = object),
      FUN = function(x) {
        cells <- WhichCells(object = object, idents = x)
        if (length(x = cells) == 1) {
          cells <- c(cells, cells)
        }
        temp <- colMeans(x = embeddings[cells, ])
      }
    )
    data.dims <- do.call(what = 'cbind', args = data.dims)
    colnames(x = data.dims) <- levels(x = object)
    data.dist <- dist(x = t(x = data.dims))
  } else {
    features <- features %||% VariableFeatures(object = object)
    features <- intersect(x = features, y = rownames(x = object))
    data.avg <- AverageExpression( # calculating centroids
      object = object,
      assays = assay,
      features = features,
      slot = slot,
      verbose = verbose
    )[[1]]
    if (dist.method == "euclidean") {
      data.dist <- dist(x = t(x = data.avg[features, ])) 
    } else if (dist.method == "correlation") {
      data.dist <- as.dist(CorDist(t(data.avg), t(data.avg), method = "pearson"))
    }
  }
  
  data.tree <- ape::as.phylo(x = hclust(d = data.dist, method = linkage))
  Tool(object = object) <- data.tree
  object@tools$DistanceMatrix <- data.dist
  
  if (reorder) {
    if (verbose) {
      message("Reordering identity classes and rebuilding tree")
    }
    old.ident.order <- levels(x = object)
    data.tree <- Tool(object = object, slot = 'BuildClusterTree')
    all.desc <- GetDescendants(tree = data.tree, node = (data.tree$Nnode + 2))
    all.desc <- old.ident.order[all.desc[all.desc <= (data.tree$Nnode + 1)]]
    Idents(object = object) <- factor(x = Idents(object = object), levels = all.desc, ordered = TRUE)
    if (reorder.numeric) {
      new.levels <- sort(x = unique(x = as.integer(x = Idents(object = object))))
      Idents(object = object) <- factor(x = as.integer(x = Idents(object = object)), levels = new.levels)
      object[['tree.ident']] <- as.integer(x = Idents(object = object))
    }
    object <- BuildClusterTree(
      object = object,
      assay = assay,
      features = features,
      dims = dims,
      graph = graph,
      slot = slot,
      reorder = FALSE,
      verbose = verbose
    )
  }
  return(object)
}

# This is to be re-named to BaseClust after preliminary tests
BuildClusterTree <- function(object, space = "pca", dims = 50, 
                      distance = "euclidean", linkage = "average") {
  time.start <- Sys.time()
  if (space == "pca") {
    d <- dist(object@reductions$pca@cell.embeddings[,1:dims], method = distance)
  } else if (space == "umap") {
    d <- dist(object@reductions$umap@cell.embeddings, method = distance)
  } else if (space == "variable_features") {
    m = GetAssayData(object, slot = "data")
    d <- dist(t(as.matrix(m[VariableFeatures(object),])), method = distance)
  }
  
  tree <- hclust(d = d, method = linkage)
  
  time.end <- Sys.time()
  runtime <- time.end - time.start
  print(runtime)
  
  Tool(object) <- tree
  object@tools$DistanceMatrix <- d
  object@tools$hclust_runtime <- runtime
  
  return(object)
}

# try to manually merge nodes 8 and 9
test = tree
# name new branch "8"
new.lab = test$labels[!(test$labels==9)]
new.order = test$order[!(test$labels==9)]

test$labels = new.lab
test$order = new.order
m = test$merge[-5,]
m[6,2] = -17 # originally at row 7, but now 6 because 5 was removed.
m[m > 5] <- m[m > 5] -length(5)
test$merge = m
test$height = test$height[-5]


# subtree function in progress
subtree <- function(tree, h = NULL , k = NULL) {
  cut = cutree(tree, h = h, k = k)
  branches.to.make = names(table(cut))[table(cut) > 1]
  
  # for each branch to make
  for (i in branches.to.make) {
    print(i)
    labels.to.combine = cut == as.numeric(i) # i.e. node 8 and 9
    branch.name = min(as.numeric(tree$labels[labels.to.combine]))
    indices = tree$order[labels.to.combine]

    dim1 = dim(tree$merge)[1]
    boolmat = matrix(tree$merge %in% -indices, nrow = dim1, ncol = 2)
    merge.rows = which(rowSums(boolmat) > 0)
    to.remove = merge.rows
    
    # Climb tree from bottom to top of branch 
    # identify indices of all merges
    while (length(merge.rows) > 1) {
      to.remove = append(to.remove, values = merge.rows, after = length(to.remove))
      boolmat = matrix(tree$merge %in% merge.rows, nrow = dim1, ncol = 2)
      merge.rows = which(rowSums(boolmat) > 0)
    }

    new.label = min(as.numeric(names(cut)[labels.to.combine]))
    keep.labels = as.logical((names(cut) == new.label) + !labels.to.combine)

    print(merge.rows)
    new.height = max(tree$height[merge.rows])
    tree$merge = tree$merge[-to.remove,]
    tree$height = tree$height[-to.remove]
    tree$labels = tree$labels[keep.labels]
    tree$order = tree$order[!(tree$order %in% names(cut)[keep.labels])]
  }
  return(tree)
}

# Calculates correlation distance and returns a matrix. 
CorDist <- function(points1, points2, method = "pearson") {
  output <- 1-(cor(t(points1), t(points2), method = method))
  return(output)
}

# Function to automatically go through typical Seurat analysis pipeline
Seurat_pipeline <- function(object, mito_filter = TRUE, filter = FALSE, 
                            n_var_feat = 2000, scale = FALSE, center = TRUE,
                            cluster_res = 0.8, nn_space = "pca", p_components = 20,
                            tsne = TRUE, umap = TRUE, vis = "pca") {
  require(Seurat)
  
  if (mito_filter) {
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    object <- subset(object, subset = percent.mt < 5)
    print("filtering out cells with high mitochondrial gene content")
  } 
  if (filter) {
    object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
    print("filtering out cells with too low or high feature counts")
  }
  
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = n_var_feat)
  all.genes <- rownames(object)
  object <- ScaleData(object, features = all.genes, do.scale = scale, do.center = center)
  object <- RunPCA(object, features = VariableFeatures(object))
  
  if (nn_space == "pca") {
    object <- FindNeighbors(object, dims = 1:p_components)
    object <- FindClusters(object, resolution = cluster_res)
  }
  if (nn_space == "variable_features") {
    object <- FindNeighbors(object, features = VariableFeatures(object))
    object <- FindClusters(object, resolution = cluster_res)
  }
  
  if (vis == "pca") {
    if (umap) object <- RunUMAP(object, dims = 1:p_components)
    if (tsne) object <- RunTSNE(object, dims = 1:p_components, check_duplicates = FALSE)
  }
  
  if (vis == "variable_features") {
    if (umap) object <- RunUMAP(object, features = VariableFeatures(object), 
                                dims = NULL)
    if (tsne) object <- RunTSNE(object, features = VariableFeatures(object), 
                                dims = NULL, check_duplicates = FALSE)
  }
  
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
  tally <- table(unlist(tmp))
  node_ids <- names(tally)[tally==1] # identify all tips
  
  output <- cut_output
  for (orig_group in unique(cut_output)) {
    output[cut_output==orig_group] <- as.numeric(node_ids[orig_group])
  }

  ts <- list()
  ts[["cut_output"]] <- cut_output
  ts[["tally"]] <- tally
  ts[["node_ids"]] <- node_ids
  ts[["tmp"]] <- tmp
  return(ts)
  # return(output)
}

# Use (node number + 1) as a numeric to determine all parent clusters
# minus one from all of the outputs to get the proper cluster/node numbers.
GetParents <- function(tree, node) {
  output <- c(node)
  while (sum(tree$edge[,2]==node)>0) {
    node <- tree$edge[,1][tree$edge[,2] == node]
    output <- c(output, node)
  }
  return(rev(output))
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

