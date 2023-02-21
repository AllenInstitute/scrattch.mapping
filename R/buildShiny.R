#' This function builds the minimum files required for Shiny
#'
#' @param counts A count matrix in sparse format: dgCMatrix.
#' @param meta.data Meta.data corresponding to count matrix. Rownames must be equal to colnames of counts.
#' @param feature.set Set of feature used to calculate dendrogram. Typically highly variable and/or marker genes.
#' @param umap.coords Dimensionality reduction coordiant data.frame with 2 columns. Rownames must be equal to colnames of counts.
#' @param shinyFolder The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"
#' @param cluster_colors An optional named character vector where the values correspond to colors and the names correspond to celltypes in celltypeColumn.  If this vector is incomplete, a warning is thrown and it is ignored. 
#' @param metadata_names An optional named character vector where the vector NAMES correspond to columns in the metadata matrix and the vector VALUES correspond to how these metadata should be displayed in Shiny. This is used for writing the desc.feather file later.
#' @param subsample The number of cells to retain per cluster
#' @param reorder.dendrogram Should dendogram attempt to match a preset order? (Default = FALSE).  If TRUE, the dendrogram attempts to match the celltype factor order as closely as possible (if celltype is a character vector rather than a factor, this will sort clusters alphabetically, which is not ideal).
#' 
#' @import scrattch.hicat
#' @import scrattch.io
#' @import feather
#' @import tibble
#' @import dplyr
#' @import Matrix
#' @import pvclust
#' 
#' @return
#'
#' @export
buildReference = function(counts,
                          meta.data,
                          feature.set,
                          umap.coords,
                          shinyFolder,
                          cluster_colors = NULL,
                          metadata_names = NULL,
                          subsample=2000,
                          reorder.dendrogram = FALSE){

  ## Checks
  if(!"cluster" %in% colnames(meta.data)){stop("cluster must be defined in the meta.data object")}
  if(is.null(feature.set)){stop("Compute variable features and supply feature.set")}
  if(!all(colnames(counts) == rownames(meta.data))){stop("Colnames of `counts` and rownames of `meta.data` do not match.")}
	
  ## Ensure directory exists, if not create it
  dir.create(shinyFolder, showWarnings = FALSE)

  ## Subsample nuclei per cluster, max
  if((subsample > 0)&(subsample < Inf)){
      kpSub = colnames(counts)[subsampleCells(meta.data$cluster, subsample)]
  }else{
      kpSub = colnames(counts)
  }

  ## Get the data and metadata matrices
  counts    = counts[,kpSub]
  meta.data = meta.data[kpSub,]

  ## Create log TPM matrix from counts
  tpm.matrix = scrattch.hicat::logCPM(counts)

  ## Next, generate and output the data matrices. Make sure matrix class is dgCMatrix and not dgTMatrix.  
  sample_id = colnames(tpm.matrix); meta.data$sample_id = sample_id
  gene      = rownames(tpm.matrix)
  cluster = meta.data$cluster; names(cluster) = rownames(meta.data) ## For get_cl_medians

  ## Counts feather - Not used for shiny, but useful for saving raw data nonetheless
  print("===== Building counts.feather =====")
  counts.tibble = as_tibble(counts)
  counts.tibble = cbind(gene, counts.tibble)
  counts.tibble = as_tibble(counts.tibble)
  write_feather(counts.tibble, file.path(shinyFolder, "counts.feather"))

  ## Data_t feather
  print("===== Building data_t.feather =====")
  data.tibble = as_tibble(tpm.matrix)
  data.tibble = cbind(gene, data.tibble)
  data.tibble = as_tibble(data.tibble)
  write_feather(data.tibble, file.path(shinyFolder, "data_t.feather"))
  
  ## Data feather
  print("===== Building data.feather =====")
  norm.data.t = t(as.matrix(tpm.matrix))
  norm.data.t = as_tibble(norm.data.t)
  norm.data.t = cbind(sample_id, norm.data.t)
  norm.data.t = as_tibble(norm.data.t)
  write_feather(norm.data.t, file.path(shinyFolder, "data.feather"))

  ## ----------
  ## Run auto_annotate, this changes sample_id to sample_name.
  print("===== Format metadata table for R shiny folder =====")
  meta.data = scrattch.io::auto_annotate(meta.data)

  ### Convert chars and factors to characters (moved to AFTER auto_annotate, so numbers can be assigned in correct order for factors)
  for (col in colnames(meta.data)){ 
      if(is.character(meta.data[,col]) | is.factor(meta.data[,col])){
          meta.data[,col] = as.character(meta.data[,col])
      }
  }
  # Shouldn't be needed as factors are okay, and actually suggested for cluster order
  
  ## Varibow color set is broken -- this will fix it
  for (col in which(grepl("_color",colnames(meta.data)))){
      kp = nchar(meta.data[,col])==5
      meta.data[kp,col] = paste0(meta.data[kp,col],"FF")
  }

  ## Adjust the cluster colors to match cluster_colors, if available. 
  if(!is.null(cluster_colors)){
    if(length(setdiff(annotations$cluster,names(cluster_colors)))>0){
      warning("cluster_colors is not a named vector with colors for every cluster and will therefore be ignored.")
    } else {
      meta.data$cluster_color <- as.character(cluster_colors[meta.data$cluster_label])
    }
  }
    
  print("===== Building dendrogram =====")
  ## Get cluster medians
  medianExpr = get_cl_medians(tpm.matrix, cluster) 

  ## Define the cluster info 
  unique.meta.data = meta.data %>% distinct(cluster_id, 
                                            cluster_label, 
                                            cluster_color)
  rownames(unique.meta.data) = unique.meta.data$cluster_label

  ## Dendrogram parameters and gene sets
  use.color = setNames(unique.meta.data$cluster_color, unique.meta.data$cluster_label)[colnames(medianExpr)]
  l.rank    = NULL
  if(reorder.dendrogram){
    l.rank = setNames(meta.data$cluster_id[match(unique.meta.data$cluster_label, meta.data$cluster_label)], unique.meta.data$cluster_label)
    l.rank = sort(l.rank)
  }
  
  ##
  invisible(capture.output({  # Avoid printing lots of numbers to the screen
    dend.result = scrattch.mapping::build_dend(
      cl.dat  = medianExpr[feature.set,],
      cl.cor  = NULL,
      l.color = use.color,
      l.rank  = l.rank, 
      nboot   = 1,
      ncores  = 1)
  }))
  
  ## Output tree
  dend = dend.result$dend
  saveRDS(dend, file.path(shinyFolder,"dend.RData"))

  ## Output tree order
  outDend = data.frame(cluster = labels(dend), order = 1:length(labels(dend)))
  write.csv(outDend, file.path(shinyFolder, "ordered_clusters.csv"))
    
  ## Write the desc file.  
  anno_desc = create_desc(meta.data, use_label_columns = TRUE)
  # Subset the desc file to match metadata_names, if provided
  if(!is.null(metadata_names)){
    desc <- anno_desc[match(names(metadata_names), as.character(as.matrix(anno_desc[,1]))),]
    desc[,2] <- as.character(metadata_names[as.character(as.matrix(desc[,1]))])
    desc <- desc[!is.na(desc$base),]  # Remove missing values
    anno_desc <- desc
  }
  write_feather(anno_desc, file.path(shinyFolder,"desc.feather"))

  ## Minor reformatting of metadata file, then write metadata file
  meta.data$cluster = meta.data$cluster_label; 
  colnames(meta.data)[colnames(meta.data)=="sample_name"] <- "sample_name_old" # Rename "sample_name" to avoid shiny crashing
  if(!is.element("sample_id", colnames(meta.data))){ meta.data$sample_id = meta.data$sample_name_old } ## Sanity check for sample_id
  meta.data$cluster_id <- as.numeric(factor(meta.data$cluster_label,levels=labels(dend))) # Reorder cluster ids to match dendrogram
  write_feather(meta.data, file.path(shinyFolder,"anno.feather"))

  ## Write the UMAP coordinates.  
  print("===== Building umap/tsne feathers (precalculated) =====")
  tsne      = data.frame(sample_id = rownames(umap.coords),
                         all_x = umap.coords[,1],
                         all_y = umap.coords[,2])
  tsne      = tsne[match(meta.data$sample_id, tsne$sample_id),]
  tsne_desc = data.frame(base = "all",
                         name = "All Cells UMAP")

  write_feather(tsne, file.path(shinyFolder,"tsne.feather"))
  write_feather(tsne_desc, file.path(shinyFolder,"tsne_desc.feather"))

  ##
  print("===== Building count, median, sum feathers =====")
  all_clusters = unique(meta.data$cluster_id)
  all_clusters = all_clusters[order(all_clusters)]
  allClust = paste0("cluster_", all_clusters)
  count_gt0 = matrix(0, ncol = length(all_clusters), nrow = nrow(data.tibble))
  count_gt1 = sums = medianmat = count_gt0

  ## Compute the number of genes with greater than 0,1 counts, gene sums and medians
  for (i in 1:length(all_clusters)) {
    cluster         = all_clusters[i]
    cluster_samples = which(meta.data$cluster_id == cluster)
    cluster_data    = tpm.matrix[,cluster_samples]
    cluster_counts  = counts[,colnames(cluster_data)]
    count_gt0[, i]  = Matrix::rowSums(cluster_counts > 0)
    count_gt1[, i]  = Matrix::rowSums(cluster_counts > 1)
    sums[, i]       = Matrix::rowSums(cluster_counts)
    medianmat[, i]  = apply(cluster_data, 1, median)
  }

  ##
  colnames(count_gt0) = colnames(count_gt1) = colnames(sums) = colnames(medianmat) = allClust
  count_gt0 = cbind(gene = gene, as.data.frame(count_gt0))
  count_gt1 = cbind(gene = gene, as.data.frame(count_gt1))
  sums = cbind(gene = gene, as.data.frame(sums))
  medianmat = cbind(gene = gene, as.data.frame(medianmat))

  ##
  count_n = meta.data %>% 
              arrange(cluster_id) %>% 
              group_by(cluster_id) %>% 
              summarise(n_cells = n())

  ##
  write_feather(count_gt0, file.path(shinyFolder,"count_gt0.feather"))
  write_feather(count_gt1, file.path(shinyFolder,"count_gt1.feather"))
  write_feather(count_n,   file.path(shinyFolder,"count_n.feather"))
  write_feather(medianmat, file.path(shinyFolder,"medians.feather"))
  write_feather(sums,      file.path(shinyFolder,"sums.feather"))
}