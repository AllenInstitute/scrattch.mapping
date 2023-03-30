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
#' @export
buildTaxonomy = function(counts,
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
  if(is.null(umap.coords)){stop("Compute UMAP dimensions and supply umap.coords")}
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
  ## norm.data.t = Matrix::t(tpm.matrix)
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

#' Add marker genes to reference dendrogram for tree mapping
#'
#' @param dend A dendrogram in R format to which marker genes will be added, or a character string with a file location of "dend.RData"
#' @param norm.data A matrix of log normalized reference data, or character string with a file location of "data_t.feather".  If a count matrix is provided, the data data will be log normalized.  This should be the data matrix used to generate the dendrogram.
#' @param metadata Data frame of metadata with rows corresponding to cells/nuclei, and either row names or a column called "sample_id" corresponding to cell names. This matrix must include entries for all cells in norm.data.  Could also be a file.  Columns can be numeric, categorical, or factors.
#' @param celltypeColumn Column name correspond to the cell type names in the dendrogram (default = "cluster_label"). At least two cells per cell type in the dendrogram must be included.
#' @param subsample The number of cells to retain per cluster (default = 100)
#' @param num.markers The maximum number of markers to calculate per pairwise differential calculation per direction (default = 20)
#' @param de.param Differential expression (DE) parameters for genes and clusters used to define marker genes.  By default the values are set to the 10x nuclei defaults from scrattch.hicat, except with min.cells=2 (see the function `de_param` in the scrattch.hicat for more details).
#' @param calculate.de.genes Default=TRUE. If set to false, the function will search for a file called "de.genes.rda" to load precalculated de genes.  
#' @param save.shiny.output Should standard output files be generated and saved to the directory (default=TRUE).  These are not required for tree mapping, but are required for building a patch-seq shiny instance.  This is only tested in a UNIX environment.  See notes.
#' @param mc.cores Number of cores to use for running this function to speed things up.  Default = 1.  Values>1 are only supported in an UNIX environment and require `foreach` and `doParallel` R libraries.
#' @param bs.num,p,low.th Extra variables for the `map_dend_membership` function in scrattch.hicat.  Defaults are set reasonably.
#' @param shinyFolder The location to save shiny output, if desired
#'
#' NOTES
#'
#' If save.shiny.output=TRUE, the following files will be generated:
#'   reference.rda, which includes a variable `reference` as follows:
#'       reference$cl.dat - These are the cluster means that are used for mapping comparisons
#'       reference$dend   - This is the dendrogram with marker genes attached
#'   membership_information_reference.rda, which includes two variables
#'       `memb.ref`   - matrix indicating how much confusion there is the mapping between each cell all of the nodes in the tree (including all cell types) when comparing clustering and mapping results with various subsamplings of the data
#'       `map.df.ref` - Result of tree mapping for each cell in the reference against the clustering tree, including various statistics and marker gene evidence.  This is the same output that comes from tree mapping.#'
#' 
#' @import feather
#' @import scrattch.hicat
#' @import MatrixGenerics
#' @import dendextend
#'
#' @return An updated dendrogram variable that is the same as `dend` except with marker genes added to each node.
#'
#' @export
addDendrogramMarkers = function(dend,
                                norm.data,
                                metadata,
                                celltypeColumn = "cluster_label",
                                subsample = 100,
                                num.markers = 20,
                                de.param=scrattch.hicat::de_param(low.th = 1,  # Recommended values for 10x Nuclei
                                                                  padj.th = 0.01,
                                                                  lfc.th = 1,
                                                                  q1.th = 0.3,
                                                                  q2.th = NULL,
                                                                  q.diff.th = 0.7,
                                                                  de.score.th = 100,
                                                                  min.cells = 2,
                                                                  min.genes = 5),
                                calculate.de.genes = TRUE,
                                save.shiny.output = TRUE,
                                mc.cores=1, 
                                bs.num=100, p=0.8, low.th=0.1,
                                shinyFolder = paste0(getwd(),"/")
){

  ## Checks and data formatting
  # dend
  if(class(dend)=="character"){
    if(file.exists(dend)){
      dend <- readRDS(dend)
    } else {
      stop(paste(dend,"is not a valid filename."))
    }
  } else {
    if(class(dend)!="dendrogram"){stop("dend must be a dendrogram format or a file name.")}
  }
  
  #norm.data
  if(class(norm.data)=="character"){
    if(file.exists(norm.data)){
      data_t <- try({feather(norm.data)})
      if(class(data_t)=="try-error"){stop("norm.data must be in feather format with a gene column + data matrix included.")}
      norm.data <- as.matrix(data_t[,colnames(data_t)!="gene"])
      if(!is.element(class(norm.data[,1]),c("numeric","integer"))){stop("norm.data doesn't have numeric entries")}
      genes <- try(data_t$gene)
      if(class(genes)=="try-error"){stop("norm.data must have a column called gene with gene names.")}
      rownames(norm.data) <- genes
    } else {
      stop(paste(norm.data,"is not a valid filename."))
    }
  } else {
    if(!grepl("atrix",as.character(class(norm.data)))){stop("norm.data must be some kind of matrix or a file name.")}
  }
  # Assume if norm.data has a value>50 it needs to be logCPM normalized
  if(!is.element(class(norm.data[,1]),c("numeric","integer"))){stop("norm.data doesn't have numeric entries")}
  if(max(norm.data)>50) norm.data <- logCPM(norm.data)
  
  #metadata 
  if(class(metadata)=="character"){
    if(file.exists(metadata)){
      metadata <- try({feather(metadata)})
      if(class(data_t)=="try-error"){stop("metadata file must be in feather format.")}
      metadata <- as.data.frame(metadata)
    } else {
      stop(paste(metadata,"is not a valid filename."))
    }
  } else {
    if(class(metadata)!="data.frame"){stop("metadata must be a data frame or a file name.")}
  }
  if(is.element("sample_id",colnames(metadata))) rownames(metadata) <- metadata$sample_id
  if(length(intersect(rownames(metadata),colnames(norm.data)))==0){stop("Metadata sample ids (or row names) do not match sample ids in the data (column names).")}
  if(length(setdiff(colnames(norm.data),rownames(metadata)))>0){stop("Some samples are missing in the metadata file based on matching sample ids between norm.data columns and metadata rows.")}
  metadata <- metadata[colnames(norm.data),]
  
  #celltype labels
  if(!is.element(celltypeColumn,colnames(metadata))){stop(paste(celltypeColumn,"is not a column in the metadata data frame."))}
  cluster.vector <- setNames(metadata[,celltypeColumn],rownames(metadata))
  if(length(intersect(labels(dend),cluster.vector))==0){stop(paste(celltypeColumn, "data does not match dendrogram labels."))}
  if(length(setdiff(labels(dend),cluster.vector))>0){stop("Some clusters in the dendrogram are not include in the reference metadata file.")}
  keep.samples   <- subsampleCells(cluster.vector,subsample)&is.element(cluster.vector,labels(dend))
  norm.data      <- norm.data[,keep.samples]
  metadata       <- metadata[keep.samples,]
  cluster.vector <- factor(cluster.vector[keep.samples],levels=labels(dend))
  
  # Ensure directory exists, if not create it
  dir.create(shinyFolder, showWarnings = FALSE)
  
  
  ## Define and collect marker genes
  print("Define some relevant variables")
  cl.df     <- as.data.frame(metadata[match(labels(dend),cluster.vector),])
  cl.df$cluster_label <- cl.df[,celltypeColumn]
  rownames(cl.df) <- 1:length(labels(dend)) 
  cl.label  <- as.factor(setNames(cl.df$cluster_label, rownames(cl.df)))
  select.cl <- droplevels(as.factor(setNames(match(metadata[,celltypeColumn],cl.label), metadata$sample_id)))
  
  
  ## CHECK IF THIS IS NEEDED
  # We might need to relabel the dendrogram from 1 to #clusters in order
  dend_ref = dend
  labels(dend) = names(cl.label)[match(labels(dend),cl.label)]
  
  print("Define marker genes and gene scores for the tree")
  if((!file.exists(file.path(shinyFolder,"de.genes.rda")))|calculate.de.genes){
    print("=== NOTE: This step can be very slow (several minute to many hours).")
    print("      To speed up the calculation (or if it crashes) try decreasing the value of subsample.")
    de.genes = select_markers(norm.dat=norm.data, cl=select.cl, n.markers=num.markers, 
                              de.param = de.param, de.genes = NULL)$de.genes
    save(de.genes, file=file.path(shinyFolder,"de.genes.rda"))  
  } else {
    load(file.path(shinyFolder,"de.genes.rda"))
  }
  min.marker.gene.count <- as.numeric(as.character(lapply(de.genes, function(x) x$num)))
  if(sum(min.marker.gene.count<2)>0)({
    stop("Marker genes could not be calculated for at least one node in the tree. Tree mapping will not work in this situation. We recommend loosening the de.param parameters and trying again, but warn that some cell types may not be well resolved.")
  })
  # Note: this is not a robust way to address this issue!
  print("Calculate gene scores")
  gene.score = get_gene_score(de.genes)
  
  
  print("Build the reference dendrogram")
  invisible(capture.output({  # Avoid printing lots of numbers to the screen
    reference = build_reference(cl=select.cl, norm.dat=norm.data, dend=dend, de.genes=de.genes, 
                                cl.label=cl.label, up.gene.score=gene.score$up.gene.score, 
                                down.gene.score=gene.score$down.gene.score, n.markers=30)
  }))
  labels(reference$dend) <- setNames(colnames(reference$cl.dat),colnames(reference$cl.dat))
  # There is an error introduced somewhere in build_reference, which the above line fixes -- I'm not sure if this is needed and if it is, I should wrap it into build_reference
  
  
  if(sum(!is.na(get_nodes_attr(reference$dend, "original_label"))>0)){
    print("This section is needed if the starting dendrogram is from the nomenclature GitHub ")
    reference$dend <- revert_dend_label(reference$dend,get_nodes_attr(reference$dend, "original_label"),"label")
  }
  
  if(save.shiny.output){
    print("Save the reference dendrogram")
    save(reference, file=file.path(shinyFolder,"reference.rda")) # Redundant
    saveRDS(dend, file.path(shinyFolder,"dend.RData"))        # Redundant

    
    print("Build membership table of reference vs. reference for use with patch-seq mapping")
    cl.dat     <- reference$cl.dat
    dend       <- reference$dend
    invisible(capture.output({  # Avoid printing lots of numbers to the screen
      memb.ref   <- map_dend_membership(dend, cl.dat, map.dat=norm.data, map.cells=names(select.cl),
                                        mc.cores=mc.cores, bs.num=bs.num, p=p, low.th=low.th)
      map.df.ref <- summarize_cl(dend, memb.ref, norm.data)
    }))
    memb.ref   <- memb.ref[metadata$sample_id,]
    map.df.ref <- map.df.ref[metadata$sample_id,]
    save(memb.ref, map.df.ref, file=file.path(shinyFolder,"membership_information_reference.rda"))
  }
  
  return(reference$dend)  
}