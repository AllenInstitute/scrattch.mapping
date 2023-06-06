#' This function builds the minimum files required for Shiny
#'
#' @param counts A count matrix in sparse format: dgCMatrix.
#' @param meta.data Meta.data corresponding to count matrix. Rownames must be equal to colnames of counts.
#' @param feature.set Set of feature used to calculate dendrogram. Typically highly variable and/or marker genes.
#' @param umap.coords Dimensionality reduction coordiant data.frame with 2 columns. Rownames must be equal to colnames of counts.
#' @param taxonomyDir The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"
#' @param taxonomyName The name to assign for the Taxonomy h5ad
#' @param cluster_colors An optional named character vector where the values correspond to colors and the names correspond to celltypes in celltypeColumn.  If this vector is incomplete, a warning is thrown and it is ignored. 
#' @param metadata_names An optional named character vector where the vector NAMES correspond to columns in the metadata matrix and the vector VALUES correspond to how these metadata should be displayed in Shiny. This is used for writing the desc.feather file later.
#' @param subsample The number of cells to retain per cluster
#' @param reorder.dendrogram Should dendogram attempt to match a preset order? (Default = FALSE).  If TRUE, the dendrogram attempts to match the celltype factor order as closely as possible (if celltype is a character vector rather than a factor, this will sort clusters alphabetically, which is not ideal).
#' @param dend Existing dendrogram associated with this taxonomy (e.g., one calculated elsewhere).  If NULL (default) a new dendrogram will be calculated based on the input `feature.set`
#' 
#' @import scrattch.hicat 
#' @import scrattch.io
#' @import feather
#' @import tibble
#' @import dplyr
#' @import Matrix
#' @import pvclust
#'
#' @return AIT anndata
#'
#' @export
buildTaxonomy = function(counts,
                         meta.data,
                         feature.set,
                         umap.coords,
                         taxonomyDir,
                         taxonomyName = "AI_taxonomy",
                         cluster_colors = NULL,
                         metadata_names = NULL,
                         subsample=2000,
                         reorder.dendrogram = FALSE,
                         dend = NULL){

  ## Checks
  if(!"cluster" %in% colnames(meta.data)){stop("cluster must be defined in the meta.data object")}
  if(is.null(feature.set)){stop("Compute variable features and supply feature.set")}
  if(is.null(umap.coords)){stop("Compute UMAP dimensions and supply umap.coords")}
  if(!all(colnames(counts) == rownames(meta.data))){stop("Colnames of `counts` and rownames of `meta.data` do not match.")}
  if(!is.data.frame(meta.data)){stop("meta.data must be a data.frame, convert using as.data.frame(meta.data)")}
  #if(!is.matrix(counts)){stop("counts must be of type matrix.")}
	if(!is.null(dend)){
	  if(!is.element("dendrogram",class(dend))){stop("If provided, dend must be of R class dendrogram.")}
	  clusters=unique(meta.data$cluster)
	  extra_labels <- setdiff(labels(dend),clusters)
	  if(length(extra_labels)>0){stop(paste("Dendrogram has labels not included in metadata:",paste(extra_labels,collapse=", ")))}
	  extra_labels <- setdiff(clusters,labels(dend))
	  if(length(extra_labels)>0){
	    warning(paste0("Metadata include cluster labels not found in dendrogram: ",paste(extra_labels,collapse=", "),
	                   ". Cells from these clusters will be EXCLUDED from all taxonomy files."))
	  }
	}
  
  ## Ensure directory exists, if not create it
  dir.create(taxonomyDir, showWarnings = FALSE)

  ## Subsample nuclei per cluster, max
  kpClusters <- rep(TRUE,length(meta.data$cluster))
  if(!is.null(dend))
    kpClusters <- is.element(meta.data$cluster,labels(dend)) # exclude missing clusters, if any
  if((subsample > 0) & (subsample < Inf)){
      kpSub = colnames(counts)[subsampleCells(meta.data$cluster, subsample)&kpClusters]
  }else{
      kpSub = colnames(counts)[kpClusters]
  }

  ## Get the data and metadata matrices
  counts    = counts[,kpSub]
  meta.data = meta.data[kpSub,]
  umap.coords = umap.coords[kpSub,]

  ## Create log TPM matrix from counts
  tpm.matrix = scrattch.hicat::logCPM(counts)

  ## Next, generate and output the data matrices. Make sure matrix class is dgCMatrix and not dgTMatrix.  
  sample_id = colnames(tpm.matrix); meta.data$sample_id = sample_id
  gene      = rownames(tpm.matrix)
  cluster   = meta.data$cluster; names(cluster) = rownames(meta.data) ## For get_cl_medians

  ## Counts feather - Not used for shiny, but useful for saving raw data nonetheless
  print("===== Building counts.feather =====")
  counts.tibble = as_tibble(counts)
  counts.tibble = cbind(gene, counts.tibble)
  counts.tibble = as_tibble(counts.tibble)
  write_feather(counts.tibble, file.path(taxonomyDir, "counts.feather"))

  ## Data_t feather
  print("===== Building data_t.feather =====")
  data.tibble = as_tibble(tpm.matrix)
  data.tibble = cbind(gene, data.tibble)
  data.tibble = as_tibble(data.tibble)
  write_feather(data.tibble, file.path(taxonomyDir, "data_t.feather"))
  
  ## Data feather
  print("===== Building data.feather =====")
  norm.data.t = t(as.matrix(tpm.matrix))
  ## norm.data.t = Matrix::t(tpm.matrix)
  norm.data.t = as_tibble(norm.data.t)
  norm.data.t = cbind(sample_id, norm.data.t)
  norm.data.t = as_tibble(norm.data.t)
  write_feather(norm.data.t, file.path(taxonomyDir, "data.feather"))

  ## ----------
  ## Run auto_annotate, this changes sample_id to sample_name.
  print("===== Format metadata table for R shiny folder =====")
  
  # Check for duplicate columns (e.g., XXXX_label and XXXX together will crash auto_annotate)  # NEW #
  column_names = colnames(meta.data)
  column_names_revised = gsub("_label$","", column_names)
  duplicates <- names(table(column_names_revised))[table(column_names_revised)>1]
  if(length(duplicates > 0)){
    warning("Duplicate entries for", paste(duplicates, collapse=" and "),"have been DELETED. Please check output carefully!")
    meta.data <- meta.data[,setdiff(column_names, duplicates)]
  }

  # Run auto_annotate
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
  if(!is.null(dend)){
    print("...using provided dendrogram.")
    # FOR FUTURE UPDATE: should check here whether dendrogram colors match what is in meta-data.
  } else {

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
  
    ## Build the dendrogram
    invisible(capture.output({  # Avoid printing lots of numbers to the screen
      dend.result = build_dend(
        cl.dat  = medianExpr[feature.set,],
        cl.cor  = NULL,
        l.color = use.color,
        l.rank  = l.rank, 
        nboot   = 1,
        ncores  = 1)
    }))
    dend = dend.result$dend
    print("...dendrogram built.")
  }
  
  ## Output tree
  saveRDS(dend, file.path(taxonomyDir,"dend.RData"))

  ## Output tree order
  outDend = data.frame(cluster = labels(dend), order = 1:length(labels(dend)))
  write.csv(outDend, file.path(taxonomyDir, "ordered_clusters.csv"))
    
  ## Write the desc file.  
  anno_desc = create_desc(meta.data, use_label_columns = TRUE)
  # Subset the desc file to match metadata_names, if provided
  if(!is.null(metadata_names)){
    desc <- anno_desc[match(names(metadata_names), as.character(as.matrix(anno_desc[,1]))),]
    desc[,2] <- as.character(metadata_names[as.character(as.matrix(desc[,1]))])
    desc <- desc[!is.na(desc$base),]  # Remove missing values
    anno_desc <- desc
  }
  write_feather(anno_desc, file.path(taxonomyDir,"desc.feather"))

  ## Minor reformatting of metadata file, then write metadata file
  meta.data$cluster = meta.data$cluster_label; 
  colnames(meta.data)[colnames(meta.data)=="sample_name"] <- "sample_name_old" # Rename "sample_name" to avoid shiny crashing
  if(!is.element("sample_id", colnames(meta.data))){ meta.data$sample_id = meta.data$sample_name_old } ## Sanity check for sample_id
  meta.data$cluster_id <- as.numeric(factor(meta.data$cluster_label,levels=labels(dend))) # Reorder cluster ids to match dendrogram
  write_feather(meta.data, file.path(taxonomyDir,"anno.feather"))

  ## Write the UMAP coordinates.  
  print("===== Building umap/tsne feathers (precalculated) =====")
  tsne      = data.frame(sample_id = rownames(umap.coords),
                         all_x = umap.coords[,1],
                         all_y = umap.coords[,2])
  tsne      = tsne[match(meta.data$sample_id, tsne$sample_id),]
  tsne_desc = data.frame(base = "all",
                         name = "All Cells UMAP")

  write_feather(tsne, file.path(taxonomyDir,"tsne.feather"))
  write_feather(tsne_desc, file.path(taxonomyDir,"tsne_desc.feather"))

  ##
  print("===== Building count, median, sum feathers =====")
  all_clusters = unique(meta.data$cluster_label) ## cluster_id
  #all_clusters = all_clusters[order(all_clusters)]
  #allClust = paste0("cluster_", all_clusters)
  count_gt0 = matrix(0, ncol = length(all_clusters), nrow = nrow(data.tibble))
  count_gt1 = sums = medianmat = count_gt0

  ## Compute the number of genes with greater than 0,1 counts, gene sums and medians
  for (i in 1:length(all_clusters)) {
    cluster = all_clusters[i]
    cluster_samples = which(meta.data$cluster_label == cluster)
    cluster_data    = tpm.matrix[,cluster_samples]
    cluster_counts  = counts[,colnames(cluster_data)]
    count_gt0[, i]  = Matrix::rowSums(cluster_counts > 0)
    count_gt1[, i]  = Matrix::rowSums(cluster_counts > 1)
    sums[, i]       = Matrix::rowSums(cluster_counts)
    medianmat[, i]  = apply(cluster_data, 1, median)
  }

  ##
  colnames(count_gt0) = colnames(count_gt1) = colnames(sums) = colnames(medianmat) = all_clusters ## allClust
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
  write_feather(count_gt0, file.path(taxonomyDir,"count_gt0.feather"))
  write_feather(count_gt1, file.path(taxonomyDir,"count_gt1.feather"))
  write_feather(count_n,   file.path(taxonomyDir,"count_n.feather"))
  write_feather(medianmat, file.path(taxonomyDir,"medians.feather"))
  write_feather(sums,      file.path(taxonomyDir,"sums.feather"))

  ##
  print("===== Building taxonomy anndata =====")
  datReference = as.matrix(norm.data.t[,names(norm.data.t)!="sample_id"])
  annoReference = meta.data; rownames(annoReference) = meta.data$sample_id

  ##
  clustersUse = labels(dend)
  ## Define clusterInfo, which is used to convert cell types to subclass / neighborhood / class
  clusterInfo = as.data.frame(annoReference[match(clustersUse, annoReference$cluster_label),])

  ## Write taxonomy anndata
  AIT.anndata = AnnData(
    X = datReference, ## logCPM
    obs = annoReference,
    var = data.frame("gene" = colnames(datReference), 
                     "highly_variable_genes" = colnames(datReference) %in% feature.set, 
                     row.names=colnames(datReference)),
    layers = list(
      counts = Matrix::t(counts) ## Count matrix
    ),
    obsm = list(
      umap = umap.coords ## A data frame with sample_id, and 2D coordinates for umap (or comparable) representation(s)
    ),
    uns = list(
      dend        = list("standard" = file.path(taxonomyDir, "dend.RData")),  # FILE NAME with dendrogram
      filter      = list("standard" = rep(FALSE, nrow(datReference))),
      QC_markers  = list("standard"),# = file.path(taxonomyDir, "QC_markers.RData")),  # REPLACED WITH CODE BELOW
      clustersUse = clustersUse, #list("standard" = clustersUse),
      clusterInfo = clusterInfo, #list("standard" = clusterInfo),
      taxonomyName = taxonomyName,
      taxonomyDir = taxonomyDir
    )
  )
  AIT.anndata$uns$QC_markers[["standard"]]$allMarkers = allMarkers
  AIT.anndata$uns$QC_markers[["standard"]]$markers    = markers
  AIT.anndata$uns$QC_markers[["standard"]]$countsQC   = countsQC
  AIT.anndata$uns$QC_markers[["standard"]]$cpmQC      = cpmQC
  AIT.anndata$uns$QC_markers[["standard"]]$classBr    = classBr
  AIT.anndata$uns$QC_markers[["standard"]]$subclassF  = subclassF
  AIT.anndata$write_h5ad(file.path(taxonomyDir, "AI_taxonomy.h5ad"))
}

#' Add marker genes to reference dendrogram for tree mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param mode Taxonomy mode to determine which version of filtering to use.
#' @param celltypeColumn Column name correspond to the cell type names in the dendrogram (default = "cluster_label"). At least two cells per cell type in the dendrogram must be included.
#' @param subsample The number of cells to retain per cluster (default = 100)
#' @param num.markers The maximum number of markers to calculate per pairwise differential calculation per direction (default = 20)
#' @param de.param Differential expression (DE) parameters for genes and clusters used to define marker genes.  By default the values are set to the 10x nuclei defaults from scrattch.hicat, except with min.cells=2 (see notes below).
#' @param calculate.de.genes Default=TRUE. If set to false, the function will search for a file called "de.genes.rda" to load precalculated de genes.  
#' @param save.shiny.output Should standard output files be generated and saved to the directory (default=TRUE).  These are not required for tree mapping, but are required for building a patch-seq shiny instance.  This is only tested in a UNIX environment.  See notes.
#' @param mc.cores Number of cores to use for running this function to speed things up.  Default = 1.  Values>1 are only supported in an UNIX environment and require `foreach` and `doParallel` R libraries.
#' @param bs.num Number of bootstrap runs for creating the dendrogram (default of 100)
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#' @param taxonomyDir The location to save shiny output, if desired
#' @param overwriteMarkers If markers already are calculated a tree, should they be overwritten (default = TRUE)
#'
#' NOTES
#' By default VERY loose parameters are set for de_param in an effort to get extra marker genes for each node.  The defaults previously proposed for 10x nuclei are the following `de_param(low.th = 1, padj.th = 0.01, lfc.th = 1, q1.th = 0.3, q2.th = NULL, q.diff.th = 0.7, de.score.th = 100, min.cells = 2, min.genes = 5)`. See the function `de_param` in the scrattch.hicat for more details.  
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
addDendrogramMarkers = function(AIT.anndata,
                                mode = AIT.anndata$uns$mode,
                                celltypeColumn = "cluster_label",
                                subsample = 100,
                                num.markers = 20,
                                de.param=scrattch.hicat::de_param(low.th = 0.1,
                                                                  padj.th = 1,
                                                                  lfc.th = 0.01,
                                                                  q1.th = 0.1,
                                                                  q2.th = NULL,
                                                                  q.diff.th = 0.1,
                                                                  de.score.th = 1,
                                                                  min.cells = 1,
                                                                  min.genes = 1),
                                calculate.de.genes = TRUE,
                                save.shiny.output = TRUE,
                                mc.cores=1, 
                                bs.num=100, 
                                p=0.8, 
                                low.th=0.1,
                                taxonomyDir = getwd(),
                                overwriteMarkers = TRUE){

  ## We should already know this? Clean up in future.
  if(!is.element(celltypeColumn, colnames(AIT.anndata$obs))){stop(paste(celltypeColumn, "is not a column in the metadata data frame."))}

  ##
  if(mode == "standard"){ taxonomyModeDir = AIT.anndata$uns$taxonomyDir } else { taxonomyModeDir = file.path(AIT.anndata$uns$taxonomyDir, mode) }
  if(!dir.exists(taxonomyModeDir)){  stop("Taxonomy version doesn't exist, please run `buildPatchseqTaxonomy()` then retry.") }

  ## Filter
  AIT.anndata = AIT.anndata[!AIT.anndata$uns$filter[[mode]]]

  ## Subsample
  keep.samples = subsampleCells(AIT.anndata$obs[[celltypeColumn]], subsample) ##  & is.element(cluster.vector, labels(dend))

  ## Checks and data formatting
  dend = readRDS(AIT.anndata$uns$dend[[mode]])

  ## norm.data
  norm.data = Matrix::t(AIT.anndata$X[keep.samples,])
  
  ## metadata 
  metadata = AIT.anndata$obs[keep.samples,] %>% as.data.frame()

  ## celltype labels
  cluster.vector = setNames(metadata[,celltypeColumn], rownames(metadata))
  cluster.vector = factor(cluster.vector, levels=labels(dend))
  
  ## Define and collect marker genes
  print("Define some relevant variables")
  cl.df = as.data.frame(metadata[match(labels(dend), cluster.vector),])
  cl.df$cluster_label = cl.df[,celltypeColumn]
  rownames(cl.df) = 1:length(labels(dend)) 
  cl.label  = as.factor(setNames(cl.df$cluster_label, rownames(cl.df)))
  select.cl = droplevels(as.factor(setNames(match(metadata[,celltypeColumn],cl.label), metadata$sample_id)))
  
  ## CHECK IF THIS IS NEEDED
  ## We might need to relabel the dendrogram from 1 to #clusters in order
  labels(dend) = names(cl.label)[match(labels(dend), cl.label)]
  
  ## Compute markers
  print("Define marker genes and gene scores for the tree")
  
  if((sum(!is.na(get_nodes_attr(tmp, "markers")))==0)|(overwriteMarkers==FALSE)){
  
    if((!file.exists(file.path(taxonomyModeDir, "de.genes.rda"))) | calculate.de.genes){
      print("=== NOTE: This step can be very slow (several minute to many hours).")
      print("      To speed up the calculation (or if it crashes) try decreasing the value of subsample.")
      de.genes = scrattch.hicat::select_markers(norm.dat=norm.data, 
                                cl=select.cl, 
                                n.markers= num.markers, 
                                de.param = de.param, 
                                de.genes = NULL)$de.genes
      save(de.genes, file=file.path(taxonomyModeDir, "de.genes.rda"))
    } else {
      load(file.path(taxonomyModeDir, "de.genes.rda"))
    }

    ## Check number of markers for each leaf
    min.marker.gene.count = as.numeric(as.character(lapply(de.genes, function(x) x$num)))
    if(sum(min.marker.gene.count<2)>0)({
      stop("Marker genes could not be calculated for at least one node in the tree. Tree mapping will not work in this situation. We recommend loosening the de.param parameters and trying again, but warn that some cell types may not be well resolved.")
    })

    ## Note: this is not a robust way to address this issue!
    print("Calculate gene scores")
    gene.score = scrattch.hicat::get_gene_score(de.genes)
  
    print("Build the reference dendrogram")
    invisible(capture.output({  # Avoid printing lots of numbers to the screen
      reference = build_reference(cl=select.cl, 
                                  norm.dat=norm.data, 
                                  dend=dend, 
                                  de.genes=de.genes, 
                                  cl.label=cl.label, 
                                  up.gene.score=gene.score$up.gene.score, 
                                  down.gene.score=gene.score$down.gene.score, 
                                  n.markers=30)
    }))
    labels(reference$dend) = setNames(colnames(reference$cl.dat), colnames(reference$cl.dat))
    print("...marker gene calculation for reference complete")
  } else { # NOTE: THIS IS UNTESTED
    print("Build the reference dendrogram")
    reference <- list(
      dend = dend,
      cl.dat = get_cl_means(norm.dat, cl)
    )
    print("...existing markers used in reference")
  }
    
  ##
  if(sum(!is.na(get_nodes_attr(reference$dend, "original_label"))>0)){
    print("This section is needed if the starting dendrogram is from the nomenclature GitHub ")
    reference$dend = revert_dend_label(reference$dend,get_nodes_attr(reference$dend, "original_label"),"label")
  }
  
  if(save.shiny.output){
    print("Save the reference dendrogram")
    ##
    dend = reference$dend
    saveRDS(dend, file.path(taxonomyModeDir, "dend.RData"))

    ##
    print("Build membership table of reference vs. reference for use with patch-seq mapping")
    invisible(capture.output({  # Avoid printing lots of numbers to the screen
      memb.ref   = map_dend_membership(reference$dend, 
                                        reference$cl.dat, 
                                        map.dat=norm.data, 
                                        map.cells=names(select.cl),
                                        mc.cores=mc.cores, 
                                        bs.num=bs.num, 
                                        p=p, 
                                        low.th=low.th)
      map.df.ref = summarize_cl(reference$dend, 
                                memb.ref, 
                                norm.data)
    }))
    memb.ref   = memb.ref[metadata$sample_id,]
    map.df.ref = map.df.ref[metadata$sample_id,]
    save(memb.ref, map.df.ref, file=file.path(taxonomyModeDir, "membership_information_reference.rda"))
  }
  return(dend)
}
