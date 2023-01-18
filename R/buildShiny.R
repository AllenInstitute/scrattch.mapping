#' Build a Seurat object for buildReferenceFolder 
#' 
#' This function builds a Seurat object from a count matrix and metadata table for use with buildReferenceFolder, along with some other optional inputs. Unless otherwise specific log2(CPM(counts)+1) are saved in the "data" slot of the Seurat object because that is what is used for most downstream analysis.
#'
#' @param counts Count matrix with genes as rows and cells/nuclei as columns
#' @param metadata Matrix of metadata with rows corresponding to cells/nuclei, and row names corresponding to cell names. IF column names of counts are unspecified, then this matrix must be the same length and ordered in the same way as the counts matrix.  In any case, this must be a superset of counts (e.g., no missing metadata)
#' @param clusterColumn Column name correspond to the clusters to include in the clustering results
#' @param celltypeColumn Column name correspond to the cell type names to include in the dendrogram.  By default (and in most cases) this is identical to clusterColumn.
#' @param normalizebyLogCPM If TRUE (default) log(CPM+1) is used for normalization. If FALSE, `NormalizeData` Seurat function is used 
#' @param generateUMAP A boolean (default=TRUE) indicating whether UMAP should be calculated
#' @param npcs Number of pcs used for building the UMAP (default = 10)
#' @param cluster.resolution Resolution for louvain clustering algorithm
#' @param feature.set An optional user supplied gene set otherwise Seurat objects variable features are used for umap generation
#' @param celltype_colors An optional named character vector where the values correspond to colors and the names correspond to celltypes in celltypeColumn.  If this vector is incomplete, a warning is thrown and it is ignored. 
#' @param metadata_names An optional named character vector where the vector NAMES correspond to columns in the metadata matrix and the vector VALUES correspond to how these metadata should be displayed in Shiny. This is used for writing the desc.feather file later.
#' @param subsample The number of cells to retain per cluster (default is to keep all of them)
#' @param gene_names Gene names corresponding to rows in the count matrix (by default, checks the rownames)
#' @param cell_names Sample names corresponding to the columns in the count matrix (by default, checks the column names)
#'
#' @improt Seurat
#' @import scrattch.hicat
#' 
#' @return
#'
#' @export
createSeuratObjectForReferenceFolder = function(counts, 
                                                metadata,
                                                clusterColumn = "cluster",
                                                celltypeColumn = clusterColumn,
                                                normalizebyLogCPM = TRUE,
                                                generateUMAP = TRUE,
                                                npcs = 10,
                                                cluster.resolution=0.3,
                                                feature.set = NULL,
                                                celltype_colors = NULL,
                                                metadata_names = setNames(colnames(metadata),colnames(metadata)),
                                                subsample = Inf,
                                                gene_names = rownames(counts),
                                                cell_names = colnames(counts)
                                                ){
  ## Checks
  if(length(cell_names)<dim(counts)[2]){stop("Too few cell names provided.")}
  if(length(gene_names)<dim(counts)[1]){stop("Too few gene names provided.")}
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  if(length(intersect(cell_names,rownames(metadata)))<length(cell_names)){
    stop("Either missing metadata, incorrectely named, or unnamed metadata rows.")
  }
  
  # Convert clusterColumn to "cluster" if needed
  if(clusterColumn!="cluster"){
    eval(parse(text=paste0("metadata$cluster <- metadata$",clusterColumn)))
  }
  # Convert celltypeColumn to "celltype" if needed
  if(celltypeColumn!="celltype"){
    eval(parse(text=paste0("metadata$celltype <- metadata$",celltypeColumn)))
  }
  # Check if cluster colors are present and complete
  if(!is.null(celltype_colors)){
    if(length(unique(metadata$celltype))>length(interesect(metadata$celltype,names(celltype_colors)))){
      celltype_colors = NULL
      warning("Cluster color vector does not include all celltypes and will be ignored.")
    }
  }
  # Reorder and rename metadata checks
  if(length(intersect(names(metadata_names),colnames(metadata))) < 1){
    metadata_names = NULL
    warning("No value metadata names are provided and therefore metadata_names will be ignored.")
  } else {
    metadata_names <- c(setNames(c("cluster","celltype"),c("cluster","celltype")),
                        metadata_names[!is.element(metadata_names,c(clusterColumn,celltypeColumn))],
                        setNames(c("nCount_RNA","nFeature_RNA"),c("nCount_RNA","nFeature_RNA")))
  }
  # Probably need more checks!
  
  ## Format the metadata correctly
  metadata <- metadata[cell_names,]
  metadata$sample_id <- cell_names
  
  ## Subsample nuclei per cluster, if requested
  if((subsample > 0)&(subsample < Inf)){
    print("===== Subsampling nuclei =====")
    kpSub = rownames(metadata)[subsampleCells(metadata$cluster, subsample)]
    metadata <- metadata[kpSub,]
    counts   <- counts[,kpSub]
  }
  
  ## Create and normalize the Seurat object
  print("===== Building Seurat object =====")
  seurat.obj <- suppressWarnings(CreateSeuratObject(counts = counts, meta.data = metadata))
  if (normalizebyLogCPM){
    seurat.obj@assays$RNA@data <- scrattch.hicat::logCPM(counts)
  } else {
    seurat.obj <- NormalizeData(object = seurat.obj, verbose=FALSE)
  }
  
  ## Generate the UMAP if desired
  if(generateUMAP){
    print("===== Generate UMAP =====")
    feature.set <- intersect(feature.set,gene_names)
    if(length(feature.set)<=npcs){
      # If feature.set is missing or too-few features are entered, then variable genes are calculated
      seurat.obj <- FindVariableFeatures(object = seurat.obj, verbose=FALSE)
    } else{
      VariableFeatures(seurat.obj) <- feature.set
    }
    seurat.obj <- ScaleData(object = seurat.obj, verbose=FALSE)
    seurat.obj <- RunPCA(object = seurat.obj, npcs = npcs, verbose=FALSE)
    seurat.obj <- suppressWarnings(RunUMAP(object = seurat.obj, dims = 1:npcs, verbose=FALSE))
  }

  ## Cluster
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:npcs, verbose=FALSE)
  seurat.obj <- FindClusters(seurat.obj, resolution = cluster.resolution, verbose=FALSE)
  seurat.obj$cluster = seurat.obj$seurat_clusters

  ## Add the metadata names and cluster colors to misc for use in buildReferenceFolder
  print("===== Adding misc components =====")
  if(!is.null(metadata_names)) seurat.obj@misc$metadata_names <- metadata_names
  if(!is.null(celltype_colors)) seurat.obj@misc$celltype_colors <- celltype_colors
  
  ##
  return(seurat.obj)
}


#' Starting from a Seurat object this function builds the minimum files required for Shiny
#'
#' @param seurat.obj A Seurat object as specified in the notes, or as output from `createSeuratObjectForReferenceFolder`
#' @param shinyFolder The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"
#' @param subsample The number of cells to retain per cluster
#' @param feature.set An optional user supplied gene set otherwise Seurat obects variable features are used
#' @param reorder.dendrogram Should dendogram attempt to match a preset order? (Default = FALSE).  If TRUE, the dendrogram attempts to match the celltype factor order as closely as possible (if celltype is a character vector rather than a factor, this will sort clusters alphabetically, which is not ideal).
#' @param save.normalized.data Should normalized data (TRUE="data" slot in Seurat object) or raw count matrix (FALSE="counts" slot in Seurat object; default) be saved?
#'
#' @import Seurat
#' @import scrattch.hicat
#' @import scrattch.io
#' @import feather
#' @import tibble
#' @import dplyr
#' @import Matrix
#'
#' ## NOTES
#'
#' Seurat object must have the following: 
#'      Expr. data is gathered from: GetAssayData(seurat.obj, "data") and "counts"
#'      Either VariableFeatures(seurat.obj) or feature.set parameter must be supplied.
#'      seurat.obj@meta.data$cluster contains clustering either via seurat or other method, used for subsampling and dendrogram
#'      seurat.obj@meta.data$sample_id is set via colnames(seurat.obj)
#'      seurat.obj@meta.data$celltype is the labels to be shown on the dendrogram for each cluster, doesn't need to be 1-1 (but usually is)
#'      UMAP/tSNE is pulled from: FetchData(seurat.obj, vars=c("UMAP_1","UMAP_2")), simply run Seurats UMAP function or add a new dimensionSlot with the UMAP_ name.
#' 
#' @return
#'
#' @export
buildReferenceFolder = function(seurat.obj,
                                shinyFolder,
                                subsample=2000,
                                feature.set=NULL,
                                reorder.dendrogram = FALSE,
                                save.normalized.data=FALSE){

  ## Checks
  if(!"cluster" %in% colnames(seurat.obj@meta.data)){stop("cluster must be defined in the seurat object, set via seurat_clusters or other clustering")}
  if(!"celltype" %in% colnames(seurat.obj@meta.data)){print("Setting 'celltype' annotation to provided clustering"); seurat.obj$celltype = seurat.obj$cluster}
  if((length(VariableFeatures(seurat.obj)) == 0) & is.null(feature.set)){stop("Compute variable features or supply feature.set")}
	
  ## Ensure directory exists, if not create it
  dir.create(shinyFolder, showWarnings = FALSE)

  ## Subsample nuclei per cluster, max
  if((subsample > 0)&(subsample < Inf)){
      kpSub = colnames(seurat.obj)[subsampleCells(seurat.obj$cluster, subsample)]
  }else{
      kpSub = colnames(seurat.obj)
  }

  ## Get the data and metadata matrices
  seurat.obj = subset(seurat.obj, cells=kpSub)
  slot = ifelse(save.normalized.data, "data", "counts")
  tpm.matrix = GetAssayData(seurat.obj, slot, assay="RNA")
  meta.data  = seurat.obj@meta.data
  meta.data$sample_id = colnames(seurat.obj)

  ## Next, generate and output the data matrices.  Make sure matrix class is dgCMatrix and not dgTMatrix.  
  sample_id = colnames(tpm.matrix)
  gene      = rownames(tpm.matrix)

  ## Data_t feather
  print("===== Building data_t.feather =====")
  data.tibble = as_tibble(tpm.matrix)
  data.tibble = cbind(gene, data.tibble)
  data.tibble = as_tibble(data.tibble)
  write_feather(data.tibble, file.path(shinyFolder, "data_t.feather"))

  ## Counts feather - Not used for shiny, but useful for saving raw data nonetheless
  if(save.normalized.data){
    print("===== Building counts.feather (if normalized matrix saved in `data_t.feather`) =====")
    counts.tibble = as_tibble(GetAssayData(seurat.obj, "counts", assay="RNA"))
    counts.tibble = cbind(gene, counts.tibble)
    counts.tibble = as_tibble(counts.tibble)
    write_feather(counts.tibble, file.path(shinyFolder, "counts.feather"))
  }
  
  ## Data feather
  print("===== Building data.feather =====")
  norm.data.t = t(as.matrix(tpm.matrix))
  norm.data.t = as_tibble(norm.data.t)
  norm.data.t = cbind(sample_id, norm.data.t)
  norm.data.t = as_tibble(norm.data.t)
  write_feather(norm.data.t, file.path(shinyFolder, "data.feather"))

  ## ----------
  ## Run auto_annotate, this changes sample_id to sample_name.
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

  ## Adjust the celltype colors (and cluster colors if they are they same) to match celltype_colors, if available
  if(is.element("celltype_colors",names(seurat.obj@misc))){
    meta.data$celltype_color <- as.character(seurat.obj@misc$celltype_colors[meta.data$celltype_label])
    if(sum(meta.data$celltype_label!=meta.data$cluster_label)==0) meta.data$cluster_color <- meta.data$celltype_color
  }
    
  ## Define the cluster info file
  unique.meta.data = meta.data %>% distinct(cluster_id, 
                                            cluster_label, 
                                            cluster_color, 
                                            celltype_label, 
                                            celltype_color)
  rownames(unique.meta.data) = unique.meta.data$cluster_label
    
  ## Get cluster medians, proportions, and expressed genes
  print("===== Building dendrogram =====")
  medianExpr = get_cl_medians(GetAssayData(seurat.obj, "data"), seurat.obj$cluster)
  propExpr   = get_cl_prop(GetAssayData(seurat.obj, "data"), seurat.obj$cluster)
  presentGn  = apply(medianExpr,1,max) > 0

  ## Dendrogram labels
  colnames(medianExpr) = unique.meta.data[colnames(medianExpr), "celltype_label"]

  if(is.null(feature.set)){feature.set = VariableFeatures(seurat.obj)}

  ## Regenerate the dendrogram with the correct color scheme.  Note that the dendrogram was generated from the entire data set, not the subset data set.  
  use.color     = setNames(unique.meta.data$celltype_color, unique.meta.data$celltype_label)
  l.rank        = NULL
  if(reorder.dendrogram){
    l.rank = setNames(meta.data$celltype_id[match(unique.meta.data$celltype_label,meta.data$celltype_label)], unique.meta.data$celltype_label)
    l.rank = sort(l.rank)
  }
  dend.result   = build_dend(medianExpr[feature.set,],
                             cl.cor = NULL,
                             l.color = use.color,
                             nboot = 1)

  ## Output tree
  dend = dend.result$dend
  saveRDS(dend, file.path(shinyFolder,"dend.RData"))

  ## Output tree order
  outDend = data.frame(cluster = labels(dend), order = 1:length(labels(dend)))
  write.csv(outDend, file.path(shinyFolder, "ordered_clusters.csv"))
    
  ## Write the desc file.  
  anno_desc = create_desc(meta.data, use_label_columns = TRUE)
  # Subset the desc file to match metadata_names, if provided
  if(is.element("metadata_names",names(seurat.obj@misc))){
    desc <- anno_desc[match(names(seurat.obj@misc$metadata_names),as.character(as.matrix(anno_desc[,1]))),]
    desc[,2] <- as.character(seurat.obj@misc$metadata_names[as.character(as.matrix(desc[,1]))])
    desc <- desc[!is.na(desc$base),]  # Remove missing values
    anno_desc <- desc
  }
  write_feather(anno_desc, file.path(shinyFolder,"desc.feather"))

  ## Minor reformatting of metadata file, then write metadata file
  meta.data$cluster = meta.data$cluster_label; 
  colnames(meta.data)[colnames(meta.data)=="sample_name"] <- "sample_name_old" # Rename "sample_name" to avoid shiny crashing
  meta.data$cluster_id <- as.numeric(factor(meta.data$cluster_label,levels=labels(dend))) # Reorder cluster ids to match dendrogram
  write_feather(meta.data, file.path(shinyFolder,"anno.feather"))

  ## Write the UMAP coordinates.  
  print("===== Building umap/tsne feathers (precalculated) =====")
  umapCoords = FetchData(seurat.obj, vars=c("UMAP_1","UMAP_2"))

  tsne      = data.frame(sample_id = colnames(seurat.obj),
                         all_x = umapCoords[,1],
                         all_y = umapCoords[,2])
  tsne      = tsne[match(meta.data$sample_id, tsne$sample_id),]
  tsne_desc = data.frame(base = "all",
                         name = "All Cells UMAP")

  write_feather(tsne, file.path(shinyFolder,"tsne.feather"))
  write_feather(tsne_desc, file.path(shinyFolder,"tsne_desc.feather"))

  ##
  print("===== Building count,median,sum feathers =====")
  all_clusters = unique(meta.data$cluster_id)
  all_clusters = all_clusters[order(all_clusters)]
  allClust = paste0("cluster_", all_clusters)
  count_gt0 = matrix(0, ncol = length(all_clusters), nrow = nrow(data.tibble))
  count_gt1 = sums = medianmat = count_gt0

  ## Compute the number of genes with greater than 0,1 counts, gene sums and medians
  data   = GetAssayData(seurat.obj, "data", assay="RNA")
  counts = GetAssayData(seurat.obj, "counts", assay="RNA")
  for (i in 1:length(all_clusters)) {
    cluster         = all_clusters[i]
    cluster_samples = which(meta.data$cluster_id == cluster)
    cluster_data    = data[,cluster_samples]
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
#' @param celltypeColumn Column name correspond to the cell type names in the dendrogram (default = "celltype_label"). At least two cells per cell type in the dendrogram must be included.
#' @param subsample The number of cells to retain per cluster (default = 100)
#' @param num.markers The maximum number of markers to calculate per pairwise differential calculation per direction (default = 20)
#' @param de.param Differential expression (DE) parameters for genes and clusters used to define marker genes.  By default the values are set to the 10x nuclei defaults from scrattch.hicat, except with min.cells=2 (see the function `de_param` in the scrattch.hicat for more details).
#' @param calculate.de.genes Default=TRUE. If set to false, the function will search for a file called "de.genes.rda" to load precalculated de genes.  
#' @param save.shiny.output Should standard output files be generated and saved to the directory (default=TRUE).  These are not required for tree mapping, but are required for building a patch-seq shiny instance.  This is only tested in a UNIX environment.  See notes.
#' @param mc.cores Number of cores to use for running this function to speed things up.  Default = 1.  Values>1 are only supported in an UNIX environment and require and `foreach` and `doParallel` R libraries.
#' @param bs.num,p,low.th Extra variables for the `map_dend_membership` function in scrattch.hicat.  Defaults are set reasonably.
#' @param shinyFolder The location to save shiny output, if desired
#' 
#' @import feather
#' @import scrattch.hicat
#' @import MatrixGenerics
#' @import dendextend
#'
#' ## NOTES
#'
#' If save.shiny.output=TRUE, the following files will be generated:
#'   reference.rda, which includes a variable `reference` as follows:
#'       reference$cl.dat - These are the cluster means that are used for mapping comparisons
#'       reference$dend   - This is the dendrogram with marker genes attached
#'   membership_information_reference.rda, which includes two variables
#'       `memb.ref`   - matrix indicating how much confusion there is the mapping between each cell all of the nodes in the tree (including all cell types) when comparing clustering and mapping results with various subsamplings of the data
#'       `map.df.ref` - Result of tree mapping for each cell in the reference against the clustering tree, including various statistics and marker gene evidence.  This is the same output that comes from tree mapping.
#'
#' @return An updated dendrogram variable that is the same as `dend` except with marker genes added to each node.
#'
#' @export
addDendrogramMarkers = function(dend,
                                norm.data,
                                metadata,
                                celltypeColumn = "celltype_label",
                                subsample = 100,
                                num.markers = 20,
                                de.param=scrattch.hicat::de_param(low.th = 1,
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
                                bs.num=100, p=0.7, low.th=0.15,
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
  if((!file.exists(paste0(shinyFolder,"de.genes.rda")))|calculate.de.genes){
    print("=== NOTE: This step is very slow (sometimes upwards of several hours to several days depending on the size of the dendrogram).  To speed up the calculation, or if it crashes, try decreasing the value of subsample.")
    de.genes = select_markers(norm.dat=norm.data, cl=select.cl, n.markers=num.markers, de.genes = NULL)$de.genes
    save(de.genes, file=paste0(shinyFolder,"de.genes.rda"))  
  } else {
    load(paste0(shinyFolder,"de.genes.rda"))
  }
  min.marker.gene.count <- as.numeric(as.character(lapply(de.genes, function(x) x$num)))
  if(sum(min.marker.gene.count<2)>0)({
    stop("Marker genes could not be calculated for at least one node in the tree. Tree mapping will not work in this situation. We recommend loosening the de.param parameters and trying again, but warn that some cell types may not be well resolved.")
  })
  # Note: this is not a robust way to address this issue!
  print("Calculate gene scores")
  gene.score = get_gene_score(de.genes)
  
  
  print("Build the reference dendrogram")
  reference = build_reference(cl=select.cl, norm.dat=norm.data, dend=dend, de.genes=de.genes, 
                              cl.label=cl.label, up.gene.score=gene.score$up.gene.score, 
                              down.gene.score=gene.score$down.gene.score, n.markers=30)
  labels(reference$dend) <- setNames(colnames(reference$cl.dat),colnames(reference$cl.dat))
  # There is an error introduced somewhere in build_reference, which the above line fixes -- I'm not sure if this is needed and if it is, I should wrap it into build_reference
  
  
  if(sum(!is.na(get_nodes_attr(reference$dend, "original_label"))>0)){
    print("This section is needed if the starting dendrogram is from the nomenclature GitHub ")
    reference$dend <- revert_dend_label(reference$dend,get_nodes_attr(reference$dend, "original_label"),"label")
  }
  
  if(save.shiny.output){
    print("Save the reference dendrogram")
    save(reference, file=paste0(shinyFolder,"reference.rda"))
    #reference$cl.dat  # These are the cluster means that are used for mapping comparisons
    #reference$dend    # This is the dendrogram with marker genes attached
    
    print("Build membership table of reference vs. reference for use with patch-seq mapping")
    cl.dat     <- reference$cl.dat
    dend       <- reference$dend
    memb.ref   <- map_dend_membership(dend, cl.dat, map.dat=norm.data, map.cells=names(select.cl),
                                      mc.cores=mc.cores, bs.num=bs.num, p=p, low.th=low.th)
    map.df.ref <- summarize_cl(dend, memb.ref, norm.data)
    memb.ref   <- memb.ref[metadata$sample_id,]
    map.df.ref <- map.df.ref[metadata$sample_id,]
    save(memb.ref, map.df.ref, file=paste0(shinyFolder,"membership_information_reference.rda"))
  }
  
  return(reference$dend)
  
}


#' Save marker genes for patchSeqQC
#'
#' This function write a file called `QC_markers.RData` that contains all the variables required for applying the patchseq QC algorithm `pathseqtools` (which is an more flexible version of the `patchSeqQC` algorithm).  This is only used for patch-seq analysis.
# ----- Subclass calls for each cell
# ----- Broad class class calls for each cell
# ----- Distinction of neuron (e.g., mappable type) vs. non-neuron (e.g., contamination type)
#'
#' @param counts A matrix of counts for the reference data, or character string with a file location of "counts.feather".  If it appears cpm or logCPM matrix is provided, an warning will be thrown.
#' @param metadata Data frame of metadata with rows corresponding to cells/nuclei, and either row names or a column called "sample_id" corresponding to cell names. This matrix must include entries for all cells in norm.data.  Could also be a file.  Columns can be numeric, categorical, or factors and must include "subclass.column" and "class.column"
#' @param subsample The number of cells to retain per cluster (default = 100).
#' @param subclass.column Column name corresponding to the moderate-resolution cell types used for the cell types of interest (default = "subclass_label").
#' @param class.column Column name corresponding to the low-resolution cell types used for the off-target cell types (default = "class_label").
#' @param off.target.types A character vector of off-target (also known as 'contamination') cell types.  This must include at least one of the cell types found in "class.column" and/or "subclass.column" (both columns are checked)
#' @param num.markers The maximum number of markers to calculate per node per direction (default = 50)
#' @param shinyFolder = The location to save shiny output (default = current working directory).
#' 
#' @import patchseqtools
#' @import scrattch.hicat
#'
#' Nothing is returned; however an R data object called "QC_markers.RData"
#' markers, countsQC, cpmQC, classBr, subclassF, allMarkers, 
#'
#' @export
save_patchseqQC_markers = function(counts,
                                   metadata,
                                   subsample = 100,
                                   subclass.column = "subclass_label",
                                   class.column = "class_label",
                                   off.target.types = c("Glia","glia","non-neuronal","Non-neuronal"),
                                   num.markers = 50,
                                   shinyFolder = paste0(getwd(),"/")
){
  
  ## Checks and data formatting
  #counts
  if(class(counts)=="character"){
    if(file.exists(counts)){
      data_t <- try({feather(counts)})
      if(class(data_t)=="try-error"){stop("counts must be in feather format with a gene column + data matrix included.")}
      counts <- as.matrix(data_t[,colnames(data_t)!="gene"])
      if(!is.element(class(counts[,1]),c("numeric","integer"))){stop("counts doesn't have numeric entries")}
      genes <- try(data_t$gene)
      if(class(genes)=="try-error"){stop("counts must have a column called gene with gene names.")}
      rownames(counts) <- genes
    } else {
      stop(paste(counts,"is not a valid filename."))
    }
  } else {
    if(!grepl("atrix",as.character(class(counts)))){stop("counts must be some kind of matrix or a file name.")}
  }
  # Check if values are numeric and likely to be counts
  if(!is.element(class(counts[,1]),c("numeric","integer"))){stop("norm.data doesn't have numeric entries")}
  if(max(counts)<50) warning("Maximum value of counts matrix is <50.  Please check counts matrix is not log normalized.")
  if(abs(sum(counts[,1])-1000000)<1) warning("Maximum value of counts sums to 1e6.  Please check counts and not CPM was input.")
  
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
  if(length(intersect(rownames(metadata),colnames(counts)))==0){stop("Metadata sample ids (or row names) do not match sample ids in the data (column names).")}
  if(length(setdiff(colnames(counts),rownames(metadata)))>0){stop("Some samples are missing in the metadata file based on matching sample ids between counts columns and metadata rows.")}
  metadata <- metadata[colnames(counts),]
  
  #celltype labels
  if(!is.element(subclass.column,colnames(metadata))){stop(paste(subclass.column,"is not a column in the metadata data frame."))}
  if(!is.element(class.column,colnames(metadata))){stop(paste(class.column,"is not a column in the metadata data frame."))}
  metadata$subclass_label = metadata[,subclass.column]  # For compatibility with existing code.
  metadata$class_label = metadata[,class.column]  # For compatibility with existing code.
  
  # Ensure directory exists, if not create it
  dir.create(shinyFolder, showWarnings = FALSE)
  
  
  # Subsample and filter metadata and data
  kpSamp2  <- subsampleCells(metadata$subclass_label, subsample)
  goodSamp <- !is.na(metadata$class_label)  # For back-compatibility; usually not used
  kpSamp2  <- kpSamp2&goodSamp              # For back-compatibility; usually not used
  annoQC   <- metadata[kpSamp2,]
  annoQC$subclass_label = make.names(annoQC$subclass_label)
  datQC    <- as.matrix(counts[,kpSamp2])

  # Define class and subclass labels
  # --- We wrap on-target types by class but retain off-target types by subclass
  offTarget <- is.element(annoQC$class_label,off.target.types)|is.element(annoQC$subclass_label,off.target.types)
  if(sum(offTarget)==0){stop("No valid off-target classes or subclasses are provided. Please update off.target.types accordingly.")}
  
  classBr   <- annoQC$subclass_label
  classBr[!offTarget] = annoQC$class_label[!offTarget]
  classBr   <- factor(classBr)
  subclassF <- factor(annoQC$subclass_label)


  ## This is the main code here
  
  print("Define and output marker genes for each broad class and off-target subclass.") 
  # -- These are selected using some reasonable approach that could probably be improved, if needed.    
  markers    <- defineClassMarkers(datQC,subclassF,classBr,numMarkers = 50)
  allMarkers <- unique(unlist(markers))
  rownames(datQC) <- make.names(rownames(datQC))
  countsQC   <- datQC[allMarkers,]
  cpmQC      <- cpm(datQC)[allMarkers,]  # Only use of scrattch.hicat in this function

  save(markers, countsQC, cpmQC, classBr, subclassF, allMarkers, file=paste0(shinyFolder,"QC_markers.RData"))
  print("Complete!")

}