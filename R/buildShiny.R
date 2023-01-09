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
#' @param feature.set An optional user supplied gene set otherwise Seurat objects variable features are used for umap generation
#' @param celltype_colors An optional named character vector where the values correspond to colors and the names correspond to celltypes in celltypeColumn.  If this vector is incomplete, a warning is thrown and it is ignored. 
#' @param metadata_names An optional named character vector where the vector NAMES correspond to columns in the metadata matrix and the vector VALUES correspond to how these metadata should be displayed in Shiny. This is used for writing the desc.feather file later.
#' @param subsample The number of cells to retain per cluster (default is to keep all of them)
#' @param gene_names Gene names corresponding to rows in the count matrix (by default, checks the rownames)
#' @param cell_names Sample names corresponding to the columns in the count matrix (by default, checks the column names)

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
                                                feature.set = NULL,
                                                celltype_colors = NULL,
                                                metadata_names = setNames(colnames(metadata),colnames(metadata)),
                                                subsample = Inf,
                                                gene_names = rownames(counts),
                                                cell_names = colnames(counts)
                                                )
{
  ## Libraries
  suppressPackageStartupMessages({
    library(Seurat)
    library(scrattch.hicat)
  })
  
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
  if(length(intersect(names(metadata_names),colnames(metadata)))<1){
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
    seurat.obj@assays$RNA@data <- logCPM(counts)
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

  ## Add the metadata names and cluster colors to misc for use in buildReferenceFolder
  print("===== Adding misc components =====")
  if(!is.null(metadata_names)) seurat.obj@misc$metadata_names <- metadata_names
  if(!is.null(celltype_colors)) seurat.obj@misc$celltype_colors <- celltype_colors
  
  seurat.obj
}


#' Starting from a Seurat object this function builds the minimum files required for Shiny
#'
#' @param seurat.obj A Seurat object as specified in the notes, or as output from `createSeuratObjectForReferenceFolder`
#' @param shinyFolder The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"
#' @param subsample The number of cells to retain per cluster
#' @param feature.set An optional user supplied gene set otherwise Seurat obects variable features are used
#' @param save.normalized.data Should normalized data (TRUE="data" slot in Seurat object) or raw count matrix (FALSE="counts" slot in Seurat object; default) be saved?
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
                                save.normalized.data=FALSE){

  ## Libraries
  suppressPackageStartupMessages({
    library(Seurat)
    library(scrattch.hicat)
    library(scrattch.io)
    library(feather)
    library(tibble)
    library(dplyr)
    library(Matrix)
  })
  
  ## Checks
  if(!"cluster" %in% colnames(seurat.obj@meta.data)){stop("cluster must be defined in the seurat object, set via seurat_clusters or other clustering")}
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
  slot = ifelse(save.normalized.data,"data","counts")
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
  ## Now let's correctly update the annotation file for use with shiny.  We'd like to correctly order things to match the tree first.    

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
  dend.result   = build_dend(medianExpr[feature.set,],
                             cl.cor = NULL,
                             l.color = use.color,
                             nboot = 1)

  ## Output tree
  dend = dend.result$dend
  saveRDS(dend, file.path(shinyFolder,"dend.RData"))

  ## Output tree order
  outDend = data.frame(cluster=labels(dend),order=1:length(labels(dend)))
  write.csv(outDend,file.path(shinyFolder,"ordered_clusters.csv"))
    
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

  ##
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

