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
    print("=== NOTE: This step can be very slow (several minute to many hours).")
    print("      To speed up the calculation (or if it crashes) try decreasing the value of subsample.")
    de.genes = select_markers(norm.dat=norm.data, cl=select.cl, n.markers=num.markers, 
                              de.param = de.param, de.genes = NULL)$de.genes
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
    save(reference, file=paste0(shinyFolder,"reference.rda")) # Redundant
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
    save(memb.ref, map.df.ref, file=paste0(shinyFolder,"membership_information_reference.rda"))
  }
  
  return(reference$dend)
  
}