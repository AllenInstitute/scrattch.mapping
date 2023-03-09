#' Read in a reference data set in Allen taxonomy format
#'
#' @param refFolder Directory containing the Shiny taxonomy.
#' @param sample_id Field in reference taxonomy that defines the sample_id.
#' @param hGenes User supplied variable gene vector.  If not provided, then all genes are used.
#' @param gene_id Field in counts.feather that defines the gene_id.
#'
#' @return Organized reference object ready for mapping against.
#'
#' @export
loadTaxonomy = function(refFolder, 
                        sample_id = "sample_id", 
                        hGenes=NULL, 
                        gene_id = "gene"){ 
  
  print("Loading reference taxonomy into memory.")
  
  ## Read in reference data and annotation files and format correctly
  annoReference   = feather(file.path(refFolder,"anno.feather")) 
  exprReference   = feather(file.path(refFolder,"data.feather"))
  
  ## Convert log2CPM-normalized data into a matrix
  datReference = as.matrix(exprReference[,names(exprReference)!=sample_id])

  ## Read in reference count data if available
  if(file.exists(file.path(refFolder,"counts.feather"))){
    countsReference = feather(file.path(refFolder,"counts.feather"))  # Note that this is transposed relative to data.feather
    ## Convert count data into a matrix
    countsReference = as.matrix(countsReference[,names(countsReference)!=gene_id])
    rownames(countsReference) = colnames(datReference)
    countsReference = t(countsReference)
  }else{
    countsReference = NULL
  }
  
  ## Match meta.data to data
  annoReference = as.data.frame(annoReference[match(exprReference[[sample_id]], annoReference[[sample_id]]),])
  rownames(annoReference) = rownames(datReference) = annoReference[[sample_id]]
  
  ## Consider only genes present in both data sets
  if(!is.null(hGenes)){ 
    varFeatures = intersect(hGenes, colnames(datReference)) 
  } else { 
    varFeatures = colnames(datReference) 
  }
  
  ## Read in the dendrogram for tree mapping and cluster order
  # NOTE: This section may need to be updated to read "dend.RData" once redundancies are sorted out
  # NOTE: Technically this section isn't used at all right now.  We need to sort our conversion between R and python file formats
  if(!file.exists(file.path(refFolder,"reference.rda"))){
    dend = NULL
    clustersUse = unique(annoReference$cluster_label)
    clusterInfo = as.data.frame(annoReference) ## No dendrogram ordering so just convert to df
  }else{
    tryCatch(
      expr = {
        load(file.path(refFolder,"reference.rda"))
        dend = reference$dend
        clustersUse = labels(dend)
        ## Define clusterInfo, which is used to convert cell types to subclass / neighborhood / class
        clusterInfo = as.data.frame(annoReference[match(clustersUse, annoReference$cluster_label),])
      },
      error = function(e){ 
        print("Error caught for dendogram loading.")
        dend = NULL
        clustersUse = unique(annoReference$cluster_label)
        clusterInfo = as.data.frame(annoReference) ## No dendrogram ordering so just convert to df
      },
      warning = function(w){
      },
      finally = {
        print("Done loading reference")
      }
    )
  }
  
  ## Read in the umap
  if(!file.exists(file.path(refFolder,"tsne.feather"))){
    umap.coords = NULL
  }else{
    tryCatch(
      expr = {
        umap.coords = as.data.frame(read_feather(file.path(refFolder,"tsne.feather")))
        rownames(umap.coords) = umap.coords[,sample_id]
        umap.coords = umap.coords[rownames(annoReference),]
      },
      error = function(e){ 
        print("Error caught for umap loading. Setting all UMAP values to 0.")
        l = length(rownames(annoReference))
        umap.coords = data.frame(sample_id=rownames(annoReference),all_x=rep(0,l),all_y=rep(0,l))
        rownames(umap.coords) = umap.coords[,sample_id]
      },
      warning = function(w){
      },
      finally = {
        print("Done loading UMAP")
      }
    )
  }
  # This next bit should NOT be needed, but the anndata crashes without it
  if(is_tibble(umap.coords)){
    umap.coords = as.data.frame(umap.coords)
    print(class(umap.coords))
    rownames(umap.coords) = umap.coords[,sample_id]
  }
  
  ## Build reference object
  AIT.anndata = AnnData(
    X = datReference, ## logCPM
    obs = annoReference,
    var = data.frame("gene" = colnames(datReference), 
                     "highly_variable_genes" = colnames(datReference) %in% varFeatures, 
                     row.names=colnames(datReference)),
    layers = list(
      counts = countsReference # Counts. We may want to keep genes as rows and not transpose this
    ),
    obsm = list(
      umap = umap.coords # A data frame with sample_id, and 2D coordinates for umap (or comparable) representation(s)
    ),
    uns = list(
      dend        = file.path(refFolder,"reference.rda"),  # FILE NAME with dendrogram
      QC_markers  = file.path(refFolder,"QC_markers.RData"), # FILE NAME with variables for patchseqQC
      clustersUse = clustersUse,
      clusterInfo = clusterInfo
    )
  )
  
  ## Return taxonomy anndata
  return(AIT.anndata)
}
