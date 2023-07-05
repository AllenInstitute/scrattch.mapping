#' Read in a reference data set in Allen taxonomy format
#'
#' @param taxonomyDir Directory containing the Shiny taxonomy -OR- a direct h5ad file name.
#' @param anndata_file File name of the anndata object to be loaded.
#' @param sample_id Field in reference taxonomy that defines the sample_id.
#' @param hGenes User supplied variable gene vector.  If not provided, then all genes are used.
#' @param gene_id Field in counts.feather that defines the gene_id.
#' @param patchseq Should the patchseq variant of this taxonomy be loaded.
#' @param force Force rebuild the anndata object for the taxonomy.
#'
#' @return Organized reference object ready for mapping against.
#'
#' @export
loadTaxonomy = function(taxonomyDir, 
                        anndata_file = "AI_taxonomy.h5ad",
                        sample_id = "sample_id", 
                        hGenes=NULL, 
                        gene_id = "gene",
                        force=FALSE){

  ## Load from directory name input 
  if(file.exists(file.path(taxonomyDir, anndata_file)) & force == FALSE){
    print("Loading reference taxonomy into memory from .h5ad")
    ## Load taxonomy directly!
    AIT.anndata = read_h5ad(file.path(taxonomyDir, "AI_taxonomy.h5ad"))
    ## Ensure anndata is in scrattch.mapping format
    ## ..
  } else if(all(file.exists(c(file.path(taxonomyDir,"anno.feather"), 
                              file.path(taxonomyDir,"data.feather"), 
                              file.path(taxonomyDir,"counts.feather"), 
                              file.path(taxonomyDir,"tsne.feather"))))){
    ##
    print("Loading reference taxonomy into memory from .feather")
    
    ## Read in reference data and annotation files and format correctly
    annoReference   = feather(file.path(taxonomyDir,"anno.feather")) 
    exprReference   = feather(file.path(taxonomyDir,"data.feather"))
    
    ## Convert log2CPM-normalized data into a matrix
    datReference = as.matrix(exprReference[,names(exprReference)!=sample_id])

    ## Read in reference count data if available
    if(file.exists(file.path(taxonomyDir,"counts.feather"))){
      print("Loading counts matrix")
      countsReference = feather(file.path(taxonomyDir,"counts.feather"))  # Note that this is transposed relative to data.feather
      ## Convert count data into a matrix
      countsReference = as.matrix(countsReference[,names(countsReference)!=gene_id])
      rownames(countsReference) = colnames(datReference)
      countsReference = Matrix::t(countsReference)
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
    
    ## Read in cluster info
    clustersUse = unique(annoReference$cluster_label)
    clusterInfo = as.data.frame(annoReference) ## No dendrogram ordering so just convert to df
 
    ## Read in the umap
    if(!file.exists(file.path(taxonomyDir,"tsne.feather"))){
      umap.coords = NULL
    }else{
      tryCatch(
        expr = {
          umap.coords = as.data.frame(read_feather(file.path(taxonomyDir,"tsne.feather")))
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
        dend        = list("standard" = file.path(taxonomyDir, "dend.RData")),  # FILE NAME with dendrogram
        filter      = list("standard" = rep(FALSE, nrow(datReference))),
        QC_markers  = list("standard"),# = file.path(taxonomyDir, "QC_markers.RData")),  # REPLACED WITH CODE BELOW
        clustersUse = clustersUse,
        clusterInfo = clusterInfo,
        taxonomyName = "",
        taxonomyDir = taxonomyDir
      )
    )
    AIT.anndata$write_h5ad(file.path(taxonomyDir, "AI_taxonomy.h5ad")) ## Save the anndata taxonomy so the next person doesn't have to build it :).
  }else{
    stop("Required files to load Allen Institute taxonomy are missing.")
  }

  ## Set scrattch.mapping to default standard mapping mode
  AIT.anndata$uns$mode = "standard"

  ## Return
  return(AIT.anndata)
}
