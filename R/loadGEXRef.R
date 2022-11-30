#' Read in a reference data set in Allen taxonomy format
#'
#' @param refFolder Directory containing the Shiny taxonomy.
#' @param sample_id Field in reference taxonomy that defines the sample_id.
#' @param nGenes Number of variable genes to compute.
#' @param hGenes User supplied variable gene vector.
#' @param sub.sample Number of cells to keep per cluster.
#'
#' @return Organized reference object ready for mapping against.
#'
#' @export
loadGEXRef = function(refFolder, sample_id = "sample_id", nGenes=2000, hGenes=NULL, sub.sample = 1000){ 
    
    print("Loading reference taxonomy into memory.")

    ## Read in reference data and annotation files and format correctly
    annoReference = feather(file.path(refFolder,"anno.feather")) 
    exprReference = feather(file.path(refFolder,"data.feather"))

    ## Match meta.data to data
    annoReference = as.data.frame(annoReference[match(exprReference[[sample_id]], annoReference[[sample_id]]),])

    ## 
    datReference = as.matrix(exprReference[,names(exprReference)!=sample_id])  
    rownames(datReference) = rownames(annoReference) = annoReference[[sample_id]]
    datReference = t(datReference)

    ## Log2 cpm to normalize the data
    datReference = logCPM(datReference)

    ## Consider only genes present in both data sets
    if(!is.null(hGenes)){ hGenes = intersect(hGenes, rownames(datReference)) } else { hGenes = rownames(datReference) }

    ## Find most variable genes by beta score
    cl          = setNames(annoReference$cluster_label,colnames(datReference))
    propExpr    = get_cl_prop(datReference[hGenes,], cl)
    betaScore   = getBetaScore(propExpr)
    betaOut     = data.frame(Gene=hGenes,BetaScore=betaScore[hGenes])
    betaOut     = betaOut[order(-betaScore),]
    varFeatures = betaOut$Gene[1:nGenes]

    ## Read in the dendrogram for tree mapping and cluster order
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

    ## Create a Seurat object for the reference data set subset to a max of 1000 cells per cluster
    kpSamp        = subsampleCells(annoReference$cluster_label, sub.sample)
    referenceData = CreateSeuratObject(counts = datReference[,kpSamp], meta.data = as.data.frame(annoReference)[kpSamp,])

    ## Return named list
    return.ref = list(annoReference, exprReference, datReference, dend, clustersUse, clusterInfo, kpSamp, referenceData, varFeatures)
    names(return.ref) = c('annoReference', 'exprReference', 'datReference', 'dend', 'clustersUse', 'clusterInfo', 'kpSamp', 'referenceData', 'varFeatures')
    return(return.ref)
}