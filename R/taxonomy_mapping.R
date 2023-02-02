#' Cell type annotation and initial QC
#'
#' Perform initial mapping using three methods: Correlation-based, tree-based, and Seurat based, and will calculate some QC metrics.   
#'
#' @param GEXRef A reference taxonomy object.
#' @param query.data A logCPM normalized matrix to be annotated.
#' @param label.cols Column names of annotations to map against
#' @param corr.map Should correlation mapping be performed?
#' @param tree.map Should tree mapping be performed?
#' @param seurat.map Should seurat mapping be performed?
#' @param dims Number of PCA dimensions for Seurat mapping.
#' @param k.weight K neighborhood weight for Seurat mapping.
#'
#' @return Mapping results from all methods.
#'
#' @export
taxonomy_mapping = function(GEXRef, query.data, corr.map=TRUE, tree.map=TRUE, seurat.map=TRUE, label.cols = c("cluster_label","subclass_label", "class_label"), dims=30, k.weight=15){

    print(paste("==============================","Mapping","======================="))
    print(date())
    mappingResults=list()

    ## Sanity check on user input and taxonomy/reference annotations
    if(!all(label.cols %in% colnames(GEXRef$clusterInfo))){
      stop("Not all label.cols exists in GEXRef$clusterInfo")
    }

    ############
    ## ----- data and annotation variables -------------------------------------------------------------

    ## Ensure variable features are in common across data
    GEXRef$varFeatures = intersect(GEXRef$varFeatures, rownames(query.data))

    ############
    ## ----- Correlation mapping ------------------------------------------------------------------------
    if(corr.map == TRUE){ mappingResults[["Corr"]] = corrMap(GEXRef, query.data) }
    
    #############
    ## ----- Tree mapping -------------------------------------------------------------------------------
    if(tree.map == TRUE & !is.null(GEXRef$dend)){ mappingResults[["Tree"]] = treeMap(GEXRef, query.data) }
    
    #############
    ## ----- Seurat mapping ------------------------------------------------------------------------------
    if(seurat.map == TRUE){ mappingResults[["Seurat"]] = seuratMap(GEXRef, query.data) }

    #############
    ## Combine mapping results
    mappingAnno = Reduce(cbind, mappingResults)

    #############
    ## ---- Convert cell type mappings to subclass, neighborhood (if available), and class -------------------------------
    methods <- sort(colnames(mappingAnno)[grepl("map",colnames(mappingAnno))])
    names(methods) <- gsub("map.", "", methods)

    ## Now map back up the tree to subclass and class based on cluster labels
    for(method in names(methods)){
        convert <- GEXRef$clusterInfo[match(mappingAnno[,methods[names(methods) == method]], GEXRef$clusterInfo$cluster_label), label.cols, drop=F]
        colnames(convert) <- gsub("label", method, colnames(convert))
        mappingAnno <- cbind(mappingAnno, convert)
    }
    mappingAnno = mappingAnno[,-which(colnames(mappingAnno) %in% methods)]
    
    return(mappingAnno)
}