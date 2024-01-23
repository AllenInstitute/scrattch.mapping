#' Cell type annotation and initial QC
#'
#' Perform initial mapping using three methods: Correlation-based, tree-based, and Seurat based, and will calculate some QC metrics.   
#'
#' @param AIT.anndata A reference taxonomy object.
#' @param query.data A logCPM normalized matrix to be annotated.
#' @param label.cols Column names of annotations to map against.  Note that this only works for metadata that represent clusters or groups of clusters (e.g., subclass, supertype, neighborhood, class)
#' @param corr.map Should correlation mapping be performed?
#' @param tree.map Should tree mapping be performed?
#' @param seurat.map Should seurat mapping be performed?
#'
#' @return Mapping results from all methods.
#'
#' @export
taxonomy_mapping = function(AIT.anndata, query.data, 
                            corr.map=TRUE, tree.map=TRUE, seurat.map=TRUE, 
                            label.cols = c("cluster_label","subclass_label", "class_label")){

    print(paste("==============================","Mapping","======================="))
    print(date())
    mappingResults=list()

    ## Sanity check on user input and taxonomy/reference annotations
    if(!all(label.cols %in% colnames(AIT.anndata$uns$clusterInfo))){
      stop("Not all label.cols exists in AIT.anndata$uns$clusterInfo")
    }

    ## Modify taxonomy based on mapping mode
    if(!is.element("mode", names(AIT.anndata$uns))){
      print("Mapping in mode: `standard`. No filtering.")
    } else if(!AIT.anndata$uns$mode %in% names(AIT.anndata$uns$filter)){
      stop(print("Mapping mode, ", AIT.anndata$uns$mode, ", is not available for current Taxonomy."))
    } else{
      ## Remove off-target cell types
      AIT.anndata = AIT.anndata[!AIT.anndata$uns$filter[[AIT.anndata$uns$mode]]]
      AIT.anndata$uns$clusterInfo = AIT.anndata$uns$clusterInfo %>% filter(cluster_label %in% unique(AIT.anndata$obs$cluster_label))
      AIT.anndata$uns$clustersUse = unique(AIT.anndata$obs$cluster_label)
    }

    ############
    ## ----- data and annotation variables -------------------------------------------------------------

    ## Ensure variable features are in common across data
    AIT.anndata$var$common_genes = AIT.anndata$var$gene %in% rownames(query.data)
    AIT.anndata$var$highly_variable_genes = AIT.anndata$var$highly_variable_genes & AIT.anndata$var$common_genes

    ############
    ## ----- Correlation mapping ------------------------------------------------------------------------
    if(corr.map == TRUE){ mappingResults[["Corr"]] = corrMap(AIT.anndata, query.data) }
    
    #############
    ## ----- Tree mapping -------------------------------------------------------------------------------
    if(tree.map == TRUE & !is.null(AIT.anndata$uns$dend)){ mappingTree = treeMap(AIT.anndata, query.data); mappingResults[["Tree"]] = mappingTree[["result"]] }
    
    #############
    ## ----- Seurat mapping ------------------------------------------------------------------------------
    if(seurat.map == TRUE){ mappingResults[["Seurat"]] = seuratMap(AIT.anndata, query.data) }

    #############
    ## Combine mapping results
    mappingAnno = Reduce(cbind, mappingResults)
 
    #############
    ## ---- Convert cell type mappings to subclass, neighborhood (if available), and class -------------------------------
    methods <- sort(colnames(mappingAnno)[grepl("map", colnames(mappingAnno))])
    names(methods) <- gsub("map.", "", methods)

    ## Now map back up the tree to subclass and class based on cluster labels
    for(method in names(methods)){
        convert <- AIT.anndata$uns$clusterInfo[match(mappingAnno[,methods[names(methods) == method]], AIT.anndata$uns$clusterInfo$cluster_label), label.cols, drop=F]
        colnames(convert) <- gsub("label", method, colnames(convert))
        mappingAnno <- cbind(mappingAnno, convert)
    }
    mappingAnno = mappingAnno[,-which(colnames(mappingAnno) %in% methods)]

    ## Build mapping class object
    resultAnno <- mappingClass(annotations = mappingAnno,
                                detailed_results = list("corr" = NA, 
                                                        "tree" = mappingTree[["detail"]], 
                                                        "seurat" = NA))
    
    ## Return annotations and detailed model results
    return(resultAnno)
}