#' Cell type annotation and initial QC
#'
#' This function performs mapping using four methods (Correlation-based, tree-based, heirarchical, and Seurat-based) and will return the top one (or in some cases more) best matching reference cell type along with some associated QC metrics.   
#'
#' @param AIT.anndata A reference taxonomy object.
#' @param query.data A logCPM normalized matrix to be annotated.
#' @param label.cols Column names of annotations to map against.  Note that this only works for metadata that represent clusters or groups of clusters (e.g., subclass, supertype, neighborhood, class) and will default to whatever is included in AIT.anndata$uns$hierarchy. This is highly related to the variable called "hierarchy" in other functions.
#' @param corr.map Should correlation mapping be performed? (see methods)
#' @param tree.map Should tree mapping be performed? (see methods)
#' @param mapmycells.hierarchical.map Should mapmycells' hierarchical mapping be performed? (see methods)
#' @param mapmycells.flat.map Should mapmycells' flat mapping be performed? (see methods)
#' @param seurat.map Should seurat mapping be performed? (see methods)
#' @param genes.to.use The set of genes to use for correlation calculation and/or Seurat integration (default is the highly_variable_genes associated with the current mode). Can either be a character vector of gene names or a TRUE/FALSE (logical) vector of which genes to include.
#'
#' Mapping methods currently available in `taxonomy_mapping` include:  
#' 
#' 1. **corr.map**: This method calculates the Pearson correlation between each cell and each cluster median, and returns the cluster with the highest correlation along with the associated correlation score. Despite being a very simple method, this works quite well in some circumstances.  
#' 2. **tree.map**: Historical implementation of tree mapping used for assigning cell types to patch-seq cells in several studies of mouse visual cortex. This method requires a dendogram and iteratively walks down the tree from the root node to the leave nodes deciding the most likely cell type based on a distict set of marker genes at each node. By subsampling genes, this method provides a bootstrapping probability/confidence. *Implementation of tree mapping herein is not fully tested, so use with caution.*  
#' 3. **mapmycells.hierarchical.map**: Current version of iterative (or hierarchical) mapping used in MapMyCells, this function imports the python `cell_type_mapper` library. It requires a leveled hierarchy (e.g., cluster columns corresponding to cell type definitions at different levels of resolutions such as "cluster" and "subclass" and "class") and performs correlation-based mapping with different marker genes for each level, iterating through the levels similar to tree mapping. Like tree mapping this method provides a bootstrapping probability/confidence by subsampling genes. We find that this method works quite well in some circumstances.  
#' 4. **mapmycells.flat.map**: A single-level implementation of hierarchical mapping. Essentially it is the same as corr.map, except that it uses a prespecified set of marker genes for calculating the correlation and that it outputs bootstrapping probabilities.
#' 5. **seurat.map**: Historical implementation of Seurat mapping for assigning cell types to patch-seq cells in a study of human temporal cortex. This method performs integration and label transfer using `FindTransferAnchors` and `TransferData` functions in Seurat (v4.4.0) with a prespecified set of variable genes and a reasonable set of parameters. *We are not maintaining this method for compatibility in Seurat versions 5.0 or higher, and therefore this function will likely fail outside of the Docker environment.*  
#'
#' @return Mapping results from all methods.
#'
#' @export
taxonomy_mapping = function(AIT.anndata, 
                            query.data, 
                            corr.map=TRUE, 
                            tree.map=TRUE, 
                            mapmycells.hierarchical.map=TRUE,
                            mapmycells.flat.map=TRUE,
                            seurat.map=TRUE, 
                            label.cols = AIT.anndata$uns$hierarchy,  # NOTE THE NEW DEFAULT
                            genes.to.use = NULL,
                            mapmycells_params_list = list()){ 

  suppressWarnings({ # wrapping the whole function in suppressWarnings to avoid having this printed a zillion times: 'useNames = NA is deprecated. Instead, specify either useNames = TRUE or useNames = FALSE.'
  
    print(paste("==============================","Mapping","======================="))
    print(date())
    mappingResults=list()
    if(sum(class(label.cols)=="list")>=1) 
      label.cols = names(AIT.anndata$uns$hierarchy[order(-as.numeric(as.character(label.cols)))]) 
      # Convert from list to character for this function, from highest to lowest resolution (opposite normal order)

    ## Sanity check on user input and taxonomy/reference annotations
    if(!all(label.cols %in% colnames(AIT.anndata$uns$cluster_info))){
      stop("Not all label.cols exists in AIT.anndata$uns$cluster_info")
    }

    ## Modify taxonomy based on mapping mode
    if(!is.element("mode", names(AIT.anndata$uns))){
      print("Mapping in mode: `standard` with no subsampling. No filtering.")
    } else if(!AIT.anndata$uns$mode %in% names(AIT.anndata$uns$filter)){
      stop(print("Mapping mode, ", AIT.anndata$uns$mode, ", is not available for current Taxonomy."))
    } else{
      ## Message if subsampling is used
      if(sum(AIT.anndata$uns$filter[[AIT.anndata$uns$mode]])>0){
        print("Removing off-target cell types and/or subsampled cells.")
      }
      ## Remove off-target cell types and/or subsampled cells
      AIT.anndata = AIT.anndata[!AIT.anndata$uns$filter[[AIT.anndata$uns$mode]]]
      AIT.anndata$uns$cluster_info = AIT.anndata$uns$cluster_info %>% filter(cluster_id %in% unique(AIT.anndata$obs$cluster_id))
      #AIT.anndata$uns$clustersUse = as.character(unique(AIT.anndata$obs$cluster_id))
      AIT.anndata$uns$filter[[AIT.anndata$uns$mode]] <- rep(FALSE,sum(!AIT.anndata$uns$filter[[AIT.anndata$uns$mode]])) # New for compatibility with Seurat mapping updates
    }

    ############
    ## ----- data and annotation variables -------------------------------------------------------------
    
    ## Determine common genes and transpose if needed
    common_genes <- AIT.anndata$var$gene %in% rownames(query.data)
    if(sum(common_genes)<10){
      warning("Too few common genes found between query and reference data set.  Attempting to transpose query data.")
      query.data <- t(query.data)
      if("dgRMatrix" %in% as.character(class(taxonomy.counts)))
        query.data <- as(query.data, "dgCMatrix")
      common_genes <- AIT.anndata$var$gene %in% rownames(query.data)
      if(sum(common_genes)<10){
        stop("Too few common genes found between query and reference data set. Ensure gene names in two data sets match (e.g., both exist and are gene symbols in row or column names of matrices).")
      }
    }

    ## Ensure variable features are in common across data
    ## -- UPDATE: pull from the mode-based highly variable gene set rather than the generic one
    AIT.anndata$var$common_genes = common_genes 
    highvar.check <- paste0("highly_variable_genes_",AIT.anndata$uns$mode)
    if(!is.element(highvar.check,colnames(AIT.anndata$var))){
      print(paste("Highly variable genes for",AIT.anndata$uns$mode,"not found! Defaulting to global highly variable genes."))
      highly_variable_genes_mode = AIT.anndata$var$highly_variable_genes
    } else {
      highly_variable_genes_mode = AIT.anndata$var[,highvar.check]
    }
    AIT.anndata$var$highly_variable_genes = highly_variable_genes_mode & AIT.anndata$var$common_genes

    ############
    ## ----- Correlation mapping ------------------------------------------------------------------------
    if(corr.map == TRUE){ 
      mappingResults[["Corr"]] = corrMap(AIT.anndata, query.data, genes.to.use = genes.to.use) 
    } else { 
      mappingResults[["Corr"]] = NULL 
    }
    
    #############
    ## ----- Tree mapping -------------------------------------------------------------------------------
    if(tree.map == TRUE & !is.null(AIT.anndata$uns$dend)){ 
      suppressWarnings({
        mappingTree = treeMap(AIT.anndata, query.data); 
      })
      mappingResults[["Tree"]] = mappingTree[["result"]] 
    } else { 
      mappingTree = NULL
    }
    
    #############
    ## ----- Seurat mapping ------------------------------------------------------------------------------
    if(seurat.map == TRUE){ 
      mappingResults[["Seurat"]] = seuratMap(AIT.anndata, query.data, genes.to.use = genes.to.use) 
    } else { 
      mappingResults[["Seurat"]] = NULL 
    }

    #############
    ## ----- Hierarchical mapping ------------------------------------------------------------------------------
    if(mapmycells.hierarchical.map == TRUE){ 
      mappingHierarchical <- hierarchicalMapMyCells(AIT.anndata, query.data, mapping_params_list=mapmycells_params_list) 
      mappingResults[["Hierarchical"]] <- mappingHierarchical[["result"]]
    } else { 
      mappingHierarchical = NULL
    }

    #############
    ## ----- Flat mapping ------------------------------------------------------------------------------
    if(mapmycells.flat.map == TRUE){ 
      mappingFlat <- hierarchicalMapMyCells(AIT.anndata, query.data, mapping_params_list=mapmycells_params_list, flat_mapping=TRUE) 
      mappingResults[["Flat"]] <- mappingFlat[["result"]]
    } else { 
      mappingFlat = NULL
    }

    #############
    ## Combine mapping results
    mappingAnno = Reduce(cbind, mappingResults)
    rownames(mappingAnno) = colnames(query.data)
 
    #############
    ## ---- Convert cell type mappings for all hierarchy levels -------------------------------
    methods <- sort(colnames(mappingAnno)[grepl("map", colnames(mappingAnno))])
    names(methods) <- gsub("map.", "", methods)

    ## Now map back up the tree to subclass and class based on cluster labels
    for(method in names(methods)){
        convert <- AIT.anndata$uns$cluster_info[match(mappingAnno[,methods[names(methods) == method]], AIT.anndata$uns$cluster_info$cluster_id), label.cols, drop=F]
        convert <- convert[,setdiff(colnames(convert),"cluster_id")] # Remove cluster_id to avoid duplication
        colnames(convert) <- paste0(colnames(convert),".",method)
        mappingAnno <- cbind(mappingAnno, convert)
    }
    colnames(mappingAnno) <- gsub("map.","cluster_id.",colnames(mappingAnno))

    ## Build mapping class object
    resultAnno <- mappingClass(annotations = mappingAnno,
                               detailed_results = list("corr" = NA, 
                                                       "tree" = mappingTree[["detail"]], 
                                                       "seurat" = NA,
                                                       "hierarchical" = mappingHierarchical[["detail"]],
                                                       "flat" = mappingFlat[["detail"]]))
    
    ## Return annotations and detailed model results
    return(resultAnno)
    
  }) # End suppressWarnings
}