#' Tree based mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param query.data A logCPM normalized matrix to be annotated.
#' @param p The proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th The minimum difference in Pearson correlation to the reference cluster mean gene expression between the top-matched cluster and others. If the difference is higher than low.th, the mapping process continues; otherwise, a random branch is chosen.
#' @param bootstrap Number of bootstrapping runs to calculate the membership from (default = 100)
#' @param genes.to.use The set of genes to use for tree mapping (default is all marker genes in the tree). Can be (1) a character vector of gene names, (2) a TRUE/FALSE (logical) vector of which genes to include, or (3) a column name in AIT.anndata$var corresponding to a logical vector of variable genes. If anything is provided the set of genes used is the intersection of all marker genes and the gene set here.
#' @param seed Value of the seed for reproducibility
#'
#' @return Tree mapping results as a data.frame.
#'
#' @export
treeMap = function(AIT.anndata, 
                   query.data, 
                   p = 0.8,
                   low.th = 0.1,
                   bootstrap = 100,
                   genes.to.use = rep(TRUE,dim(AIT.anndata)[2]),
                   seed = 1)
  {
    print("Tree-based mapping")
    ## Attempt Tree mapping
    mappingTarget = tryCatch(
        expr = {
            ## Load dendrogram
            dend = json_to_dend(AIT.anndata$uns$dend[[AIT.anndata$uns$mode]])
            ## Check if AIT.anndata is set up for tree mapping by seeing if marker genes are on the root node
            if (is.null(attr(dend, "markers"))){
              stop("AIT.anndata does not appear to be set up for tree mapping. Try running AIT.anndata = addDendrogramMarkers(AIT.anndata) prior to mapping.")
            }
            ##
            
            ## Find and reformat inputted genes.to.use, for use in gathering marker genes below
            genes.to.use = .convert_gene_input_to_vector(AIT.anndata,genes.to.use)
            
            ## Filter cells
            retainCells <- !AIT.anndata$uns$filter[[AIT.anndata$uns$mode]]
            
            ## Determine the cluster column, then set clReference
            hierarchy = AIT.anndata$uns$hierarchy
            hierarchy = hierarchy[order(unlist(hierarchy))]
            if(is.null(hierarchy)) stop("Hierarchy must be included in the standard AIT mode in proper format to create a mode.  Please run checkTaxonomy().")
            celltypeColumn = names(hierarchy)[length(hierarchy)][[1]]
            
            clReference  = setNames(factor(AIT.anndata$obs[,celltypeColumn], levels= AIT.anndata$uns$clusterStatsColumns[[AIT.anndata$uns$mode]]), AIT.anndata$obs_names)
            ## Gather marker genes
            commonGenes = intersect(AIT.anndata$var_names[genes.to.use],rownames(query.data))
            allMarkers = unique(unlist(get_dend_markers(dend)))
            allMarkers = Reduce(intersect, list(allMarkers, commonGenes)) 
            #old_allMarkers = Reduce(intersect, list(allMarkers, AIT.anndata$var_names[AIT.anndata$var[,paste0("highly_variable_genes_",AIT.anndata$uns$mode)]]))  # NOT used, but saved for historic reasons
            ## Perform tree mapping
            invisible(capture.output({  # Avoid printing lots of numbers to the screen
              referenceMatrix = Matrix::t(AIT.anndata$X[,allMarkers])
              referenceMatrix = as.matrix(referenceMatrix)[,retainCells]
              membNode = scrattch.mapping::rfTreeMapping(dend, 
                                       referenceMatrix, 
                                       droplevels(clReference[retainCells]), 
                                       query.data[allMarkers,],
                                       p=p, 
                                       low.th=low.th, 
                                       bootstrap=bootstrap, 
                                       seed=seed)
            },type="message"))
            ## Gather results
            topLeaf = getTopMatch(membNode[, AIT.anndata$uns$clusterStatsColumns[[AIT.anndata$uns$mode]]])
            topLeaf = topLeaf[colnames(query.data),]
            ## Create results data.frame
            mappingTarget = data.frame(map.Tree=as.character(topLeaf$TopLeaf), 
                                        score.Tree=topLeaf$Value)
            mappingTarget
        },
        error = function(e){ 
            print("Error caught for Tree mapping.")
            print(e)
            return(NULL)
        },
        finally = {
            return(list("result"=mappingTarget, "detail"=membNode))
        }
    )
}