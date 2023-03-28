#' Tree based mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param query.data A logCPM normalized matrix to be annotated.
#' @param p The proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th The minimum difference in Pearson correlation to the reference cluster mean gene expression between the top-matched cluster and others. If the difference is higher than low.th, the mapping process continues; otherwise, a random branch is chosen.
#' @param bootstrap Number of bootstrapping runs to calculate the membership from (default = 100)
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
                   seed = 1)
  {
    print("Tree-based mapping")
    ## Attempt Tree mapping
    mappingTarget = tryCatch(
        expr = {
            ## Load dendrogram
            load(AIT.anndata$uns$dend)
            dend = reference$dend
            clReference  = setNames(factor(AIT.anndata$obs$cluster_label, levels=AIT.anndata$uns$clustersUse),
                                    AIT.anndata$obs_names)
            ## Gather marker genes
            allMarkers = unique(unlist(get_dend_markers(dend)))
            allMarkers = Reduce(intersect, list(allMarkers, AIT.anndata$var_names[AIT.anndata$var$common_genes]))
            ## Perform tree mapping
            invisible(capture.output({  # Avoid printing lots of numbers to the screen
              membNode = rfTreeMapping(dend, 
                                       t(AIT.anndata$X[,allMarkers]), 
                                       clReference, 
                                       query.data[allMarkers,],
                                       p=p, low.th=low.th, bootstrap=bootstrap, seed=seed)
            },type="message"))
            ## Gather results
            topLeaf = getTopMatch(membNode[,AIT.anndata$uns$clustersUse])
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
        }
    )
    return(mappingTarget)
}