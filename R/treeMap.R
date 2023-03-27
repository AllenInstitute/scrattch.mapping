#' Tree based mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param query.data A logCPM normalized matrix to be annotated.
#'
#' @return Tree mapping results as a data.frame.
#'
#' @export
treeMap = function(AIT.anndata, query.data){
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
                                       query.data[allMarkers,])
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