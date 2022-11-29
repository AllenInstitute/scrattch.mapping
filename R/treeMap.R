#' Tree based mapping
#'
#' @param GEXRef A reference taxonomy object.
#' @param query.data A logCPM normalized matrix to be annotated.
#'
#' @return Tree mapping results as a data.frame.
#'
#' @export
treeMap = function(GEXRef, query.data){
    print("Tree-based mapping")
    ## Attempt Tree mapping
    tryCatch(
        expr = {
            ## Gather marker genes
            allMarkers = unique(unlist(get_dend_markers(GEXRef$dend)))
            allMarkers = intersect(allMarkers,rownames(query.data))
            ## Perform tree mapping
            membNode = rfTreeMapping(GEXRef$dend,GEXRef$datReference[allMarkers,GEXRef$kpSamp], 
                                        clReference, 
                                        query.data[allMarkers,])
            ## Gather results
            topLeaf = getTopMatch(membNode[,GEXRef$clustersUse])
            topLeaf = topLeaf[colnames(query.data),]
            ## Create results data.frame
            mappingTarget = data.frame(map.Tree=as.character(topLeaf$TopLeaf), 
                                       score.Tree=topLeaf$Value)
        },
        error = function(e){ 
            print("Error caught for Tree mapping.")
            print(e)
        },
        warning = function(w){
        },
        finally = {
            print("Tree mapping complete")
            return(mappingTarget)
        }
    )
}