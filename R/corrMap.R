#' Correlation based mapping
#'
#' @param GEXRef A reference taxonomy object.
#' @param query.data A logCPM normalized matrix to be annotated.
#'
#' @return Correlation mapping results as a data.frame.
#'
#' @export
corrMap = function(GEXRef, query.data){
    print("Correlation-based mapping")
    ## Attempt Tree mapping
    tryCatch(
        expr = {
            clReference  = setNames(factor(GEXRef$annoReference$cluster_label,levels=GEXRef$clustersUse),
                                    colnames(GEXRef$datReference))[GEXRef$kpSamp]
            corMapTarget = map_by_cor(GEXRef$datReference[GEXRef$varFeatures,GEXRef$kpSamp], 
                                        clReference, 
                                        query.data[GEXRef$varFeatures,]) 
            mappingTarget = data.frame(map.Corr=as.character(corMapTarget[[1]]$pred.cl), 
                                    score.Corr=corMapTarget[[1]]$pred.score)
        },
        error = function(e){ 
            print("Error caught for Correlation mapping.")
            print(e)
        },
        warning = function(w){
        },
        finally = {
            print("Correlation mapping complete")
            return(mappingTarget)
        }
    )
}

#' Correlation based cell type mapping wrapper
#'
#' @param input_logcpm
#' @param reference_medians
#' @param nGenes
#'
#' @return Correlation mapping results
#'
#' @keywords internal
cor <- function(...) WGCNA::cor(...)
cor_mapping_wrapper <- function(input_logcpm, reference_medians, nGenes = 1000){
  tmpGenes <- intersect(rownames(input_logcpm), rownames(reference_medians))
  varGn    <- apply(reference_medians[tmpGenes,],1,function(x) diff(range(x))/(max(x)+1))
  varGenes <- names(sort(-varGn))[1:nGenes]
  mapping  <- corTreeMapping(input_logcpm[varGenes,], reference_medians[varGenes,])
  getTopMatch(memb.cl = mapping)
}