#' Correlation based mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param query.data A logCPM normalized matrix to be annotated.
#'
#' @import scrattch.hicat
#'
#' @return Correlation mapping results as a data.frame.
#'
#' @export
corrMap = function(AIT.anndata, query.data){
    print("Correlation-based mapping")
    ## Attempt Correlation mapping
    mappingTarget = tryCatch(
        expr = {
            clReference  = setNames(factor(AIT.anndata$obs$cluster_label, levels=AIT.anndata$uns$clustersUse),
                                    AIT.anndata$obs_names)
            corMapTarget = map_by_cor(Matrix::t(AIT.anndata$X[,AIT.anndata$var$highly_variable_genes]), 
                                        clReference, 
                                        query.data[AIT.anndata$var_names[AIT.anndata$var$highly_variable_genes],]) 
            mappingTarget = data.frame(map.Corr=as.character(corMapTarget[[1]]$pred.cl), 
                                        score.Corr=corMapTarget[[1]]$pred.score)
            mappingTarget
        },
        error = function(e){ 
            print("Error caught for Correlation mapping.")
            print(e)
            return(NULL)
        },
        finally = {
        }
    )
    return(mappingTarget)
}

#' Correlation based cell type mapping wrapper
#'
#' @param input_logcpm Inputted log2(CPM+1) values.
#' @param reference_medians Median expression levels (log2(CPM+1)) of reference cell types.
#' @param nGenes Number of high variance genes to include.
#'
#' @import WGCNA
#' @import mfishtools
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