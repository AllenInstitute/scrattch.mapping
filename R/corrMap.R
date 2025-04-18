#' Correlation based mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param query.data A logCPM normalized matrix to be annotated.
#' @param genes.to.use The set of genes to use for correlation calculation (default is the highly_variable_genes associated with the current mode). Can either be a character vector of gene names or a TRUE/FALSE (logical) vector of which genes to include.
#'
#' @import scrattch.hicat
#'
#' @return Correlation mapping results as a data.frame.
#'
#' @export
corrMap = function(AIT.anndata, 
                   query.data,
                   genes.to.use = NULL){
    print("Correlation-based mapping")
    ## Attempt Correlation mapping
    mappingTarget = tryCatch(
        expr = {
            ## Create the vector of genes to use.
            if(is.null(genes.to.use)){
              # If null, default to correct set of highly variable genes
              genes.to.use.vector <- AIT.anndata$var[,paste0("highly_variable_genes_",AIT.anndata$uns$mode)]
            } else if (is.logical(genes.to.use)) {
              if (length(genes.to.use)!=dim(AIT.anndata)[2]) stop("If genes.to.use is logical it must be the same length as the total number of genes in AIT.anndata.")
              genes.to.use.vector = genes.to.use
            } else if (is.character(genes.to.use)){
              genes.to.use <- intersect(genes.to.use,AIT.anndata$var_names)
              if (length(genes.to.use)==0) stop("No valid gene names provided in genes.to.use.")
              if (length(genes.to.use)<=2) stop("More than 2 valid gene names must be provided in genes.to.use to calculate correlation.")
              genes.to.use.vector <- is.element(AIT.anndata$var_names,genes.to.use)
            } else {
              stop("genes.to.use must be a character or logical vector.")
            }
          
            ##
            medianExpr = AIT.anndata$varm[[paste0("cluster_id_median_expr_",AIT.anndata$uns$mode)]]
            rownames(medianExpr) = rownames(AIT.anndata$var)
            colnames(medianExpr) = AIT.anndata$uns$clusterStatsColumns[[AIT.anndata$uns$mode]]
            ##
            mapping.genes = intersect(rownames(query.data),AIT.anndata$var_names[genes.to.use.vector])
            if(length(mapping.genes)<=2) stop("Too few intersecting variable genes with query genes to perform mapping.")
            corMapTarget  = cor_mapping_wrapper(query.data[mapping.genes,], medianExpr)
            mappingTarget = data.frame(map.Corr=as.character(corMapTarget$TopLeaf), 
                                        score.Corr=corMapTarget$Value)
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
#'
#' @import WGCNA
#' @import mfishtools
#'
#' @return Correlation mapping results
#'
#' @keywords internal
cor <- function(...) WGCNA::cor(...)
cor_mapping_wrapper <- function(input_logcpm, reference_medians){
  # tmpGenes <- intersect(rownames(input_logcpm), rownames(reference_medians))
  # varGn    <- apply(reference_medians[tmpGenes,],1,function(x) diff(range(x))/(max(x)+1))
  # varGenes <- names(sort(-varGn))[1:nGenes]
  mapping  <- corTreeMapping(input_logcpm, reference_medians)
  getTopMatch(memb.cl = mapping)
}