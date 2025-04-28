#' Correlation based mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param query.data A logCPM normalized matrix to be annotated.
#' @param genes.to.use The set of genes to use for correlation calculation and/or Seurat integration (default is the highly_variable_genes associated with the current mode). Can be (1) a character vector of gene names, (2) a TRUE/FALSE (logical) vector of which genes to include, or (3) a column name in AIT.anndata$var corresponding to a logical vector of variable genes.
#' @param normalize.if.needed Should query.data be automatically log-normalized if it contains exceedingly large values (>30). Default = TRUE.
#'
#' @import scrattch.hicat
#'
#' @return Correlation mapping results as a data.frame.
#'
#' @export
corrMap = function(AIT.anndata, 
                   query.data,
                   genes.to.use = NULL,
                   normalize.if.needed=TRUE){
    print("Correlation-based mapping")
    ## Attempt Correlation mapping
    mappingTarget = tryCatch(
        expr = {
            ## Find and reformat inputted genes.to.use
            genes.to.use = .convert_gene_input_to_vector(AIT.anndata,genes.to.use)
            
            ## Check values of query data and if it seems to be a count matrix, log normalize
            if ((max(query.data)>25)&normalize.if.needed){
              warning("Data does not appear to be log-normalized. Trying to log-normalize.")
              query.data <- logCPM(query.data)
            }
          
            ##
            medianExpr = AIT.anndata$varm[["cluster_id_median_expr_standard"]]
            rownames(medianExpr) = rownames(AIT.anndata$var)
            colnames(medianExpr) = AIT.anndata$uns$clusterStatsColumns[["standard"]]
            medianExpr = medianExpr[,AIT.anndata$uns$clusterStatsColumns[[AIT.anndata$uns$mode]]]
            ##
            mapping.genes = intersect(rownames(query.data),AIT.anndata$var_names[genes.to.use])
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