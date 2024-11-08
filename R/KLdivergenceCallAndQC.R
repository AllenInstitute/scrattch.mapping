
#' Compute the KL Divergence
#' 
#' From Wikipedia: "The Kullback–Leibler (KL) divergence (also called relative entropy and I-divergence) is a type of statistical distance that measures how one probability distribution P is different from a second, reference probability distribution Q. A simple interpretation of the KL divergence of P from Q is the expected excess surprise from using Q as a model when the actual distribution is P."  In scrattch.mapping this distance is required for computing categorical confidence calls used for assessing cell quality.
#'
#' @param query_probabilities A matrix of query cell mapping probabilities where rows represents cells, columns represent clusters and values rowsum to 1. 
#' @param reference_probabilities A matrix of reference cluster probabilities where rows and columns both represents clusters values represent confusion between cluster mappings/clusterings and rowsum to 1. 
#' @param select.cl An (optional) vector of cluster ids (e.g., row/col names of reference_probabilities) representing clusters for which KL divergence should be calculated using
#' @param select.cells An (optional) vector of cell ids (e.g., row names of query_probabilities) representing cells for which KL divergence should be calculated on
#'
#' @return A matrix of KL divergenes for each requested cell (row) in each requsted cluster (column)
#' 
#' @keywords external
compute_KLdiv <- function(query_probabilities, reference_probabilities, select.cl=NULL, select.cells=NULL){
  
  # variables prep and checks. THIS NEEDS BETTER CHECKS!
  library("LaplacesDemon")
  if(is.null(select.cells)){
    select.cells = rownames(query_probabilities)
  }
  if(is.null(select.cl)){
    select.cl <- colnames(query_probabilities)
  }
  query_probabilities <- query_probabilities[select.cells,select.cl]
  reference_probabilities <- reference_probabilities[select.cl, select.cl]
  
  # Get KL divergences
  KLdiv <- matrix(nrow = length(select.cells), ncol = length(select.cl))
  reference_probabilities[reference_probabilities==0] <- 0.0000001
  reference_probabilities <- reference_probabilities/rowSums(reference_probabilities)
  query_probabilities[query_probabilities==0] <- 0.0000001
  query_probabilities <- query_probabilities/rowSums(query_probabilities)
  
  for (cell in 1:length(select.cells)) {
    for (cl in 1:length(select.cl)){
      P <- query_probabilities[cell,select.cl]
      Q <- reference_probabilities[cl, select.cl]
      KLdiv[cell, cl] <- KLD(px= as.matrix(P),py =as.matrix(Q))$sum.KLD.px.py
    }
  }
  
  rownames(KLdiv) <- rownames(query_probabilities)
  colnames(KLdiv) <- as.character(select.cl)
  KLdiv
}


#' Calculate tree mapping quality call 
#' 
#' This function returns the "Good", "I1", "I2", "I3", "PoorQ" calls used to help assess Patch-seq quality.  These scores are based on a combination of KL divergence, correlation to the best matching type, and tree mapping to the top two types.
#'
#' @param AIT.anndata A reference taxonomy object on which buildPatchseqTaxonomy has already been run
#' @param query.mapping The output from taxonomy_mapping with both corr.map and tree.map set to TRUE, and with the same mode as currently set in AIT.anndata
#'
#' @return A matrix with "Tree_call" and some other useful columns, one per cell in query.mapping
#' 
#' @keywords external
tree_quality_call <- function(AIT.anndata, query.mapping){
  
  ## NEED TO ADD ALL THE TESTS HERE 
  
  # NEED TO EDIT location of "memb.ref" FROM separate file to specified location in uns
  load(file.path(AIT.anndata$uns$mode,"membership_information_reference.rda"))
  select.cl = intersect(colnames(memb.ref),unique(AIT.anndata$obs$cluster_label))
  memb.ref  = memb.ref[,select.cl]
  cls <- as.character(AIT.anndata$obs[rownames(memb.ref),"cluster_label"])
  reference_probability = NULL
  for (cl in select.cl){
    reference_probability <- rbind(reference_probability,colMeans(memb.ref[cls==cl,]))
  }
  rownames(reference_probability) = select.cl
  tree_memb <- query.mapping@detailed_results$tree[,select.cl]
  kldiv   <- compute_KLdiv(tree_memb,reference_probability, select.cl)
  results <- as.data.frame(cbind(query.mapping@annotations$score.Corr,
                                 t(apply(tree_memb,1,function(x) -sort(-x)[1:2])),
                                 t(apply(kldiv,1,function(x) sort(x)[1:2]))))
  colnames(results) <- c("Tree_first_cor","Tree_first_bt","Tree_second_bt","Tree_first_KL","Tree_second_KL")
  
  results <- results %>% 
    rownames_to_column("id") %>%
    mutate(Tree_not_finall_call = ifelse(Tree_first_cor > 0.5  & Tree_first_KL < 2, "Good", "PoorQ")) %>%
    mutate(Tree_call = case_when(Tree_not_finall_call == "Good" & Tree_first_bt >= 0.9 ~ "Core",
                                 Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                   Tree_first_bt + Tree_second_bt >= 0.7 &
                                   Tree_first_bt / Tree_second_bt >= 2 ~ "I1", 
                                 Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                   Tree_first_bt + Tree_second_bt >= 0.7 &
                                   Tree_first_bt / Tree_second_bt < 2 ~ "I2",
                                 Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                   Tree_first_bt + Tree_second_bt < 0.7 ~ "I3",
                                 Tree_not_finall_call == "PoorQ" ~ "PoorQ",
                                 TRUE ~ "Other")) %>%
    column_to_rownames("id") 
  results <- results[,c("Tree_call","Tree_first_cor","Tree_first_bt","Tree_second_bt","Tree_first_KL","Tree_second_KL")]
  
  results
}