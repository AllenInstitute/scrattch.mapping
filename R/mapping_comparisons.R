#' Compare and plot two sets of cluster assignments
#' 
#' This is the subset of the `compare_annotate` function that does the plotting.  It is identical to the scrattch.hicat implementation, but a bit more flexible on input formats.
#' 
#' @param cl A cluster factor object to compare to a reference
#' @param ref.cl A cluster factor object for the reference clusters
#'
#' @return g A ggplot2 dot plot object for the comparison.
#'
#' @keywords internal
compare_plot <- function(cl, ref.cl){
  
  if(!is.factor(cl)) cl = factor(cl)
  if(!is.factor(ref.cl)) cl = factor(ref.cl)
  
  if(is.null(names(cl))) names(cl) <- 1:length(cl)
  if(is.null(names(ref.cl))) names(ref.cl) <- 1:length(ref.cl)
  
  scrattch.hicat::compare_plot(cl, ref.cl)
}

#' Compare and plot two sets of cluster assignments as a heatmap
#' 
#' This is a wrapper for the function `heatmap.3` in scrattch.hicat
#' 
#' @param cl A cluster factor object to compare to a reference
#' @param ref.cl A cluster factor object for the reference clusters
#' @param threshold Maximum value to show in heatmap.  Lower values will highlight off-target expression more.
#' @param cexLab Size of the label names to display on the screen.  The function will attempt to guess this if not inputted, adjustments may be needed if not all labels are shown.
#' @param ... Other inputs to the function `heatmap.3`
#'
#' @return a list with output from heatmap.3, after displaying the heatmap to the screen.
#'
#' @keywords internal
compare_heatmap <- function(cl, 
                            ref.cl,
                            threshold=0.2,
                            cexLab=NULL,
                            Rowv=NA, Colv=NA, ylab=NULL, xlab=NULL, main=NULL, 
                            margins = c(6,6), scale="none",trace="none", dendrogram="none", ...){
  
  tab  <- pmin((t(table(clust,map[[i]]))/apply(table(clust,map[[i]]),1,max)),threshold); 
  mn   <- paste("Score:",signif(sum(tab)-threshold*length(clusters),3),"(low=good)")
  if(!is.null(main)) mn <- paste0(main,": ",mn)
  if(is.null(cexLab)) cexLab <- sqrt(1/sqrt(dim(tab)[1]))+0.2
  cexRow <- cexCol <- cexLab
  heatmap.3(tab, Rowv=Rowv, Colv=Colv, ylab=ylab, xlab=xlab, main=mn, margins = margins,
            scale=scale,trace=trace,dendrogram=dendrogram,cexRow=cexRow,cexCol=cexCol, ...)
}

#' Calculate a compactness score
#' 
#' This function calculates the compactness score, defined as the the average (Pearson) correlation-based distance between each cell and the assigned group centroid (median) using the variable genes.  If a secondary group is provided (e.g., transgenic line, cortical layer, etc.), the function first sets gene expression values for each cell as expression values for the cluster median and then returns compactness per cell summarized by the secondary group.   
#' 
#' @param query.data A logCPM normalized matrix of the data
#' @param query.group A group vector to calculate compactness distance over (e.g., cluster assignments)
#' @param query.secondary An optional secondary group vector for comparison with the primary group vector (e.g., cortical layer or transgenic line)
#' @param variable.features A precomputed set of variable features.  If not provided, all genes are used.
#'
#' @return A vector of compactness scores for each cell
#'
#' @keywords internal
compactness_distance <- function(query.data,
                                 query.group,
                                 query.secondary = NULL,
                                 variable.features = rownames(query.data))
{
  # Formatting and checks
  if(length(query.group)!=dim(query.data)[2]){stop("group must have the same length as the number of cells in the data")}
  if(!is.factor(query.group))     query.group <- factor(query.group,levels=sort(unique(query.group)))
  if(is.null(names(query.group))) names(query.group) <- colnames(query.data)
  if(sum(names(query.group)==colnames(query.data))<1){
    warning("Renaming query.group to match query.data.  Confirm these variables are in the same order.")
    names(query.group) <- colnames(query.data)        
  }
  variable.features <- intersect(variable.features,rownames(query.data))
  if(length(variable.features)<=1){stop("Not enough (<2) valid variable features in the data set are provided.")}
  
  # Define group medians and correlate against them
  medianExpr <- get_cl_medians(query.data,query.group)[rownames(query.data),levels(query.group)]

  # Perform correlation of each cell to matched cluster
  if(is.null(query.secondary)){
    corData    <- rbind(query.data[variable.features,],medianExpr[variable.features,query.group])
    corResults <- apply(corData,2,function(x) cor(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)]))
    return(setNames(1-corResults,names(query.group)))
  } 
  
  # The remainder of the script is for cases with a secondary group.
  
  # Formatting and checks
  if(length(query.secondary)!=dim(query.data)[2]){stop("secondary group, if provided, must have the same length as the number of cells in the data")}
  if(!is.factor(query.secondary))     query.secondary <- factor(query.secondary,levels=sort(unique(query.secondary)))
  if(is.null(names(query.secondary))) names(query.secondary) <- colnames(query.data)
  if(sum(names(query.secondary)==colnames(query.data))<1){
    warning("Renaming query.secondary to match query.data.  Confirm these variables are in the same order.")
    names(query.secondary) <- colnames(query.data)        
  }
    
  # Calculate and return group distances of cells within secondary group 
  cor.medians <- cor(medianExpr[variable.features,])
  correlation <- rep(0,length(query.secondary))
  names(correlation) <- names(query.secondary)
  for (g in levels(query.secondary)){
    gt  <- query.secondary==g
    dat <- cor.medians[query.group[gt],query.group[gt]]
    diag(dat) <- NA
    correlation[gt] <- as.numeric(colMeans(rbind(dat,dat),na.rm = TRUE))
  }
  1-correlation

}
