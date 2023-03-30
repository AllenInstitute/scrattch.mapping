#' Function to subsample cells
#'
#' @param cluster.names A vector of cluster names in the reference taxonomy.
#' @param subSamp Number of cells to keep per cluster.
#' @param seed Random seed used for subsampling.
#'
#' @return Boolean vector of cells to keep (TRUE) and cells to remove (FALSE)
#'
#' @keywords internal
#' 
#' @export
subsampleCells <- function(cluster.names, subSamp=25, seed=5){
  # Returns a vector of TRUE false for choosing a maximum of subsamp cells in each cluster
  # cluster.names = vector of cluster labels in factor format
  kpSamp = rep(FALSE,length(cluster.names))
  for (cli in unique(as.character(cluster.names))){
    set.seed(seed)
    seed   = seed+1
    kp     = which(cluster.names==cli)
    kpSamp[kp[sample(1:length(kp),min(length(kp),subSamp))]] = TRUE
  }
  return(kpSamp)
}

#' Get top genes by beta (binary) score
#'
#' @param data A count (or CPM or logCPM) matrix
#' @param cluster.names A vector of cluster names in the reference taxonomy.
#' @param gene.count The number of top genes to return (Default=2000)
#'
#' @return Boolean vector of cells to keep (TRUE) and cells to remove (FALSE)
#'
#' @keywords external
#' 
#' @export
top_binary_genes <- function(data, cluster.names, gene.count=2000){
  cluster.names <- setNames(as.factor(cluster.names),colnames(data))
  propExpr  <- get_cl_prop(data,cluster.names)
  betaScore <- getBetaScore(propExpr,returnScore=FALSE)
  betaScore <- sort(betaScore)
  top.genes <- names(betaScore)[1:gene.count]
  return(top.genes)
}

##################################################################################################################
## The functions below are mapping function from scrattch.hicat dev_zy branch that are required for tree mapping

# Libraries required for these functions
#library(scrattch.hicat)
#library(MatrixGenerics)
#library(randomForest)
#library(doMC)     # for parallelization in Unix environments
#library(foreach)  # for parallelization in Unix environments


#' Function for building the standard reference format, including adding marker genes to the clustering tree
#'
#' @param cl Factor vector where values are cluster ids (e.g., a numeric vector of corresponding to cell type order in the tree) and values are sample ids for cells (e.g., this vector has length = number of cells) 
#' @param norm.dat log normalized expression data
#' @param dend Input dendrogram 
#' @param de.genes output from `display_cl` function
#' @param cl.label Factor vector here values are cluster ids (e.g., a numeric vector of corresponding to cell type order in the tree) and values are dendrogram labels (e.g., this vector has length = number of clusters) 
#' @param up.gene.score Output from `get_gene_score`
#' @param down.gene.score Output from `get_gene_score`
#' @param n.markers Number of marker genes to return per comparison (default=30)
#'
#' @return A list where `dend` is the updated dendrogram with markers attached and `cl.dat` is a matrix of cluster means
#'
#' @keywords internal
build_reference <- function(cl, norm.dat, dend, de.genes, cl.label=NULL, up.gene.score=NULL, down.gene.score=NULL, n.markers=30)
{
  suppressPackageStartupMessages({
    library(randomForest)
    library(scrattch.hicat)
  })
  
  cl.dat = get_cl_means(norm.dat, cl)
  if(is.null(up.gene.score)){
    de.gene.score = get_gene_score(de.genes)
    up.gene.score = de.gene.score[[1]]
    down.gene.score = de.gene.score[[2]]
  }    
  select.genes = intersect(row.names(norm.dat), row.names(up.gene.score))
  dend = select_dend_markers(dend, norm.dat=norm.dat, cl=cl, de.genes=de.genes,
                             up.gene.score=up.gene.score[select.genes,], 
                             down.gene.score=down.gene.score[select.genes,], n.markers=n.markers)
  dend = select_pos_dend_markers(dend= dend, norm.dat = norm.dat, cl = cl)
  if(!is.null(cl.label)){
    colnames(cl.dat) = cl.label[colnames(cl.dat)]
    labels(dend) = cl.label[labels(dend)]
  }
  dend = label_dend(dend)[[1]]
  labels(dend) <- setNames(colnames(cl.dat),colnames(cl.dat)) # This line might not be needed
  
  if(sum(!is.na(get_nodes_attr(dend, "original_label"))>0)){
    print("This section is needed if the starting dendrogram includes ccn nomenclature.")
    dend <- revert_dend_label(dend,get_nodes_attr(dend, "original_label"),"label")
  }
  
  
  return(list(cl.dat=cl.dat, dend=dend))
}


#' Strip extra annotation information from dendrogram
#'
#' @param dend R dendrogram object
#' @param value Vector of values pulled from the dendrogram
#' @param attribute Which attribute should be overwritten
#'
#' @return R dendrogram object with updated attributes
#'
#' @keywords internal
revert_dend_label <- function(dend, value, attribute="label")
{
  if(attr(dend, attribute)=="")
    attr(dend, attribute) <- value[attr(dend,"original_label")]
  if (length(dend)>1) for(i in 1:length(dend))
    dend[[i]]=revert_dend_label(dend[[i]], value=value, attribute)
  return(dend)
}


#' map_dend_membership
#'
#' @param dend R dendrogram in a specific format
#' @param cl A cluster factor object to compare to a reference
#' @param cl.dat gene by cell type matrix (I think?)
#' @param map.dat normalized data of the MAPPING data set.
#' @param map.cells names of cells to map (e.g., the column names of the cell x gene matrix)
#' @param mc.cores number of cores to run the mapping on 
#' @param bs.num Number of bootstrapping runs to calculate the membership from (default = 100)
#' @param seed = random seed
#' @param ... other variables to pass to map_dend
#'
#' @import foreach
#'
#' @return membership table
#' 
#' @keywords internal
map_dend_membership <-
  function(dend,
           cl.dat,
           map.dat,
           map.cells,
           mc.cores = 10,
           bs.num = 100,
           seed = 42,
           ...)
  {
    
    # Optional libraries for UNIX parallel implementation (likely will crash in Windows)
    if (mc.cores>1) {
      suppressPackageStartupMessages({
        library(doMC)
        library(parallel)
      })
    } 
    
    if(mc.cores ==1){
      registerDoSEQ()
    }
    else{
      registerDoMC(cores=mc.cores)
    }
    mem = foreach(i = 1:bs.num, .combine = 'c') %dopar% {
      print(i)
      map_dend(dend, cl.dat, map.dat, map.cells, seed=i, ...)
    }
    memb = data.frame(cell = names(mem), cl = mem)
    memb = table(memb$cell, memb$cl)
    memb = memb / bs.num
    tmp = get_nodes_attr(dend, "label")
    tmp = tmp[tmp %in% colnames(memb)]
    memb = memb[, tmp]
    return(memb)
  }



#' map_dend
#'
#' @param dend A dendrogram in R format 
#' @param cl A cluster factor object to compare to a reference
#' @param cl.dat gene by cell type matrix (I think?)
#' @param map.dat normalized data of the MAPPING data set.
#' @param select.cells names of cells to map (e.g., the column names of the cell x gene matrix)
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#' @param default.markers What genes to include in every bootstrap run (default is none)
#' @param seed = random seed
#' 
#' @return tree mapping to the dendrogram table (cells x nodes with values as probabilities)
#' 
#' @keywords internal
map_dend <-
  function(dend,
           cl.dat,
           map.dat,
           select.cells=colnames(map.dat),
           p = 0.8,
           low.th = 0.1,
           default.markers = NULL,
           seed = 42)
  {
    final.cl = c(setNames(rep(
      attr(dend, "label"), length(select.cells)
    ), select.cells))
    if (length(dend) <= 1) {
      return(final.cl)
    }
    markers = attr(dend, "markers")
    markers = markers[names(markers) %in% row.names(map.dat)]
    cl.g = sapply(dend, labels, simplify = F)
    names(cl.g) = 1:length(cl.g)
    genes = names(markers)
    genes = union(genes, default.markers)
    mapped.cl = resolve_cl(cl.g,
                           cl.dat,
                           markers,
                           map.dat,
                           select.cells,
                           p = p,
                           low.th = low.th,
                           seed = seed+1)
    if (length(mapped.cl) > 0) {
      for (i in unique(mapped.cl)) {
        select.cells = names(mapped.cl)[mapped.cl == i]
        if (length(select.cells) > 0) {
          final.cl = c(
            final.cl,
            map_dend(
              dend[[as.integer(i)]],
              cl.dat,
              map.dat,
              select.cells,
              p = p,
              low.th = low.th,
              seed = seed+2
            )
          )
        }
      }
      return(cl = final.cl)
    }
    
  }


#' resolve_cl
#'
#' @param cl.g Cluster labels in some format
#' @param cl.med Cluster medians
#' @param markers Genes to use as markers for this function
#' @param map.dat normalized data of the MAPPING data set.
#' @param select.cells names of cells to map (e.g., the column names of the cell x gene matrix)
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#' @param seed - random seed for reproducibility
#'
#' @return mapped.cl output
#' 
#' @keywords internal
resolve_cl <-
  function(cl.g,
           cl.dat,
           markers,
           map.dat,
           select.cells,
           p = 0.8,
           low.th = 0.1,
           seed = 42)
  {
    ##
    genes = names(markers)[markers > 0]
    tmp.cl = unlist(cl.g)
    
    ###For each branch point, find the highest expression cluster.
    tmp.med = sapply(cl.g, function(g)
      rowMaxs(cl.dat[genes, g, drop = F]))
    row.names(tmp.med) = genes
    ###Make sure the genes are discriminative between all the branches.
    genes = genes[rowMaxs(tmp.med) - rowMins(tmp.med) > 1]
    
    ###Sample the markers based on the weigts.
    ##TO DO: randomforest sometimes give importance value of 0. adjust for that.
    set.seed(seed)
    seed  = seed+1
    genes = sample(genes, round(length(genes) * p), prob = markers[genes])
    
    ###Compute the correlation with the median cluster profile.
    ###add drop=F
    cl.cor = cor(as.matrix(map.dat[genes, select.cells, drop = F]), cl.dat[genes, tmp.cl, drop =
                                                                             F])
    cl.cor[is.na(cl.cor)] = 0
    ###Compute the best match in each branch.
    tmp.score = do.call("cbind", sapply(cl.g, function(x)
      rowMaxs(cl.cor[, x, drop = F]), simplify = F))
    row.names(tmp.score) = row.names(cl.cor)
    ####Determine the best match.
    best.score = setNames(rowMaxs(tmp.score), row.names(tmp.score))
    ###determine the difference from the best match.
    diff.score = best.score - tmp.score
    
    ####Give up on cells can't be discriminated,choose one branch randomly.
    unresolved.cl = row.names(tmp.score)[rowSums(diff.score < low.th) ==
                                           ncol(diff.score)]
    set.seed(seed)
    seed  = seed+1
    mapped.cl = setNames(sample(colnames(tmp.score), length(unresolved.cl), replace =
                                  T), unresolved.cl)
    
    ###Cells mapped to one or more branches.
    mapped.cells = setdiff(row.names(cl.cor), unresolved.cl)
    ###For binary branch, done already
    if (length(cl.g) == 2) {
      mapped.cl = c(mapped.cl, setNames(colnames(diff.score)[apply(diff.score[mapped.cells, , drop =
                                                                                F], 1, which.min)], mapped.cells))
      return(mapped.cl)
    }
    ##The remaining options for mapped cells
    tmp.cl = sapply(mapped.cells, function(x)
      colnames(diff.score)[which(diff.score[x,] < low.th)], simplify = F)
    ###cells with multiple options
    resolve.cells = names(tmp.cl)[sapply(tmp.cl, length) > 1]
    ###cells with only one option. Not further job.
    mapped.cells = setdiff(mapped.cells, resolve.cells)
    if (length(mapped.cells) > 0) {
      mapped.cl = c(mapped.cl, setNames(unlist(tmp.cl[mapped.cells]), mapped.cells))
    }
    ###Resolve further options.
    if (length(resolve.cells) > 0) {
      tmp.cat = sapply(tmp.cl[resolve.cells], function(x)
        paste(x, collapse = " "))
      for (cat in unique(tmp.cat)) {
        tmp.cl = unlist(strsplit(cat, " "))
        select.cells = names(tmp.cat)[tmp.cat == cat]
        mapped.cl = c(
          mapped.cl,
          resolve_cl(
            cl.g[tmp.cl],
            cl.dat,
            markers,
            map.dat,
            select.cells,
            p = p,
            low.th = low.th
          )
        )
      }
    }
    return(mapped.cl)
  }

#' Build dend (updated to specify dendextend version of "set")
#'
#' @param cl.dat Normalized data of the REFERENCE data set
#' @param cl.cor Matrix of cell x cell correlations (calculated if not provided) 
#' @param l.rank Factor of cluster order (in a specific format)
#' @param l.color Factor of clluster colors (in a specific format)
#' @param nboot Number of bootstrapping runs to calculate the membership from (default = 100)
#' @param ncores Number of cores for performing calculations
#'
#' @return dendrogram and a couple of related things
#'
#' @import dendextend
#' @import pvclust
#' 
#' @keywords internal
build_dend <- function(cl.dat, cl.cor=NULL, l.rank=NULL, l.color=NULL, nboot=100, ncores=1)
{
  if(is.null(cl.cor)){
    cl.cor = cor(cl.dat)
  }
  pvclust.result=NULL
  if(nboot > 0){
    require(pvclust)
    parallel= FALSE
    if(ncores > 1){
      parallel = as.integer(ncores)
    }
    pvclust.result <- pvclust::pvclust(cl.dat, method.dist = "cor" ,method.hclust = "average", nboot=nboot, parallel=parallel)
    dend = as.dendrogram(pvclust.result$hclust)
    dend = label_dend(dend)$dend
    dend = dend %>% pvclust_show_signif_gradient(pvclust.result, signif_type = "bp", signif_col_fun=colorRampPalette(c("white","gray","darkred","black")))
  }
  else{
    cl.hc = hclust(as.dist(1-cl.cor),method="average")      
    dend = as.dendrogram(cl.hc)
  }
  dend = dend %>% dendextend::set("labels_cex", 0.7)
  if(!is.null(l.color)){
    dend = dend %>% dendextend::set("labels_col", l.color[labels(dend)])
  }
  dend = dend %>% dendextend::set("leaves_pch", 19) %>% dendextend::set("leaves_cex", 0.5)
  if(!is.null(l.color)){
    dend = dend %>% dendextend::set("leaves_col", l.color[labels(dend)])
  }
  if(!is.null(l.rank)){
    dend =reorder_dend(dend,l.rank)
  }
  return(list(dend=dend, cl.cor=cl.cor, pvclust.result=pvclust.result))
}

#' Compute cluster medians for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with medians for each cluster
#' @export
#' 
get_cl_medians <- function(mat, cl)
{
  library(Matrix)
  library(matrixStats)
  
  cl.med <- do.call("cbind",
                    tapply(names(cl), 
                           cl, 
                           function(x){
                             matrixStats::rowMedians(as.matrix(mat[,x]))
                           }
                    )
  )
  
  rownames(cl.med) <- rownames(mat)
  
  return(cl.med)
}