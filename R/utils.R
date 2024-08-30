#' Map samples to a training dataset by correlation
#' 
#' @param train.dat Training data matrix, usually log-transformed CPM
#' @param train.cl Training cluster factor object
#' @param test.dat Data for cells to map to the training set. Should have the same genes as train.dat.
#' @param method Which statistic to compare. "median" or "mean". Default is "median".
#' 
#' @return a list object containing two objects:
#' \itemize{
#' \item pred.df: a data.frame with two columns, pred.cl and pred.score with the predicted cluster and correlation scores.
#' \item cor.matrix: a matrix object with correlation scores for each cluster.
#' }
#' 
#' @export
#' 
map_by_cor <- function(cl.meds, 
                       test.dat,
                       method = "median") {
  
  method <- match.arg(arg = method, 
                      choices = c("mean","median"))
  
  # Get medians or means for each cluster
  if(method == "median"){
    cl.meds <- tapply(names(train.cl), 
                      train.cl, 
                      function(x) {
                        train.mat <- train.dat[, x, drop = F]
                        train.mat <- as.matrix(train.mat)
                        matrixStats::rowMedians(train.mat)
                      }
    )
    
    cl.dat <- do.call("cbind", cl.meds)
  } else {
    cl.dat <- get_cl_means(train.dat, train.cl)
  }
  row.names(cl.dat) <- row.names(train.dat)
  
  # Perform correlations
  if(!is.matrix(test.dat) & nrow(test.dat)*ncol(test.dat) > 1e8){
    test.cl.cor <- qlcMatrix::corSparse(test.dat, cl.dat)
    colnames(test.cl.cor) = colnames(cl.dat)
    row.names(test.cl.cor) = colnames(test.dat)
  } else{
    test.cl.cor <- cor(as.matrix(test.dat), cl.dat)
  }
  test.cl.cor[is.na(test.cl.cor)] <- 0
  
  # Find maximum correlation
  max.cl.cor <- apply(test.cl.cor, 1, which.max)
  pred.cl <- colnames(test.cl.cor)[max.cl.cor]
  pred.cl <- setNames(pred.cl, row.names(test.cl.cor))
  
  # Get maximum correlation values
  pred.score <- apply(test.cl.cor, 1, max)
  
  # Convert to factor if train.cl was a factor and match levels.
  if(is.factor(train.cl)){
    pred.cl <- setNames(factor(pred.cl, levels = levels(train.cl)), names(pred.cl))
  }
  
  # Output results
  pred.df <- data.frame(pred.cl = pred.cl,
                        pred.score = pred.score)
  
  out_list <- list(pred.df = pred.df,
                   cor.matrix = test.cl.cor)
  
  return(out_list)    
}

#' Function to set mapping mode
#'
#' @param AIT.anndata A vector of cluster names in the reference taxonomy.
#' @param mode Number of cells to keep per cluster.
#'
#' @return AIT anndata with mode set for mapping
#' 
#' @export
mappingMode <- function(AIT.anndata, mode){
  if(!mode %in% names(AIT.anndata$uns$filter)){ stop(paste0(mode, " is invalid. Choose from: ", names(AIT.anndata$uns$filter))) }
  AIT.anndata$uns$mode = mode
  ## Backwards compatibility for new QC markers storage
  # if(is.null(AIT.anndata$uns$QC_markers[[mode]]$allMarkers)){
  #   QC_marker_file = file.path(AIT.anndata$uns$taxonomyDir, mode, "QC_markers.rda")
  #   if(file.exists(QC_marker_file)){
  #     ## 
  #     print("Converting taxonomy .h5ad to new format and saving. This should only happen once per taxonomy.")
  #     ##
  #     load(QC_marker_file)
  #     AIT.anndata$uns$QC_markers[[mode]]$allMarkers = allMarkers
  #     AIT.anndata$uns$QC_markers[[mode]]$markers    = markers
  #     AIT.anndata$uns$QC_markers[[mode]]$countsQC   = countsQC
  #     AIT.anndata$uns$QC_markers[[mode]]$cpmQC      = cpmQC
  #     AIT.anndata$uns$QC_markers[[mode]]$classBr    = classBr
  #     AIT.anndata$uns$QC_markers[[mode]]$subclassF  = subclassF
  #     AIT.anndata$uns$QC_markers[[mode]]$qc_samples = colnames(countsQC) # since colnames are lost
  #     AIT.anndata$uns$QC_markers[[mode]]$qc_genes   = rownames(countsQC) # since rownames are lost
  #     AIT.anndata$write_h5ad(file.path(AIT.anndata$uns$taxonomyDir, "AI_taxonomy.h5ad"))
  #   }else{
  #     stop("Could not find QC marker files required for taxonomy mode. Please run buildPatchseqTaxonomy and try again.")
  #   }
  # }
  return(AIT.anndata)
}

#' Function to subsample cells
#'
#' @param cluster.names A vector of cluster names in the reference taxonomy.
#' @param subSamp Number of cells to keep per cluster.
#' @param seed Random seed used for subsampling.
#'
#' @return Boolean vector of cells to keep (TRUE) and cells to remove (FALSE)
#' 
#' @keywords external
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
#' @export
top_binary_genes <- function(data, cluster.names, gene.count=2000){
  cluster.names <- setNames(as.factor(cluster.names),colnames(data))
  propExpr  <- get_cl_prop(data,cluster.names)
  betaScore <- getBetaScore(propExpr,returnScore=FALSE)
  betaScore <- sort(betaScore)
  top.genes <- names(betaScore)[1:gene.count]
  return(top.genes)
}

#' Tree-based mapping
#'
#' Returns the mapping membership of each cell to each node and leaf using a
#'   tree-based method.  This is a wrapper function for map_dend.  Includes
#'   Minor adjustments from the function of the same name in `mfishtools`.
#'
#' @param dend dendrogram for mapping
#' @param refDat normalized data of the REFERENCE data set
#' @param clustersF factor indicating which cluster each cell type is actually assigned to
#'   in the reference data set
#' @param mapDat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#'   to map to. otherwise, a random branch is chosen.
#' @param bootstrap Number of bootstrapping runs to calculate the membership from (default = 100)
#' @param seed added for reproducibility
#'
#' @return a matrix of confidence scores (from 0 to 100) with rows as cells and columns
#'   as tree node/leafs.  Values indicate the fraction of permutations in which the cell
#'   mapped to that node/leaf using the subset of cells/genes in map_dend
#'
#' @keywords internal
rfTreeMapping <- function (dend, refDat, clustersF, mapDat = refDat, p = 0.7, 
                           low.th = 0.15, bootstrap = 100, seed = 1) 
{
  genes <- intersect(rownames(refDat), rownames(mapDat))
  refDat <- as.matrix(refDat)[genes, ]
  mapDat <- as.matrix(mapDat)[genes, ]
  pseq.cells <- colnames(mapDat)
  
  pseq.mem <- sapply(1:bootstrap, function(i) {
    j <- i
    if (i %% 25 == 0) print(paste0("Bootstrap: ", i))
    go <- TRUE
    while (go) {
      j <- j + 1000
      set.seed(j + seed)
      tmp <- try(mfishtools::map_dend(dend, clustersF, refDat, mapDat, pseq.cells, 
                                      p = p, low.th = low.th), silent=FALSE)
      if(length(tmp) > 1 | j > 20000) ## Give it 20 tries then fail out.
        go <- FALSE
    }
    tmp
  }, simplify = F)
  memb <- unlist(pseq.mem)
  memb <- data.frame(cell = names(memb), cl = memb)
  memb$cl <- factor(memb$cl, levels = get_nodes_attr(dend, "label"))
  memb <- table(memb$cell, memb$cl)
  memb <- memb/bootstrap
  return(memb)
}

#' Tree-based mapping
#'
#' Returns the mapping membership of each cell to each node and leaf using a
#'   tree-based method.  This is a wrapper function for map_dend.
#'
#' @param dend dendrogram for mapping
#' @param cl factor indicating which cluster each cell type is actually assigned to
#'   in the reference data set
#' @param dat normalized data of the REFERENCE data set
#' @param map.dat normalized data of the MAPPING data set.  Default is to map the
#'   data onto itself.
#' @param p proportion of marker genes to include in each iteration of the mapping
#'   algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide
#'   on which branch to map to. otherwise, a random branch is chosen.
#' @param default.markers not used
#'
#' @return a matrix of confidence scores (from 0 to 100) with rows as cells and columns
#'   as tree node/leafs.  Values indicate the fraction of permutations in which the cell
#'   mapped to that node/leaf using the subset of cells/genes in map_dend
#'
#' @export
map_dend <- function(dend,
                     cl,
                     dat,
                     map.dat,
                     select.cells,
                     p = 0.8,
                     low.th = 0.2,
                     default.markers = NULL) {
  final.cl <- c(setNames(rep(attr(dend2, "label"), length(select.cells)), select.cells))
  if (length(dend) <= 1) {
    return(final.cl)
  }
  markers <- attr(dend, "markers")
  markers <- markers[names(markers) %in% row.names(map.dat)]
  cl.g <- sapply(dend, labels, simplify = F)
  names(cl.g) <- 1:length(cl.g)
  select.cl <- cl[cl %in% unlist(cl.g)]
  ### Sampling the cells from the reference cluster
  cells <- unlist(tapply(names(select.cl), select.cl, function(x) sample(x, round(length(x) * p))))
  genes <- names(markers)
  genes <- union(genes, default.markers)
  ### Compute reference cluster median based on
  ### subsampled cells
  cl.med <- do.call("cbind", tapply(
    cells, droplevels(cl[cells]),
    function(x) rowMedians(dat[genes, x, drop = F])
  ))
  row.names(cl.med) <- genes
  ### determine which branch to take.
  mapped.cl <- resolve_cl(cl.g, cl.med, markers, dat,
    map.dat, select.cells,
    p = p, low.th = low.th
  )
  if (length(mapped.cl) > 0) {
    for (i in unique(mapped.cl)) {
      select.cells <- names(mapped.cl)[mapped.cl == i]
      if (length(select.cells) > 0) {
        final.cl <- c(final.cl, map_dend(dend[[as.integer(i)]],
          cl, dat, map.dat, select.cells,
          p = p, low.th = low.th
        ))
      }
    }
  }
  return(cl = final.cl)
}

#' Tree-based mapping (internal)
#'
#' Returns the mapped cluster call of each cell to each leaf. This function is called by map_dend
#'
#' @param cl.g all clusters
#' @param cl.med cluster medians
#' @param markers gene markers
#' @param dat normalized data of the REFERENCE data set
#' @param map.dat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param select.cells which cells to use?
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#'   to map to. otherwise, a random branch is chosen.
#'
#' @return a vector of the mapped cluster
#'
#' @export
resolve_cl <- function(cl.g,
                       cl.med,
                       markers,
                       dat,
                       map.dat,
                       select.cells,
                       p = 0.7,
                       low.th = 0.2) {
  library(matrixStats)
  ##
  genes <- names(markers)[markers > 0]
  tmp.cl <- unlist(cl.g)

  ### For each branch point, find the highest
  ### expression cluster.
  tmp.med <- sapply(cl.g, function(g) rowMaxs(cl.med[genes, g, drop = F]))
  row.names(tmp.med) <- genes
  ### Make sure the genes are discriminative between
  ### all the branches.
  genes <- genes[rowMaxs(tmp.med) - rowMins(tmp.med) > 1]

  ### Sample the markers based on the weigts. TO DO:
  ### randomforest sometimes give importance value of
  ### 0. adjust for that.
  genes <- sample(genes, round(length(genes) * p), prob = markers[genes])

  ### Compute the correlation with the median cluster
  ### profile. add drop=F
  cl.cor <- WGCNA::cor(map.dat[genes, select.cells, drop = F], cl.med[genes, tmp.cl, drop = F])
  cl.cor[is.na(cl.cor)] <- 0
  ### Compute the best match in each branch.
  tmp.score <- do.call("cbind", sapply(cl.g, function(x) rowMaxs(cl.cor[,
      x,
      drop = F
    ]), simplify = F))
  row.names(tmp.score) <- row.names(cl.cor)
  #### Determine the best match.
  best.score <- setNames(rowMaxs(tmp.score), row.names(tmp.score))
  ### determine the difference from the best match.
  diff.score <- best.score - tmp.score

  #### Give up on cells can't be discriminated,choose
  #### one branch randomly.
  unresolved.cl <- row.names(tmp.score)[rowSums(diff.score < low.th) == ncol(diff.score)]
  mapped.cl <- setNames(sample(colnames(tmp.score), length(unresolved.cl), replace = T), unresolved.cl)

  ### Cells mapped to one or more branches.
  mapped.cells <- setdiff(row.names(cl.cor), unresolved.cl)
  ### For binary branch, done already
  if (length(cl.g) == 2) {
    mapped.cl <- c(mapped.cl, setNames(colnames(diff.score)[apply(diff.score[mapped.cells,
      ,
      drop = F
    ], 1, which.min)], mapped.cells))
    return(mapped.cl)
  }
  ## The remaining options for mapped cells
  tmp.cl <- sapply(mapped.cells, function(x) colnames(diff.score)[which(diff.score[x, ] < low.th)], simplify = F)
  ### cells with multiple options
  resolve.cells <- names(tmp.cl)[sapply(tmp.cl, length) > 1]
  ### cells with only one option. Not further job.
  mapped.cells <- setdiff(mapped.cells, resolve.cells)
  if (length(mapped.cells) > 0) {
    mapped.cl <- c(mapped.cl, setNames(unlist(tmp.cl[mapped.cells]), mapped.cells))
  }
  ### Resolve further options.
  if (length(resolve.cells) > 0) {
    tmp.cat <- sapply(tmp.cl[resolve.cells], function(x) paste(x, collapse = " "))
    for (cat in unique(tmp.cat)) {
      tmp.cl <- unlist(strsplit(cat, " "))
      select.cells <- names(tmp.cat)[tmp.cat == cat]
      mapped.cl <- c(mapped.cl, resolve_cl(cl.g[tmp.cl],
        cl.med, markers, dat, map.dat, select.cells,
        p = p, low.th = low.th
      ))
    }
  }
  return(mapped.cl)
}


#' Compute cluster sums for each row in a matrix
#'
#' This is the scrattch.hicat version of this function (the scrattch.bigcat version crashes the code).
#'
#' @import scrattch.hicat
#' @import MatrixGenerics
#' @import randomForest
#' @import doMC
#' @import foreach
#'
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with sums for each cluster
#' 
#' @keywords internal
get_cl_sums <- function(mat, 
                        cl) {
  
  cl.mat <- get_cl_mat(cl)
  if(all(names(cl) %in% colnames(mat))){
    cl.sums <- Matrix::tcrossprod(mat[,rownames(cl.mat)], Matrix::t(cl.mat))
  }
  else{
    cl.sums <- Matrix::crossprod(mat[rownames(cl.mat),], cl.mat)
  }
  cl.sums <- as.matrix(cl.sums)
  return(cl.sums)
}

#' Compute cluster means for each row in a matrix
#' 
#' This is the scrattch.hicat version of this function (the scrattch.bigcat version crashes the code).
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with means for each cluster
#' 
#' @keywords internal
get_cl_means <- function (mat, cl) 
{
  cl.sums <- get_cl_sums(mat, cl)
  cl.size <- table(cl)
  cl.means <- as.matrix(Matrix::t(Matrix::t(cl.sums)/as.vector(cl.size[colnames(cl.sums)])))
  return(cl.means)
}

#' Compute cluster medians for each row in a matrix
#' 
#' @param mat A gene (rows) x samples (columns) sparse matrix
#' @param cl A cluster factor object
#' 
#' @return a matrix of genes (rows) x clusters (columns) with medians for each cluster
#'
#' @keywords external
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

# #' Strip extra annotation information from dendrogram
# #'
# #' @param dend R dendrogram object
# #' @param value Vector of values pulled from the dendrogram
# #' @param attribute Which attribute should be overwritten
# #'
# #' @return R dendrogram object with updated attributes
# #'
# #' @keywords internal
# revert_dend_label <- function(dend, value, attribute="label")
# {
#   if(attr(dend, attribute)=="")
#     attr(dend, attribute) <- value[attr(dend,"original_label")]
#   if (length(dend)>1) for(i in 1:length(dend))
#     dend[[i]]=revert_dend_label(dend[[i]], value=value, attribute)
#   return(dend)
# }

# #' map_dend_membership
# #'
# #' @param dend R dendrogram in a specific format
# #' @param cl.dat gene by cell type matrix (I think?)
# #' @param map.dat normalized data of the MAPPING data set.
# #' @param map.cells names of cells to map (e.g., the column names of the cell x gene matrix)
# #' @param mc.cores number of cores to run the mapping on 
# #' @param bs.num Number of bootstrapping runs to calculate the membership from (default = 100)
# #' @param seed = random seed
# #' @param ... other variables to pass to map_dend
# #'
# #' @import foreach
# #'
# #' @return membership table
# #' 
# #' @keywords internal
# map_dend_membership <-
#   function(dend,
#            cl.dat,
#            map.dat,
#            map.cells,
#            mc.cores = 10,
#            bs.num = 100,
#            seed = 42,
#            ...)
#   {
    
#     # Optional libraries for UNIX parallel implementation (likely will crash in Windows)
#     if (mc.cores>1) {
#       suppressPackageStartupMessages({
#         library(doMC)
#         library(parallel)
#       })
#     } 
    
#     if(mc.cores ==1){
#       registerDoSEQ()
#     }else{
#       registerDoMC(cores=mc.cores)
#     }
#     mem = foreach(i = 1:bs.num, .combine = 'c') %dopar% {
#       print(i)
#       mfishtools::map_dend(dend, cl.dat, map.dat, map.cells, seed=i)
#     }
#     memb = data.frame(cell = names(mem), cl = mem)
#     memb = table(memb$cell, memb$cl)
#     memb = memb / bs.num
#     tmp = get_nodes_attr(dend, "label")
#     tmp = tmp[tmp %in% colnames(memb)]
#     memb = memb[, tmp]
#     return(memb)
#   }

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
# map_dend <-
#   function(dend,
#            cl.dat,
#            map.dat,
#            select.cells=colnames(map.dat),
#            p = 0.8,
#            low.th = 0.1,
#            default.markers = NULL,
#            seed = 42)
#   {
#     final.cl = c(setNames(rep(
#       attr(dend, "label"), length(select.cells)
#     ), select.cells))
#     if (length(dend) <= 1) {
#       return(final.cl)
#     }
#     markers = attr(dend, "markers")
#     markers = markers[names(markers) %in% row.names(map.dat)]
#     cl.g = sapply(dend, labels, simplify = F)
#     names(cl.g) = 1:length(cl.g)
#     select.cl <- cl[cl %in% unlist(cl.g)]
#     ### Sampling the cells from the reference cluster
#     cells <- unlist(tapply(names(select.cl), select.cl, function(x) sample(x, round(length(x) * p))))
#     genes = names(markers)
#     genes = union(genes, default.markers)
#     ### Compute reference cluster median based on
#     ### subsampled cells
#     cl.med <- do.call("cbind", tapply(
#       cells, droplevels(cl[cells]),
#       function(x) rowMedians(dat[genes, x, drop = F])
#     ))
#     row.names(cl.med) <- genes
#     ### determine which branch to take.
#     mapped.cl = resolve_cl(cl.g,
#                            cl.dat,
#                            markers,
#                            map.dat,
#                            select.cells,
#                            p = p,
#                            low.th = low.th,
#                            seed = seed+1)
#     if (length(mapped.cl) > 0) {
#       for (i in unique(mapped.cl)) {
#         print(i)
#         select.cells = names(mapped.cl)[mapped.cl == i]
#         if (length(select.cells) > 0) {
#           final.cl = c(
#             final.cl,
#             map_dend(
#               dend[[as.integer(i)]],
#               cl.dat,
#               map.dat,
#               select.cells,
#               p = p,
#               low.th = low.th,
#               seed = seed+2
#             )
#           )
#         }
#       }
#       return(cl = final.cl)
#     }
# }

# #' resolve_cl
# #'
# #' @param cl.g Cluster labels in some format
# #' @param cl.med Cluster medians
# #' @param markers Genes to use as markers for this function
# #' @param map.dat normalized data of the MAPPING data set.
# #' @param select.cells names of cells to map (e.g., the column names of the cell x gene matrix)
# #' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
# #' @param low.th the minimum difference in Pearson correlation required to decide on which branch
# #' @param seed - random seed for reproducibility
# #'
# #' @return mapped.cl output
# #' 
# #' @keywords internal
# resolve_cl <-
#   function(cl.g,
#            cl.dat,
#            markers,
#            map.dat,
#            select.cells,
#            p = 0.8,
#            low.th = 0.1,
#            seed = 42)
#   {
#     ##
#     genes = names(markers)[markers > 0]
#     tmp.cl = unlist(cl.g)
    
#     ###For each branch point, find the highest expression cluster.
#     tmp.med = sapply(cl.g, function(g) rowMaxs(cl.dat[genes, g, drop = F]))
#     row.names(tmp.med) = genes
#     ###Make sure the genes are discriminative between all the branches.
#     genes = genes[rowMaxs(tmp.med) - rowMins(tmp.med) > 1]
    
#     ###Sample the markers based on the weights.
#     ##TO DO: randomforest sometimes give importance value of 0. adjust for that.
#     set.seed(seed)
#     seed  = seed+1
#     genes = sample(genes, round(length(genes) * p), prob = markers[genes])
    
#     ###Compute the correlation with the median cluster profile.
#     ###add drop=F
#     cl.cor = cor(as.matrix(map.dat[genes, select.cells, drop = F]), cl.dat[genes, tmp.cl, drop =
#                                                                              F])
#     cl.cor[is.na(cl.cor)] = 0
#     ###Compute the best match in each branch.
#     tmp.score = do.call("cbind", sapply(cl.g, function(x)
#       rowMaxs(cl.cor[, x, drop = F]), simplify = F))
#     row.names(tmp.score) = row.names(cl.cor)
#     ####Determine the best match.
#     best.score = setNames(rowMaxs(tmp.score), row.names(tmp.score))
#     ###determine the difference from the best match.
#     diff.score = best.score - tmp.score
    
#     ####Give up on cells can't be discriminated,choose one branch randomly.
#     unresolved.cl = row.names(tmp.score)[rowSums(diff.score < low.th) ==
#                                            ncol(diff.score)]
#     set.seed(seed)
#     seed  = seed+1
#     mapped.cl = setNames(sample(colnames(tmp.score), length(unresolved.cl), replace =
#                                   T), unresolved.cl)
    
#     ###Cells mapped to one or more branches.
#     mapped.cells = setdiff(row.names(cl.cor), unresolved.cl)
#     ###For binary branch, done already
#     if (length(cl.g) == 2) {
#       mapped.cl = c(mapped.cl, setNames(colnames(diff.score)[apply(diff.score[mapped.cells, , drop =
#                                                                                 F], 1, which.min)], mapped.cells))
#       return(mapped.cl)
#     }
#     ##The remaining options for mapped cells
#     tmp.cl = sapply(mapped.cells, function(x)
#       colnames(diff.score)[which(diff.score[x,] < low.th)], simplify = F)
#     ###cells with multiple options
#     resolve.cells = names(tmp.cl)[sapply(tmp.cl, length) > 1]
#     ###cells with only one option. Not further job.
#     mapped.cells = setdiff(mapped.cells, resolve.cells)
#     if (length(mapped.cells) > 0) {
#       mapped.cl = c(mapped.cl, setNames(unlist(tmp.cl[mapped.cells]), mapped.cells))
#     }
#     ###Resolve further options.
#     if (length(resolve.cells) > 0) {
#       tmp.cat = sapply(tmp.cl[resolve.cells], function(x)
#         paste(x, collapse = " "))
#       for (cat in unique(tmp.cat)) {
#         tmp.cl = unlist(strsplit(cat, " "))
#         select.cells = names(tmp.cat)[tmp.cat == cat]
#         mapped.cl = c(
#           mapped.cl,
#           resolve_cl(
#             cl.g[tmp.cl],
#             cl.dat,
#             markers,
#             map.dat,
#             select.cells,
#             p = p,
#             low.th = low.th
#           )
#         )
#       }
#     }
#     return(mapped.cl)
#   }
