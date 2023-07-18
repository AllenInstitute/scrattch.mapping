#' Function to set mapping mode
#'
#' @param AIT.anndata A vector of cluster names in the reference taxonomy.
#' @param mode Number of cells to keep per cluster.
#'
#' @return AIT anndata with mode set for mapping
#' 
#' @export
file.path <- function(AIT.anndata, mode){
  if(!mode %in% names(AIT.anndata$uns$filter)){ stop(paste0(mode, " is invalid. Choose from: standard, patchseq")) }
  AIT.anndata$uns$mode = mode
  return(AIT.anndata)
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
  if(is.null(AIT.anndata$uns$QC_markers[[mode]]$allMarkers)){
    QC_marker_file = file.path(AIT.anndata$uns$taxonomyDir, mode, "QC_markers.rda")
    if(file.exists(QC_marker_file)){
      ## 
      print("Converting taxonomy .h5ad to new format and saving. This should only happen once per taxonomy.")
      ##
      load(QC_marker_file)
      AIT.anndata$uns$QC_markers[[mode]]$allMarkers = allMarkers
      AIT.anndata$uns$QC_markers[[mode]]$markers    = markers
      AIT.anndata$uns$QC_markers[[mode]]$countsQC   = countsQC
      AIT.anndata$uns$QC_markers[[mode]]$cpmQC      = cpmQC
      AIT.anndata$uns$QC_markers[[mode]]$classBr    = classBr
      AIT.anndata$uns$QC_markers[[mode]]$subclassF  = subclassF
      AIT.anndata$uns$QC_markers[[mode]]$qc_samples = colnames(countsQC) # since colnames are lost
      AIT.anndata$uns$QC_markers[[mode]]$qc_genes   = rownames(countsQC) # since rownames are lost
      AIT.anndata$write_h5ad(file.path(AIT.anndata$uns$taxonomyDir, "AI_taxonomy.h5ad"))
    }else{
      stop("Could not find QC marker files required for taxonomy mode. Please run buildPatchseqTaxonomy and try again.")
    }
  }
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
rfTreeMapping <- function (dend, refDat, clustersF, mapDat = refDat, p = 0.8, 
                           low.th = 0.1, bootstrap = 100, seed = 1) 
{
  genes <- intersect(rownames(refDat), rownames(mapDat))
  refDat <- as.matrix(refDat)[genes, ]
  mapDat <- as.matrix(mapDat)[genes, ]
  pseq.cells <- colnames(mapDat)
  pseq.mem <- sapply(1:bootstrap, function(i) {
    j <- i
    go <- TRUE
    while (go) {
      j <- j + 1000
      set.seed(j + seed)
      tmp <- try(mfishtools::map_dend(dend, clustersF, refDat, mapDat, pseq.cells, 
                                      p = p, low.th = low.th), silent=TRUE)
      if (length(tmp) > 1) 
        go <- FALSE
    }
    tmp
  }, simplify = F)
  memb <- unlist(pseq.mem)
  memb <- data.frame(cell = names(memb), cl = memb)
  memb$cl <- factor(memb$cl, levels = get_nodes_attr(dend, 
                                                     "label"))
  memb <- table(memb$cell, memb$cl)
  memb <- memb/bootstrap
  return(memb)
}


####################################################################
## Functions for reversing '\' and '/'

#' Sets the default leading_string for file.path()
#' 
#' @param leading_string Default (NULL) sets to "\\\\" for Windows and "/" otherwise; or can provide any character vector
#'
#' @return leading_string for use as default in file.path().
#'
#' @export
setLeadingString <- function(leading_string=NULL){
  if (is.null(leading_string)){
    leading_string="/"
    if(get_os()=="windows") leading_string="\\\\" 
  }
  options("leading_string"=leading_string)
}


#' Sets the default path_separator for file.path()
#' 
#' @param path_separator Default (NULL) sets to .Platform$file.sep or can provide any character vector
#'
#' @return path_separator for use as default in file.path().
#'
#' @export
setPathSeparator <- function(path_separator=NULL){
  if (is.null(path_separator)){
    path_separator = .Platform$file.sep
  }
  options("path_separator"=path_separator)
}


#' Detect the operating system
#' 
#' This function was taken directly from https://conjugateprior.org/2015/06/identifying-the-os-from-r/ and all credit goes to Will Lowe from "conjugateprior".
#' 
#' @keywords internal
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

#' Construct Path to File across platforms
#'
#' Construct the path to a file from components in a platform-independent way. This version is a wrapper of the base `file.path` function which also reverses '/' direction and can attempt to add double slashes if needed.
#'
#' @param ... character vectors. Long vectors are not supported.
#' @param path_separator the path separator to use (assumed to be ASCII).
#' @param leading_string what is the leading character(s) (e.g., '/' or '\\\\'). A string can be provided, or by default if the "leading_string" global variable is set it takes that variable, otherwise this is guessed at based on operating system and/or existence of a file at the file path
#' @param change_chars which characters should be changes to the path_separator value (default NULL = "none"). If you want to change all slashes in the file path to the correct direction set change_chars = c("/","\\"))
#'
#' @return A file path with slashes going the correct direction.
#'
#' @export
file.path <- function (...,path_separator = getOption("path_separator"),leading_string=getOption("leading_string"),change_chars=NULL){
  if(is.null(path_separator))
    path_separator = .Platform$file.sep
  path <- base::file.path(...,path_separator)
  first_character <- grep('[^[:punct:]]', strsplit(path,"")[[1]])[1]
  path <- substring(path,first_character,nchar(path))
  if (!is.null(change_chars)){
    path  <- strsplit(path,"")[[1]]
    slash <- which(is.element(path,change_chars)) 
    path[slash] <- path_separator
    path <- paste0(path,collapse="")
  }
  if (is.null(leading_string)){
    leading_string="/"
    if(get_os()=="windows") leading_string="\\\\" 
  }
  #remove trailing slashes
  while((substring(path,nchar(path),nchar(path)))==path_separator){
    path = substring(path,1,nchar(path)-1)
  }
  paste0(leading_string,path)
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

#' Compute cluster sums for each row in a matrix
#' 
#' This is the scrattch.hicat version of this function (the scrattch.bigcat version crashes the code).
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
    }else{
      registerDoMC(cores=mc.cores)
    }
    mem = foreach(i = 1:bs.num, .combine = 'c') %dopar% {
      print(i)
      map_dend(dend, cl.dat, map.dat, map.cells, seed=i)
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
    tmp.med = sapply(cl.g, function(g) rowMaxs(cl.dat[genes, g, drop = F]))
    row.names(tmp.med) = genes
    ###Make sure the genes are discriminative between all the branches.
    genes = genes[rowMaxs(tmp.med) - rowMins(tmp.med) > 1]
    
    ###Sample the markers based on the weights.
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