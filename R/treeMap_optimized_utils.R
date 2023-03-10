#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
l2norm <- function(X, by="column")
{
  if (by=="column") {
    l2norm <- sqrt(Matrix::colSums(X^2))
    if (!any(l2norm==0)) {
      X=sweep(X, 2, l2norm, "/", check.margin=FALSE)
    }
    else{
      warning("L2 norms of zero detected for distance='Cosine, no transformation")
    }
    X = X 
  } else {
    l2norm <- sqrt(Matrix::rowSums(X^2))
    if (!any(l2norm==0)) {      
      X= X/l2norm
    }
    else{
      warning("L2 norms of zero detected for distance='Cosine'")
      X = X/ pmax(l2norm,1)
    }
  }
}

#' Get KNN
#'
#' @param dat 
#' @param ref.dat 
#' @param k 
#' @param method 
#' @param dim 
#'
#' @return
#'
#' @keywords internal
get_knn <- function(dat, ref.dat, k, method ="cor", dim=NULL,index=NULL, build.index=FALSE, transposed=TRUE, return.distance=FALSE)
  {
    if(transposed){
      cell.id = colnames(dat)
    }
    else{
      cell.id= row.names(dat)
    }    
    if(transposed){
      if(is.null(index)){
        ref.dat = Matrix::t(ref.dat)
      }
      dat = Matrix::t(dat)
    }
    if(method=="RANN"){
      library(RANN)
      knn.result = RANN::nn2(ref.dat, dat, k=k)
    }
    else if(method %in% c("Annoy.Euclidean", "Annoy.Cosine","cor")){
      library(BiocNeighbors)      
      if(is.null(index)){
        if(method=="cor"){
          ref.dat = ref.dat - rowMeans(ref.dat)
          ref.dat = l2norm(ref.dat,by = "row")
        }
        if (method=="Annoy.Cosine"){
          ref.dat = l2norm(ref.dat,by = "row")
        }
        if(build.index){
          index= buildAnnoy(ref.dat)
        }
      }
      if (method=="Annoy.Cosine"){
        dat = l2norm(dat,by="row")
      }
      if (method=="cor"){
        dat = dat - rowMeans(dat)
        dat = l2norm(dat,by = "row")
      }
      knn.result = queryAnnoy(X= ref.dat, query=dat, k=k, precomputed = index)
    }
    else{
      stop(paste(method, "method unknown"))
    }    
    knn.index= knn.result[[1]]
    knn.distance = knn.result[[2]]
    row.names(knn.index) = row.names(knn.distance)=cell.id
    if(!return.distance){
      return(knn.index)
    }
    else{
      list(index=knn.index, distance=knn.distance)
    }
  }


#' get knn batch
#'
#' @param dat 
#' @param ref.dat 
#' @param k 
#' @param method 
#' @param dim 
#' @param batch.size 
#' @param mc.cores 
#'
#' @return
#'
#' @keywords internal
get_knn_batch <- function(dat, ref.dat, k=1, method="cor", dim=NULL, batch.size, mc.cores=1, return.distance=FALSE,...)
  {
    if(return.distance){
      fun = "knn_combine"
    }
    else{
      fun = "rbind"
    }
    results <- batch_process(x=1:ncol(dat), batch.size=batch.size, mc.cores=mc.cores, .combine=fun, FUN=function(bin){
      get_knn(dat=dat[row.names(ref.dat),bin,drop=F], ref.dat=ref.dat, k=k, method=method, dim=dim,return.distance=return.distance, ...)
    })
    return(results)
  }

knn_combine <- function(result.1, result.2)
{
  knn.index = rbind(result.1[[1]], result.2[[1]])
  knn.distance = rbind(result.1[[2]], result.2[[2]])
  return(list(knn.index, knn.distance))
}

#' Batch process
#'
#' @param x 
#' @param batch.size 
#' @param FUN 
#' @param mc.cores 
#' @param .combine 
#' @param ... 
#'
#' @return
#'
#' @keywords internal
batch_process <- function(x, batch.size, FUN, mc.cores=1, .combine="c",...)
  {
    require(foreach)
    require(doMC)
    if (mc.cores == 1) {
      registerDoSEQ()
    }
    else {
      registerDoMC(cores=mc.cores)
      #on.exit(parallel::stopCluster(), add = TRUE)
    }
    bins = split(x, floor((1:length(x))/batch.size))
    results= foreach(i=1:length(bins), .combine=.combine) %dopar% FUN(bins[[i]],...)
    return(results)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
build_train_index <- function(cl.dat, method= c("Annoy.Cosine","cor","Annoy.Euclidean"),fn=tempfile(fileext=".idx"))
  {
    library(BiocNeighbors)
    method = method[1]
    ref.dat = Matrix::t(cl.dat)
    if(method=="cor"){
      ref.dat = ref.dat - rowMeans(ref.dat)
      ref.dat = l2norm(ref.dat,by = "row")
    }
    if (method=="Annoy.Cosine"){
      ref.dat = l2norm(ref.dat,by = "row")
    }
    index= buildAnnoy(ref.dat, fname=fn)
    return(index)    
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
build_train_index_bs <- function(cl.dat, method= c("Annoy.Cosine","cor","Annoy.Euclidean"),sample.markers.prop=0.8, iter=100, mc.cores=10,fn=tempfile(fileext=".idx"))
  {
    library(BiocNeighbors)
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)

    if (is.na(sample.markers.prop)) {
       index = build_train_index(cl.dat = cl.dat, method=method, fn=paste0(fn, ".",1))
       return(list(list(cl.dat=cl.dat, index=index)))
    } else {
       ###for each cluster, find markers that discriminate it from other types
       train.dat <- foreach(i=1:iter, .combine="c") %dopar% {
         train.markers = sample(row.names(cl.dat), round(nrow(cl.dat) * sample.markers.prop))
         train.cl.dat = cl.dat[train.markers,]
         index = build_train_index(cl.dat = train.cl.dat, method=method, fn = paste0(fn, ".",i))
         return(list(list(cl.dat=train.cl.dat, index=index)))
       }   
    }
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
map_cells_knn <- function(topk=1,test.dat, cl.dat, train.index=NULL, method = c("Annoy.Cosine","cor"), batch.size=5000, mc.cores=1)
  {

    cl.knn = get_knn_batch(test.dat, cl.dat, k=topk, index=train.index, method=method, transposed=TRUE, batch.size=batch.size, mc.cores=mc.cores,return.distance=TRUE)
    knn.index = cl.knn[[1]]
    knn.dist = cl.knn[[2]]
    map.df = data.frame(sample_id=colnames(test.dat), cl = colnames(cl.dat)[knn.index], dist = knn.dist)
    return(map.df)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
map_cells_knn_big <- function(big.dat, cl.dat, select.cells, train.index=NULL, method = c("Annoy.Cosine","cor"), batch.size=10000, mc.cores=10)
  {    
    cl.knn =  get_knn_batch_big(big.dat, cl.dat, select.cells=select.cells, k=1, index=train.index, method=method, transposed=TRUE, batch.size=batch.size, mc.cores=mc.cores,return.distance=TRUE)
    knn.index = cl.knn[[1]]
    knn.dist = cl.knn[[2]]
    map.df = data.frame(sample_id=select.cells, cl = colnames(cl.dat)[knn.index], dist = knn.dist)
    return(map.df)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
map_cells_knn_bs <- function(topk=1, test.dat, iter=100,cl.dat=NULL,train.index.bs=NULL, method = c("Annoy.Cosine","cor"), mc.cores=20, ...)
  {
    require(doMC)
    require(foreach)
    mc.cores = min(mc.cores, length(train.index.bs))
    registerDoMC(cores=mc.cores)
    ###for each cluster, find markers that discriminate it from other types
    if(!is.null(train.index.bs)){
      iter = length(train.index.bs)
    } else{
      idx = match(rownames(cl.dat), rownames(test.dat))
      if (any(is.na(idx))) {
         print("some genes in the train data set are missing from query data")
         common.gene = intersect(rownames(cl.dat), rownames(test.dat))
         cl.dat = cl.dat[common.gene,]
      }
      train.index.bs = build_train_index_bs(cl.dat, method=method,iter=iter, ...)
    }
    library(data.table)
    map.list <- foreach(i=1:iter, .combine="c") %dopar% {
      train.index = train.index.bs[[i]]$index
      cl.dat = train.index.bs[[i]]$cl.dat
      map.df=map_cells_knn(topk, test.dat, cl.dat, train.index, method = c("Annoy.Cosine","cor"))
      map.df = list(map.df)
    }
    map.df = rbindlist(map.list)
    map.df = map.df %>% group_by(sample_id, cl) %>% summarize(freq=n(),dist = mean(dist))
    map.df$freq = map.df$freq/iter
    best.map.df = map.df %>% group_by(sample_id) %>% summarize(best.cl= cl[which.max(freq)],prob=max(freq), avg.dist = dist[which.max(freq)])
    if(method=="cor"){
      best.map.df = best.map.df%>% mutate(avg.cor = 1 - avg.dist^2/2)
    }
    return(list(map.freq=map.df, best.map.df = best.map.df))    
  }

## cl.list is cluster membership at different levels, finest at the beginning.
## val is a vector associated with sample_id
## compute z_score aggregate at different levels of clustering, start with finest level of clustering, and resort to higher level if not enough sample size
#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
z_score <- function(cl.list, val, min.samples =100)
  {

    sample_id = names(cl.list[[1]])
    z.score = c()
    for(i in 1:length(cl.list)){
      cl=cl.list[[i]][sample_id]
      cl.size = table(cl)
      if(i !=length(cl.list)){
        select.cl = names(cl.size)[cl.size > min.samples]
      }
      else{
        select.cl = names(cl.size)
      }
      df = data.frame(sample_id = names(cl),cl=cl, val=val[names(cl)])      
      df = df %>% filter(cl %in% select.cl) %>% group_by(cl) %>% mutate(z = (val - mean(val))/sd(val))
      z.score[df$sample_id] = df$z
      sample_id = setdiff(sample_id, df$sample_id)      
    }
    return(z.score)
  }

  get_gene_score_ds <- function(ds, to.add, genes, cl.bin, de=NULL, max.num=1000,mc.cores=20)
  {
    require(doMC)
    registerDoMC(cores=mc.cores)
    if(!is.null(de)){
      tmp.de = suppressMessages(de %>% right_join(to.add[,c("P1","P2")]))
      gene.score = tmp.de %>% group_by(gene) %>% summarize(score = sum(as.numeric(max.num- rank))) %>% filter(gene %in% genes) %>% arrange(-score)   
    }
    else{
      to.add = suppressMessages(to.add %>% left_join(cl.bin,by=c("P1"="cl")) %>% left_join(cl.bin,by=c("P2"="cl")))
      cl.bin.x = to.add %>% pull(bin.x) %>% unique
      cl.bin.y = to.add %>% pull(bin.y) %>% unique  
      tmp=foreach::foreach(bin1=cl.bin.x,.combine="c")%:%
        foreach::foreach(bin2=cl.bin.y,.combine="c")%dopar% {
          cat(bin1, bin2, '\r')
          tmp.pairs = to.add %>%  filter(bin.x == bin1 & bin.y ==bin2)
          tmp.de = ds %>% filter(bin.x == bin1 & bin.y ==bin2 & gene %in% genes & rank < max.num & P1 %in% tmp.pairs$P1 & P2 %in% tmp.pairs$P2) %>% collect
          tmp.de = suppressMessages(tmp.de %>% right_join(tmp.pairs,by=c("P1","P2")))
          gene.score = tmp.de %>% group_by(gene) %>% summarize(rank.sum = sum(max.num- rank))
          list(gene.score)
        }
      gene.score=rbindlist(tmp)
      gene.score = gene.score %>% group_by(gene) %>% summarize(score=sum(as.numeric(rank.sum))) %>% arrange(-score)
    }
    return(gene.score)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
update_gene_score_ds <- function(gene.score, ds, to.remove, cl.bin, de=NULL, max.num=1000,mc.cores=20)
  {
    require(doMC)
    registerDoMC(cores=mc.cores)
    to.remove = suppressMessages(to.remove %>% left_join(cl.bin,by=c("P1"="cl")) %>% left_join(cl.bin,by=c("P2"="cl")))
    cl.x = to.remove %>% pull(P1) %>%unique
    cl.y = to.remove %>% pull(P2) %>%unique
    
    cl.bin.x = to.remove %>% pull(bin.x) %>% unique
    cl.bin.y = to.remove %>% pull(bin.y) %>% unique
    if(!is.null(de)){
      tmp.de = suppressMessages(de %>% right_join(to.remove))
      rm.gene.score = tmp.de %>% group_by(gene) %>% summarize(rm.score = sum(max.num- rank))
    }
    else{
      if(length(cl.bin.x)*length(cl.bin.y) * nrow(gene.score) < 10^6){
        de=  ds %>% filter(bin.x %in% cl.bin.x & bin.y %in% cl.bin.y & gene %in% gene.score$gene & rank < max.num & P1 %in% cl.x & P2 %in% cl.y) %>% collect
        tmp.de = suppressMessages(de %>% right_join(to.remove))
        rm.gene.score = tmp.de %>% group_by(gene) %>% summarize(rm.score = sum(max.num- rank))        
      }
      else{
        tmp=foreach::foreach(bin1=cl.bin.x,.combine="c")%:%
          foreach::foreach(bin2=cl.bin.y,.combine="c")%dopar% {
            tmp.pairs = to.remove %>%  filter(bin.x == bin1 & bin.y ==bin2)
            tmp.de = ds %>% filter(bin.x == bin1 & bin.y ==bin2 & gene %in% gene.score$gene & rank < max.num & P1 %in% tmp.pairs$P1 & P2 %in% tmp.pairs$P2) %>% collect
            tmp.de = suppressMessages(tmp.de %>% right_join(tmp.pairs))
            gene.score = tmp.de %>% group_by(gene) %>% summarize(rank.sum = sum(max.num- rank))
            list(gene.score)
          }
        rm.gene.score=rbindlist(tmp)      
        if(is.null(rm.gene.score)|nrow(rm.gene.score)==0){
          return(gene.score)
        }
        rm.gene.score = rm.gene.score %>% group_by(gene) %>% summarize(rm.score=sum(as.numeric(rank.sum)))
      }
    }
    tmp = suppressMessages(gene.score %>% left_join(rm.gene.score) %>% mutate(rm.score = ifelse(is.na(rm.score), 0, rm.score)) %>% mutate(new.score = score - rm.score))
    tmp = tmp%>% select(gene, new.score) %>% rename(score=new.score) %>% filter(score > 0) %>% arrange(-score)
    return(tmp)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
check_pairs_ds <- function(de.dir, to.add, genes,cl.bin, de=NULL, mc.cores=10,max.num=1000)
  {
    require(doMC)
    registerDoMC(cores=mc.cores)
    to.add = to.add[,c("P1","P2")]
    to.add = suppressMessages(to.add %>% left_join(cl.bin,by=c("P1"="cl")) %>% left_join(cl.bin,by=c("P2"="cl")))
    select.P1 = to.add %>% pull(P1) %>% unique
    select.P2 = to.add %>% pull(P2) %>% unique    
    if(!is.null(de)){
      de.checked = suppressMessages(to.add %>% left_join(de %>% filter(rank < max.num & gene %in% genes)) %>% group_by(P1,P2) %>% summarize(checked=sum(!is.na(gene))))
    }
    else{
      cl.bin.x = to.add %>% pull(bin.x) %>% unique
      cl.bin.y = to.add %>% pull(bin.y) %>% unique
      if(length(cl.bin.x)*length(cl.bin.y)*length(genes) < 10^5){
        ds = open_dataset(de.dir)
        de = ds %>% filter(bin.x  %in% cl.bin.x & bin.y %in% cl.bin.y & gene %in% genes & P1 %in% select.P1 & P2 %in% select.P2 & rank < max.num ) %>% collect
        de.checked = suppressMessages(to.add %>% left_join(de %>% filter(rank < max.num & gene %in% genes)) %>% group_by(P1,P2) %>% summarize(checked=sum(!is.na(gene))))
      }
      else{
        de.checked=foreach::foreach(bin1=cl.bin.x,.combine="c")%:%
          foreach::foreach(bin2=cl.bin.y,.combine="c")%dopar% {
            cat("bin", bin1, bin2, "\r")
            tmp.pairs = to.add %>%  filter(bin.x == bin1 & bin.y ==bin2)
            ###
            #tmp.de =ds %>% filter(bin.x == bin1 & bin.y ==bin.y) %>% collect %>% filter(gene %in% genes & rank < max.num)
            #Directly read from parquet file is faster
            d = file.path(de.dir, paste0("bin.x=",bin1), paste0("bin.y=",bin2))
            fn = dir(d, pattern="parquet")
            tmp.de = read_parquet(file.path(d, fn))
            tmp.de = tmp.de %>% filter(gene %in% genes & rank < max.num)
            tmp.de = suppressMessages(tmp.de %>% right_join(tmp.pairs))
            de.checked = suppressMessages(tmp.de %>% group_by(P1,P2) %>% summarize(checked=n()))
            list(de.checked)
          }      
        de.checked = rbindlist(de.checked)
      }
    }
    return(de.checked)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_markers_ds <- function(de.dir, cl.bin, select.cl=NULL, top.n=20,mc.cores=10)
  {
    ds = open_dataset(de.dir)

    if(!is.null(select.cl)){
      cl.bin = cl.bin %>% filter(cl %in% select.cl)
    }
    select.bin = cl.bin %>% pull(bin) %>% unique
    mc.cores=min(mc.cores, length(select.bin))
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)    
    tmp=foreach::foreach(bin1=select.bin,.combine="c")%:%
      foreach::foreach(bin2=select.bin,.combine="c")%dopar% {
        de = ds %>% filter(bin.x %in% bin1 & bin.y %in% bin2)
        if(is.null(select.cl)){
          de = de %>% filter(P1 %in% select.cl & P2 %in% select.cl)            
        }
        de %>% filter(rank <= top.n) %>% pull(gene) %>% unique
      }
    select.markers=unique(tmp)
    return(select.markers)
  }


#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_markers_pair_direction_ds <- function(de.dir, add.num, genes, cl.bin, de=NULL, mc.cores=mc.cores,max.genes=1000,byg=10, ...)
  {
    ds = open_dataset(de.dir)
    select.genes=c()
    if(!is.null(de)){
      de = de %>% filter(gene %in% genes)
    }
    add.num = suppressMessages(add.num %>% left_join(cl.bin,by=c("P1"="cl")) %>% left_join(cl.bin,by=c("P2"="cl")))
    gene.score= get_gene_score_ds(ds, to.add=add.num, genes=genes,cl.bin=cl.bin, de=de,mc.cores=mc.cores,...)
    print(dim(gene.score))
    while(nrow(add.num)>0 & length(genes)>0 & length(select.genes)< max.genes){      
      if(is.null(gene.score) | nrow(gene.score)==0){
        break
      }
      if (length(gene.score$gene)> byg) g = as.character(gene.score$gene[1:byg])
      else g = as.character(gene.score$gene)

      new.checked = check_pairs_ds(de.dir, to.add=add.num %>% select(P1,P2), genes=g,cl.bin=cl.bin,de=de, ...)
      add.num$checked=NULL
      add.num = suppressMessages(add.num %>% left_join(new.checked, by=c("P1","P2")))      
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)
      to.remove = add.num %>% filter(num <=0) %>% select(P1,P2)
      add.num = add.num %>% filter(num > 0)
      genes = setdiff(genes,g)
      select.genes=c(select.genes,g)
      cat('gene_score', dim(gene.score), length(select.genes), '\r')
      gene.score = update_gene_score_ds(gene.score, ds, to.remove, cl.bin, de=de,...)
      gene.score = gene.score %>% filter(gene %in% genes)
      if(!is.null(de)){
        de = de %>% filter(!(gene %in% g))
        if(is.null(de)| nrow(de)==0){
          break
        }
      }
    }
    return(list(select.genes=select.genes, de=de))
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_markers_pair_group_top_ds <- function(g1,g2,ds, genes, cl.bin, select.sign="up", n.markers=20,mc.cores=1, ...)
{
  require(matrixStats)
  require(data.table)
  require(arrow)
  require(dplyr)

  require(doMC)
  registerDoMC(cores=mc.cores)
  up.to.add = down.to.add=NULL
  up.genes=down.genes=NULL
  if("up" %in% select.sign){
    up.to.add = as.data.frame(create_pairs(g1, g2, direction="unidirectional"))
    up.gene.score  = get_gene_score_ds(ds, to.add=up.to.add, genes=genes, cl.bin=cl.bin,...)
    up.genes = head(up.gene.score$gene, n.markers)
  }
  if("down" %in% select.sign){
    down.to.add = as.data.frame(create_pairs(g2, g1,direction="unidirectional"))
    down.gene.score  = get_gene_score_ds(ds, to.add=down.to.add, genes=genes, cl.bin=cl.bin, ...)
    down.genes = head(down.gene.score$gene, n.markers)
  }
  to.add=rbind(up.to.add, down.to.add)
  return(list(up.genes=up.genes, down.genes=down.genes,to.add=to.add))
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_markers_pair_group_ds <- function(g1,g2,de.dir, genes, cl.bin, n.markers=20,select.sign=c("up","down"),max.genes=50,...)
  {
    ds = open_dataset(de.dir)
    result = select_markers_pair_group_top_ds( g1,g2,ds=ds, genes=genes, cl.bin=cl.bin, n.markers=n.markers,select.sign=select.sign,...)
    markers=c(result$up.genes, result$down.genes)
    add.num = result$to.add
    add.num$num = n.markers
    new.checked = check_pairs_ds(de.dir, to.add=add.num %>% select(P1,P2),genes=markers,cl.bin=cl.bin,...)
    add.num = suppressMessages(add.num %>% left_join(new.checked, by=c("P1","P2")))      
    add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)            
    add.num = add.num %>% filter(num > 0)
    max.genes= max.genes - length(markers)
    if(nrow(add.num)>0 & max.genes > 0){
      add.num$checked=NULL
      genes = setdiff(genes,markers)
      more.markers <- select_markers_pair_direction_ds(de.dir, add.num=add.num, genes=genes,cl.bin=cl.bin,max.genes=max.genes,...)
      more.markers <- more.markers$select.genes
      markers= c(markers, more.markers)
    }
    return(markers)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_N_markers_ds<- function(de.dir, select.cl=NULL,pair.num=1, add.num=NULL, genes, cl.bin, default.markers=NULL,...)
  {
    if(is.null(add.num)){
      add.num = as.data.frame(create_pairs(select.cl, direction="directional"))
      add.num$num = pair.num
    }
    if(!is.null(default.markers)){      
      de.checked.num = check_pairs_ds(de.dir, add.num, genes=default.markers, cl.bin=cl.bin,...)
      add.num = suppressMessages(add.num %>% left_join(de.checked.num,by=c("P1","P2")))
      add.num = add.num %>% mutate(checked=ifelse(is.na(checked),0,checked)) %>% mutate(num = num-checked)
      add.add = add.num %>% filter(num > 0)
      genes = setdiff(genes, default.markers)
      add.num$checked=NULL
    }
    markers <- select_markers_pair_direction_ds(de.dir, add.num=add.num, genes=genes,cl.bin=cl.bin,...)$select.genes    
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_pos_markers_ds <- function(de.dir, cl, select.cl, genes, cl.bin, n.markers=1,  mc.cores=1,out.dir="cl.markers",...)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    ds = open_dataset(de.dir)
    if (!dir.exists(out.dir)) {
        dir.create(out.dir)
    }
    
    ###for each cluster, find markers that discriminate it from other types
    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
      #print(x)
      g1=x
      g2 = setdiff(cl, x)
      cl.bin.x = cl.bin %>% filter(cl==g1) %>% pull(bin)      
      select.de = ds %>% filter(bin.x==cl.bin.x & P1==g1 & P2 %in% g2) %>% collect()
      markers <- select_markers_pair_group_ds(g1,g2, de.dir=de.dir, de=select.de, genes=genes, cl.bin=cl.bin, n.markers=n.markers,select.sign="up",...)
      save(markers, file=file.path(out.dir, paste0(x, ".markers.rda")))
      tmp=list(markers)
      names(tmp)=x
      tmp
    }
    return(cl.markers)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_top_pos_markers_ds<- function(ds, cl, select.cl, genes, cl.bin, n.markers=3, mc.cores=10,...)
  {
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    ###for each cluster, find markers that discriminate it from other types

    cl.markers <- foreach(x=select.cl, .combine="c") %dopar% {
      #print(x)
      g1=x
      g2 = setdiff(cl, x)
      cl.bin.x = cl.bin %>% filter(cl==g1) %>% pull(bin)      
      select.de = ds %>% filter(bin.x==cl.bin.x & P1==g1 & P2 %in% g2) %>% collect()
      markers= select_markers_pair_group_top_ds(g1,g2,ds, genes=genes,cl.bin=cl.bin, de=select.de,n.markers=n.markers,select.sign="up",...)$up.genes
      tmp=list(markers)
      names(tmp)=x
      tmp
    }
    return(cl.markers)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_markers_groups_top_ds <- function(ds, cl.group, select.groups=names(cl.group), n.markers=3,mc.cores=1,...)
  {

    library(parallel)
    library(parallel)    
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)
    
    all.cl = unlist(cl.group)
    group.markers <- foreach(x=select.groups, .combine="c") %dopar% {
      #print(x)
      g1 = cl.group[[x]]
      g2 = setdiff(all.cl, g1)              
      markers=select_markers_pair_group_top_ds(g1,g2,ds=ds, select.sign="up",n.markers=n.markers, ...)$up.genes
      list(markers)
    }
    names(group.markers) = select.groups
    return(group.markers)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
select_markers_groups <- function(de.dir, cl.group, genes, cl.bin, select.groups= unique(cl.group$group), n.markers=20,mc.cores=1,byg=10, ...)
  {
    ds = open_dataset(de.dir)
    cl.group$cl = as.character(cl.group$cl)
    group_pair=create_pairs(unique(cl.group$group))
    library(parallel)
    require(doMC)
    require(foreach)
    registerDoMC(cores=mc.cores)

    group.markers <- foreach(i=1:nrow(group_pair), .combine="c") %dopar% {
      x= group_pair[i,1]
      y= group_pair[i,2]
      cat(i, nrow(group_pair), x,' vs ', y,"\n")
      g1 = cl.group %>% filter(group==x) %>% pull(cl)
      g2 = cl.group %>% filter(group==y) %>% pull(cl)
      result=select_markers_pair_group_top_ds(g1,g2,ds=ds, n.markers=n.markers, genes=genes,cl.bin=cl.bin,select.sign=c("up","down"),...)
      list(c(result$up.genes, result$down.genes))
    }
    group.markers=unique(unlist(group.markers))   
    pairs = as.data.frame(create_pairs(cl.group$cl,direction="directional"))
    pairs$pair = row.names(pairs)
    pairs = suppressMessages(pairs %>% left_join(cl.group, by=c("P1"="cl")))
    pairs = suppressMessages(pairs %>% left_join(cl.group, by=c("P2"="cl")))
    pairs = pairs %>% filter(group.x!=group.y) 
   
    registerDoMC(cores=10)
    de.checked.num = check_pairs_ds(de.dir, pairs[,c("P1","P2")], genes=group.markers,cl.bin=cl.bin, mc.cores=10, ...)    
    add.num = suppressMessages(pairs %>% left_join(de.checked.num) %>% mutate(num=n.markers - checked))
    
    select.markers=group.markers
    genes = setdiff(genes, select.markers)
    #save(group.markers, pairs, add.num, de.dir, genes, cl.bin, file="Debug.rda")
    more.markers <- select_markers_pair_direction_ds(de.dir, add.num=add.num,  genes=genes, cl.bin=cl.bin, mc.cores=10, byg=byg, ...)    
    select.markers = c(select.markers, more.markers$select.genes)
    return(select.markers)
  }

#' INFO -- PLEASE ADD --
#'
#' @param in.df
#' @param cl.df
#'
#' @return ???
#'
#' @keywords internal
get_de_genes <- function(ds, cl.bin, cl1, cl2)
  {
    bin1 = cl.bin %>% filter(cl==cl1) %>% pull(bin)
    bin2 = cl.bin %>% filter(cl==cl2) %>% pull(bin)
    select.genes = ds %>% filter(bin.x==bin1 & P1==cl1 & bin.y==bin2 & P2==cl2 | bin.x==bin2 & P2==cl1 & bin.y==bin1 & P1==cl2) %>% collect()
    return(select.genes)
  }


