
#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
update_select.markers_cl.dat <- function(train.list) {
   markers = c()
   Nnodes = length(train.list$marker_index)
   for (i in 1:Nnodes) {
      i.markers = train.list$marker_index[[i]]$marker_tree
      i.Nmarkers = length(i.markers)
      markers = union(markers, i.markers)
      Nmarkers = length(markers)
   }
   print(paste0("select.markers are updated :", Nmarkers, "/", length(train.list$select.markers)))
   train.list$select.markers = markers
   train.list$cl.dat = train.list$cl.dat[match(markers, rownames(train.list$cl.dat)),]
   return(train.list)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
create_initial_df_old <- function(in_sample_name, in_best.cl=NA) {
   N = length(in_sample_name)
   out.df = data.frame(sample_id= in_sample_name, 
                best.cl  = rep(in_best.cl, N),
                prob     = rep(1.0, N),
                avg.dist = rep(0.0, N),
                avg.cor  = rep(0.0, N),
                cl       = rep(in_best.cl, N))
   rownames(out.df) = in_sample_name
   return(out.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
create_initial_df <- function(in_sample_name, in_best.cl=NA) {
   N = length(in_sample_name)
   out.df = data.frame(sample_id= in_sample_name, 
                cl  = rep(in_best.cl, N),
                dist = rep(0.0, N),
                path.cor = rep(0.0, N))
   rownames(out.df) = in_sample_name
   return(out.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
create_initial_membership <- function(in_sample_name, ancestor) {
   Nlevel = ncol(ancestor) 
   Ng = rep(0, Nlevel)
   for (i in 2:Nlevel) Ng[i] = length(unique(ancestor[, i]))
   cumNg = cumsum(Ng)
   names(cumNg) = colnames(ancestor)

   Nsample = length(in_sample_name)
   membership.df = data.frame( cor  = matrix(0, Nsample, cumNg[Nlevel]), 
                            prob = matrix(1, Nsample, cumNg[Nlevel]) )
   return(membership.df)
}   

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
select_ancestor_markers <- function(nlvl, marker_index, ancestor) {

   out <- list()
   for (ll in 1:nrow(ancestor)) {
      markers = c()
      nn.str = ''
      i = 1
      while (i < nlvl) {
         nn = paste0(nn.str, ancestor[ll, i])
         if (nn %in% names(marker_index)) {
            markers = union(markers, marker_index[[nn]]$marker_tree)
         } else {
            print(paste0("markers are not calculated for ", nn))
         }
         nn.str = paste0(nn, ":")
         i=i+1
      }
      out[[as.character(ancestor[ll, nlvl])]] = markers
   }
   return(out) 
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
cor_ancestor_markers <- function( qdat, assigned.df, train.cl.dat, ancestor_markers) {
   cl = assigned.df$cl
   names(cl) = assigned.df$sample_id
   out <- rep(0, length(cl))
   names(out) = names(cl)
   for (i in 1:length(cl)) {
      i.cl = as.character(cl[i])
      i.idx = match(i.cl, colnames(train.cl.dat))
      if (is.na(i.idx)) {
         print(paste("cor_ancestor_marker", i, i.cl, i.idx))
         out[i] = NA
         browser()
      } else {
      cl.markers = intersect(ancestor_markers[[i.cl]], 
                             intersect(rownames(qdat), rownames(train.cl.dat)))
      out[i] = cor(qdat[cl.markers, i], train.cl.dat[cl.markers, i.idx])
      }
   }

   return(out) 
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
call_ANN_cl <- function(iter_i, lvl, index.bs, qdat, prev.df, ancestor, blocksize=50000, mc.cores=10) {
   if (length(index.bs)==1 && is.na(index.bs)) { #no index
      out.df = prev.df
      # propagate best.cl to leaf level 
      # out.df$cl = out.df$best.cl
   } else {
      if (ncol(index.bs$cl.dat)==0) { # no index
         out.df = prev.df
         # propagate best.cl to leaf level 
         # out.df$cl = out.df$best.cl
      } else {
         # knn mapping
         out.df = map_cells_knn( topk=1, qdat, cl.dat=index.bs$cl.dat, train.index=index.bs$index, method="cor", batch.size=blocksize, mc.cores=1 )
         #out.df$prob = prev.df$prob * out.df$prob
      }
   }
   ancestor.df = as.data.frame(ancestor)
   out.df$cl = ancestor.df[match(out.df$cl, ancestor.df[, ncol(ancestor.df)]), (lvl+1)]
   return(out.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
call_membership_cl <- function(lvl, index.bs, qdat, membership.df, ancestor, blocksize=50000, mc.cores=10, topk=3) {
   
   print("call_membership_cl")
   if (length(index.bs)==1 && is.na(index.bs)) { #no index
      out.df = membership.df
      # propagate best.cl to leaf level 
      # out.df$cl = out.df$best.cl
   } else {
      if (ncol(index.bs[[1]]$cl.dat)==0) { # no index
         out.df = membership.df
         # propagate best.cl to leaf level 
         # out.df$cl = out.df$best.cl
      } else {
         # knn mapping
         print(paste0("calling map_cells_knn_bs, with topk",  topk))
         tmp = map_cells_knn( topk=topk, qdat, cl.dat=index.bs$cl.dat, train.index=index.bs, method="cor", batch.size=blocksize, mc.cores=mc.cores )
         tmp.df = tmp$map.freq
         print("check tmp and tmpdf")
         rownames(tmp.df) = tmp.df$sample_id

         # order by qdat
         out.df = tmp.df[colnames(qdat),] 
         # hierarchical level by cl
         # prob that reaches to the node in level lvl. (indepent assumption)
         out.df$prob = prev.df$prob * out.df$prob
      }
   }
   ancestor.df = as.data.frame(ancestor)
   out.df$cl = ancestor.df[match(out.df$cl, ancestor.df[, ncol(ancestor.df)]), (lvl+1)]
   #out.df$cl = ancestor[match(out.df$best.cl, ancestor[, ncol(ancestor)]), (lvl+1)]
   return(out.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
predict_HKNN_cl <- function(lvl, iter_i, marker_index, query, assigned.df, children, ancestor, nodestr="", topk) {

   # initialize output
   out.df = assigned.df 

   assigned = assigned.df$cl
   names(assigned) = rownames(assigned.df)

   groups = unique(assigned)
   # for each group
   for (gg in groups) {
      nn.lvl = lvl
      nn = paste0(nodestr, gg)
      nn.sample   = names(assigned)[which(assigned == gg)]
      nn.children = children[[nn]]
      if ( nn %in% names(marker_index) && 
           length(marker_index[[nn]]$marker_tree)>0 && 
           !is.na(marker_index[[nn]]$marker_tree) ) {
         if (is.null(children[[nn]]) || is.na(children[[nn]])) {       # no children
            nn.assigned.df = assigned.df[nn.sample, ] 
         } else {                             # more than 1 children
            nn.indices   = marker_index[[nn]]$index_tree[[iter_i]]
            nn.markers   = marker_index[[nn]]$marker_tree
       
            gdx = match(nn.markers, rownames(query))
            sdx = match(nn.sample, colnames(query))
            
            nn.query     = query[gdx, sdx]

            if (length(nn.sample)==1) {
               nn.query = matrix(nn.query, ncol=1, dimnames=list(nn.markers, nn.sample)) 
            }
            nn.assigned.df  = call_ANN_cl (iter_i, lvl=nn.lvl, nn.indices, nn.query, assigned.df[nn.sample,], ancestor, topk) 
            nn.nodestr  = paste0(nn, ":")
            nn.assigned.df  = predict_HKNN_cl (lvl=nn.lvl+1, iter_i, marker_index, query, nn.assigned.df, children, ancestor, nn.nodestr, topk) 
         }
         out.df[nn.sample,] = nn.assigned.df[nn.sample,]
      } else {
         out.df[nn.sample,] = assigned.df[nn.sample,]
      }
   }
   return(out.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
predict_HKNN_cl_bs <- function (query, train.list, marker_index, iter=100, mc.cores=5, method="cor", topk=1) 
{
   #require(doMC)
   #require(foreach)

   iter = min(iter, length(marker_index[[1]]$index_tree))
   mc.cores = min(mc.cores, iter)
   registerDoMC(cores=mc.cores)

   #library(data.table)
   #map.list <- list()
   map.list <- foreach(i=1:iter, .combine="c") %dopar% {
   #for(i in 1:iter) {#, .combine="c") %dopar% {
      cat('\r', 'iter', i, '/', iter)
      assigned.df = create_initial_df(colnames(query), "root") 
      assigned.df = predict_HKNN_cl( lvl=1, iter_i=i, marker_index,  
                              query, 
                              assigned.df, 
                              children = train.list$children, 
                              ancestor = train.list$ancestor, 
                              nodestr = "", topk )
      assigned.df$cl = replace_subclass_by_cl (assigned.df$cl, train.list$cl.df)
      assigned.df$cl = replace_class_by_cl (assigned.df$cl, train.list$cl.df)

      if ("ancestor_markers" %in% names(train.list) && train.list$nlvl>2) {
         assigned.df$path.cor = cor_ancestor_markers(query, assigned.df, train.list$cl.dat, train.list$ancestor_markers)
      } else {
         assigned.df$path.cor = NA
      }
      assigned.df = add_labels (assigned.df, train.list$cl.df)
      list(assigned.df)
   #   map.list[[i]] =list(assigned.df)
   }
   map.df = rbindlist(map.list)
   map.df = map.df %>% group_by(sample_id, cl) %>% summarize(freq=n(),dist = mean(dist), path.cor=mean(path.cor))
   map.df$freq = map.df$freq/iter
   best.map.df = map.df %>% group_by(sample_id) %>% summarize(best.cl= cl[which.max(freq)],prob=max(freq), avg.dist = dist[which.max(freq)], avg.path.cor=path.cor[which.max(freq)])

   if(method=="cor"){
      best.map.df = best.map.df%>% mutate(avg.cor = 1 - avg.dist^2/2)
   }

   best.cl = best.map.df$best.cl
   names(best.cl) = best.map.df$sample_id
   avg.cor = best.map.df$avg.cor
   names(avg.cor) = best.map.df$sample_id
   best.map.df$cor.zscore = z_score(list(best.cl), avg.cor)

   # calculate zscore based on reference data
   if ("z_mean_sd" %in% names(train.list)) {
      best.map.df = cal_prob_ref_zscore (best.map.df, train.list$z_mean_sd)
   }

   return(list(map.freq=map.df, best.map.df = best.map.df))
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
cal_prob_ref_zscore <- function(best.map.df, z_mean_sd) 
{
    best.map.df = as.data.frame(best.map.df)
    cl = best.map.df$best.cl
    idx = match(cl, rownames(z_mean_sd))
    for (str in c("avg.path.cor", "avg.cor")) {
       mycor = best.map.df[, str]
       if (str %in% colnames(best.map.df)) {
          myz = (mycor-z_mean_sd[idx, paste0(str,".mean")])/z_mean_sd[idx,paste0(str, ".sd")]
          myp = pnorm(myz, lower.tail=T)
          best.map.df[, paste0(gsub("avg.", "", str), ".ref.zcore")] = myz
          best.map.df[, paste0(str, ".ref.prob")] = myp
       }
    }
    return(best.map.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
predict_fuzzy_HKNN_cl <- function(lvl, query, membership.df, children, ancestor, nodestr="") 
{

if (0) {
   Nlevel = ncol(ancestor) 
   Ng = rep(0, Nlevel)
   for (i in 2:Nlevel) Ng[i] = length(unique(ancestor[, i]))
   cumNg = cumsum(Ng)
   names(cumNg) = colnames(ancestor)

   Nsample = ncol(query)
   membership.df = data.frame( sample_id, cor, prob )
}   

   print(paste("==== ", nodestr))
   node.str = nodestr 
   groups = unique(ancestor[, lvl])
   for (gg in groups) {
      gg.lvl = lvl
      print(gg)
      gg.str = paste0(nodestr, gg)
      gg.sibling = unique(ancestor[, gg.lvl+1]) 
      if (gg.str %in% names(marker_index) && !is.na(marker_index[[gg.str]]$marker_tree)) {
         if (is.null(children[[gg.str]]) || is.na(children[[gg.str]])) {       # no children
            gg.membership.df = membership.df[, gg.sibling] 
         } else {                             # more than 1 children
            gg.indices   = marker_index[[gg.str]]$index_tree
            gg.markers   = marker_index[[gg.str]]$marker_tree
            gg.query     = query[gg.markers, ]
         }

         gg.membership.df  = call_membership_cl (lvl=gg.lvl, gg.indices, gg.query, membership.df, ancestor, topk=3) 
         gg.nodestr  = paste0(gg.str, ":")
         gg.membership.df  = predict_fuzzy_HKNN_cl (lvl=gg.lvl+1, query, membership.df, children, ancestor, gg.nodestr) 

         out.df[, gg.sibling] = gg.membership.df[, gg.sibling]
      } else {
         out.df[, gg.sibling] = gg.membership.df[, gg.sibling]
      }
   }
   return(out.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
run_predict_fuzzy_HKNN_cl <- function (query, train.list) {

   #initialize output as all 'node' and run
   membership.df = create_initial_membership(colnames(query), train.list$ancestor) 

   assigned.df = predict_fuzzy_HKNN_cl( lvl=1,
                              query, 
                              membership.df, 
                              children = train.list$children, 
                              ancestor = train.list$ancestor, 
                              nodestr="" )

   #assigned.df$cl = replace_subclass_by_cl (assigned.df$cl, train.list$cl.df)

   if ("ancestor_marker" %in% names(train.list)) {
      assigned.df$cor = cor_ancestor_markers(query, assigned.df, train.list$cl.dat, train.list$ancestor_markers)
   }

   assigned.df = add_labels (assigned.df, train.list$cl.df)
   return(assigned.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
mapping_summary <- function(qdat.mapped, FN=NA) {

   Nblocks = nrow(qdat.mapped)
   map.freq.all = c()
   for (i in 1:Nblocks) map.freq.all = rbind(map.freq.all, data.frame(qdat.mapped[i,1]))

   best.map.df.all = c()
   for (i in 1:Nblocks) best.map.df.all = rbind(best.map.df.all, data.frame(qdat.mapped[i,2]))

   summary = list()
   summary$map.freq = map.freq.all
   summary$best.map.df = best.map.df.all

   return(summary)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
mapping_by_block_serial <- function ( query.dat, train.list, blocksize=10000, mc.cores=7, method="cor", iter, topk=1, flag.fuzzy=FALSE ) 
{
 
   ##require(doMC)
   #require(foreach)
   mc.cores=mc.cores
   registerDoMC(cores=mc.cores)

   print(paste0(" getting marker_index from", train.list$MI_FN))
   print("....")
   load(train.list$MI_FN)
   print(".... marker_index ready!")
         
   qdat = query.dat[intersect(rownames(query.dat), train.list$all.markers), ]
   Nsample = ncol(qdat)
   nblk = ceiling(Nsample/blocksize)
   assigned <- c()
   for (i in 1:nblk) {
      istart = blocksize*(i-1) + 1
      istop  = blocksize*i
      if (istop > Nsample) istop=Nsample
      print(paste0("block ", i, "/", nblk ))
      idat = as.matrix(qdat[, istart:istop])
         
      iassigned = predict_HKNN_cl_bs(idat, train.list, marker_index, iter=iter, mc.cores=mc.cores, method=method, topk)
      assigned = rbind(assigned, iassigned)
   }
   return(assigned)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
mapping_by_block_parallel <- function ( query, train.list, blocksize=10000, iter=100, mc.cores=7, method="cor", topk=1,  flag.fuzzy=FALSE, tmpdir="./mapping_tmp" ) 
{

   print(paste0(" getting marker_index from", train.list$MI_FN))
   print("....")
   load(train.list$MI_FN)
   print(".... marker_index ready!")

   #require(doMC)
  # #require(foreach)
   mc.cores=mc.cores
   registerDoMC(cores=mc.cores)
         
   qdat = query.dat[intersect(rownames(query.dat), train.list$all.markers), ]
   Nsample = ncol(qdat)
   nblk = ceiling(Nsample/blocksize)
   #print(mc.cores)
   #print(Nsample)
   #print(nblk)


   #assigned <- c()
   assigned <- foreach (i=1:nblk, .combine='rbind') %dopar% {
      istart = blocksize*(i-1) + 1
      istop  = blocksize*i
      if (istop > Nsample) istop=Nsample
      print(paste("block", i, istart, istop))
      idat = as.matrix(qdat[, istart:istop])
         
      if (flag.fuzzy) {
         iassigned = predict_fuzzy_HKNN_cl_bs(idat, train.list, marker_index, mc.cores=mc.cores, topk)
      } else {
         iassigned = predict_HKNN_cl_bs(idat, train.list, marker_index, iter=iter, mc.cores=mc.cores, method=method, topk)
      }
      iassigned
   }
   return(assigned)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
mapping_by_block_parallel_tmp <- function ( query, train.list, blocksize=10000, iter=100, mc.cores=7, method="cor",
                                        flag.fuzzy=FALSE, topk=1, tmpdir=paste0("./mapping_tmp_", Sys.Date()) ) 
{
   #require(doMC)
  # #require(foreach)
   mc.cores=mc.cores
   registerDoMC(cores=mc.cores)

   print(paste0(" getting marker_index from", train.list$MI_FN))
   print("....")
   load(train.list$MI_FN)
   print(".... marker_index ready!")

   if (file.exists(tmpdir)) system(paste("rm -r", tmpdir)) 
   system(paste("mkdir", tmpdir))

   # run by blocks in parallel
   qdat = query.dat[intersect(rownames(query.dat), train.list$all.markers), ]
   Nsample = ncol(qdat)
   nblk = ceiling(Nsample/blocksize)
   foreach (i=1:nblk, .combine="c") %dopar% {
      print(i)
      istart = blocksize*(i-1) + 1
      istop  = blocksize*i
      if (istop > Nsample) istop=Nsample
      print(paste("block", i, istart, istop))
      idat = as.matrix(qdat[, istart:istop])

      iassigned = predict_HKNN_cl_bs(idat, train.list, marker_index, iter=iter, mc.cores=mc.cores, method=method, topk=topk)
      save(iassigned, file=file.path(tmpdir, paste0("iassigned_", istart, "_", istop, ".rda")))
   }

   # put together
   assigned = c()
   for (i in 1:nblk) {
      istart = blocksize*(i-1) + 1
      istop  = blocksize*i
      if (istop > Nsample) istop=Nsample
      load(file.path(tmpdir, paste0("iassigned_", istart, "_", istop, ".rda")))
      assigned = rbind(assigned, iassigned)
   }
   if (exists("assigned")) system(paste("rm -r", tmpdir))
   return(assigned)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
call_ANN <- function(index.bs, qdat, prev) {

   if (length(index.bs)==1 && is.na(index.bs)) { #no index
      out = prev
   } else {
      if (ncol(index.bs[[1]]$cl.dat)==0) { # no index
         out = prev
      } else {
         tmp = map_cells_knn_bs( qdat, train.index.bs=index.bs, method="cor" )
         tmp.cl = tmp$best.map.df$best.cl
         names(tmp.cl) = tmp$best.map.df$sample_id
         out = tmp.cl[colnames(qdat)] # order by qdat
      }
   }
   return(out)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
look_up_ancestor <- function(cl.df, nlvl=3) {


   #library(dplyr)

   if (nlvl == 3) {
      if ("class_label" %in% colnames(cl.df) && "subclass_label" %in% colnames(cl.df)) {
         ancestor = cl.df %>% select(class_label, subclass_label, cl) %>% 
                        mutate(level1=class_label, level2=subclass_label, level3=cl) %>% 
                        select (level1, level2, level3)
      }
      ancestor$level0 = rep("root", nrow(ancestor))
      out = ancestor[, c(4,1,2,3)]
   }
   if (nlvl == 2) {
      if ("class_label" %in% colnames(cl.df)) {
         ancestor = cl.df %>% select(neighborhood, cl) %>% 
                        mutate(level1=neighborhood, level2=cl) %>% 
                        select (level1, level2)
      }
      ancestor$level0 = rep("root", nrow(ancestor))
      out = ancestor[, c(3,1,2)]
   }
   if (nlvl == 1) {
      ancestor = cl.df %>% select(cl) %>% 
                        mutate(level1=cl) %>% 
                        select (level1)
      ancestor$level0 = rep("root", nrow(ancestor))
      out = ancestor[, c(2,1)]
   }

   ancestor = as.data.frame(ancestor)
   print("ancestor  LUT is set")
   return(out) 
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
look_up_children <- function(cl.df, nlvl=3) {

   children <- list()

   if (nlvl == 3) {
     if ("class_label" %in% colnames(cl.df) && "subclass_label" %in% colnames(cl.df)) {
       # root / class_label / subclass_label / cluster 
       nodestr = "root" 
       classes = unique(cl.df$class_label)
       children[[nodestr]] = classes
       for (nn in classes) {
         nn.str = paste0(nodestr, ":", nn)
         #print(nn.str)
         nn.cl.df = cl.df %>% filter(class_label == nn)
         nn.subclasses = unique(nn.cl.df$subclass_label)
         children[[nn.str]] = nn.subclasses
   
         for (ss in nn.subclasses) {
            ss.str = paste0(nn.str, ":", ss)
            #print(ss.str)
            ss.cl.df = nn.cl.df %>% filter(subclass_label == ss)
            ss.cls = unique(nn.cl.df$cl)
            children[[ss.str]] =ss.cls
         
            for (cc in ss.cls) {
               cc.str = paste0(ss.str, ":", cc)
               children[[cc.str]] = NA
            }
         }
       }
     }
   }

   if (nlvl == 2) {
      if ("class_label" %in% colnames(cl.df)) {
         # root / class(group) / cluster(cl)
         nodestr = "root" 
         groups = unique(cl.df$class_label)
         children[[nodestr]] = groups

         for (ss in groups) {
            ss.str = paste0(nodestr, ":", ss)
            #print(ss.str)
            ss.cl.df = cl.df %>% filter(class_label == ss)
            ss.cls = unique(ss.cl.df$cl)
            children[[ss.str]] =ss.cls
      
            for (cc in ss.cls) {
               cc.str = paste0(ss.str, ":", cc)
               children[[cc.str]] = NA
            }
         }
      }
   }

   if (nlvl == 1) {
      # root / cluster(cl)
      nodestr = "root" 
      cls = unique(cl.df$cl)
      children[[nodestr]] = cls

      for (cc in cls) {
         cc.str = paste0(nodestr, ":", cc)
         children[[cc.str]] = NA
      }
   }

   print("children LUT is set")
   return(children) 
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_marker_index_cl <- function( nn.str, lvl, nlevel, pre.marker_index=NA, query.genes, train.dat,
                                   cl.df, marker_index, mc.cores, div_thr=5, de.dir, subsample_pct, top.n.genes=15, 
                                   n.group.genes=3000, select.markers, cl.bin, outdir )
{

   nn = get_last_field(nn.str, ":")
   nn.cl.df = cl.df[cl.df[,lvl] == nn, ]
   if (lvl==nlevel) {
      nn.cl.group = nn.cl.df %>% select(lvl+1) 
      colnames(nn.cl.group) = c("cl")
      nn.cl.group = nn.cl.group %>% mutate(group = cl)
   } else {
      nn.cl.group = nn.cl.df %>% select(nlevel+1, lvl+1) 
      colnames(nn.cl.group) = c("cl", "group")
   }
   nn.group = unique(nn.cl.group$group)
   nn.cl = unique(nn.cl.group$cl)
   print(paste0("@", lvl, "/", nlevel, "   ngroup=",length(nn.group), ", ncluster=", length(nn.cl)))

   if (length(nn.group)==1  && length(nn.cl)==1) {
      marker_index[[nn.str]] = list(index_tree=NA, marker_tree=NA)
   } else {
      nn.tmp = gsub("/", "__", gsub(" ", "+", gsub("  ", "++", nn)))
      nn.markers.FN = paste0(outdir, "/marker.", nn.tmp,".rda")
      nn.markers.level.FN = paste0(outdir, "/marker.", lvl, ".", nn.tmp,".rda")
      if (file.exists(nn.markers.FN)) {
         if (file.info(nn.markers.FN)$size > 0) {
            load(nn.markers.FN)
            save(nn.markers, file=nn.markers.level.FN)
         } else nn.markers = NA
      } else {
         system(paste("cat /dev/null >", nn.markers.FN))
         if (is.na(pre.marker_index)) {
            if (length(nn.cl) == length(nn.group)) {
               nn.markers = select_markers_ds(de.dir, cl.bin, select.cl=nn.cl, top.n=top.n.genes)

            } else if (length(nn.group) > div_thr) {
               gc()
               nn.markers = select_markers_groups(de.dir, nn.cl.group, genes=select.markers,
                                                  cl.bin, select.groups=nn.group, 
                                                  n.markers=top.n.genes)#, byg=10, mc.cores=mc.cores)
            } else {
               nn.markers = select_markers_ds(de.dir, cl.bin, select.cl=nn.cl, top.n=top.n.genes)
            }
         } else {
            nn.markers = pre.marker_index[[nn.str]]$marker_tree
         }
         if (!is.na(query.genes)) nn.markers = intersect(nn.markers, query.genes )
         save(nn.markers, file=nn.markers.FN)
         save(nn.markers, file=nn.markers.level.FN)
      }
      nn.markers = intersect(rownames(train.dat), nn.markers)
      nn.cl.idx  = match(intersect(colnames(train.dat), nn.cl), colnames(train.dat))
   
      nn.dat = train.dat[nn.markers, nn.cl.idx]
      nn.index = build_train_index_bs( nn.dat, method="cor", 
                                       fn=paste0(outdir, "/index.", nn.tmp), 
                                       sample.markers.prop=subsample_pct, mc.cores=mc.cores )

      marker_index[[nn.str]] = list(index_tree=nn.index, marker_tree=nn.markers)
   }
   return(marker_index)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_marker_index_node_cl <- function( lvl, nlevel, nodestr, cl.df, marker_index, children, 
                                        pre.marker_index=NA, query.genes, train.dat, mc.cores=10, div_thr=5, 
                                        de.dir, subsample_pct, top.n.genes=15, n.group.genes=3000,
                                        select.markers, cl.bin, outdir)
{
   if (lvl==1) groups = "root"
   else groups = children[[substr(nodestr, 1, nchar(nodestr)-1)]]
   
   #print(groups)
   for (gg in groups) {
      nn = paste0(nodestr, gg)
      nn.lvl = lvl
      if (!(length(children[[nn]])==1 && is.na(children[[nn]]))) {
         print(paste("group :", nn))
         marker_index = build_marker_index_cl ( nn.str=nn, nn.lvl, nlevel, pre.marker_index, query.genes, train.dat, cl.df, marker_index, 
                                                mc.cores, div_thr=div_thr, de.dir, subsample_pct, top.n.genes, n.group.genes,
                                                select.markers, cl.bin, outdir)
         nn.nodestr  = paste0(nn, ":")
         marker_index = build_marker_index_node_cl ( lvl=nn.lvl+1, nlevel, nodestr=nn.nodestr, cl.df, marker_index, 
                                                     children, pre.marker_index, query.genes, train.dat, mc.cores, div_thr=div_thr, de.dir, 
                                                     subsample_pct, top.n.genes, n.group.genes, select.markers, cl.bin, outdir)
      }
   }
   return(marker_index)
}   

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_marker_index_tree_cl <- function (train.list, pre.train.list=NA, query.genes=NA, outdir, 
                                        div_thr=3, subsample_pct=0.9, top.n.genes=15, 
                                        n.group.genes=3000, mc.cores=10) 

{

   # train.list$nlvl=4 : root / class / subclass / cluster 
   # train.list$nlvl=3 : root / subclass / cluster 
   # train.list$nlvl=2 : root / cluster 
   nlevel   = train.list$nlvl-1
   children = look_up_children(train.list$cl.df, nlevel)
   ancestor = look_up_ancestor(train.list$cl.df, nlevel)

   de.dir  = train.list$dsFN
   cl.bin = train.list$cl.bin
   if (length(query.genes)==1 && is.na(query.genes)) {
      select.markers = train.list$select.markers 
   } else {
      select.markers = intersect(query.genes, train.list$select.markers)
   }
   train.dat      = train.list$cl.dat 
   cl.df.clean    = ancestor
   # create folder for parquet files
   dir.create(outdir)
   save(children, ancestor, file=file.path(outdir, "Train.children.ancestor.rda"))

   # load
   if (is.list(pre.train.list)) {
      tmp = load(pre.train.list$MI_FN)
      pre.marker_index = marker_index
      rm(marker_index)
   } else {
      pre.marker_index=NA
   }
   marker_index <- list()
   marker_index = build_marker_index_node_cl (lvl=1, nlevel, nodestr="", cl.df=cl.df.clean,
                                              marker_index, children, pre.marker_index, query.genes, train.dat,
                                              mc.cores=mc.cores, div_thr=div_thr, de.dir, 
                                              subsample_pct, top.n.genes, n.group.genes, 
                                              select.markers, cl.bin, outdir)
   all.markers = setdiff(unique(unlist(lapply(marker_index, function(x){x$marker_tree}))), NA)

   train.list$all.markers = all.markers
   train.list$children    = children
   train.list$ancestor    = ancestor
   train.list$MI_FN       = gsub("train.list.", "", train.list$TaxFN)
   save(marker_index, file=train.list$MI_FN)

   train.list$ancestor_markers = select_ancestor_markers( train.list$nlvl,
                                                          marker_index,
                                                          ancestor)

   return(train.list)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
prepare_train_mmean_list <- function( anal_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/forebrain_new",
                                modality = "10X_cells_v3" ) {

   #library(dplyr)

   # comb.dat  & imputed data
   tmp = load(file.path(anal_dir, "comb.dat.rda"))
   tmp = load(file.path(anal_dir, "impute.dat.refine.rda"))
   tmp = load(file.path(anal_dir, "cl.clean.rda")) # cl.df.clean & cl.clean
   genename.sel = rownames(impute.dat)

   # for each dat.list ("Smartseq_cells", "10X_cells_v3", "10X_cells")
   #require(doMC)
   #require(foreach)
   mc.cores=10
   registerDoMC(cores=mc.cores)
   if (is.na(modality)) mynames = c(names(comb.dat$dat.list), "imputed")
   else mynames = c(modality)
   foreach (i=1:length(mynames)) %dopar% {
      dd=mynames[i]; print(dd)
      if (dd == "imputed") {
         samplename.sel = intersect(colnames(impute.dat), names(cl.clean))
         mydat     = impute.dat[, samplename.sel]
      } else {
         dd.dat = comb.dat$dat.list[[dd]]
         gdx = match(genename.sel, dd.dat$row_id)
         samplename.sel = intersect(grep(dd, names(cl.clean), value=T), dd.dat$col_id)
         sdx = match(samplename.sel, dd.dat$col_id)
         mydat = dd.dat$fbm[gdx, sdx]
      }
      idx = match(cl.clean[samplename.sel], cl.df.clean$cl)
      myanno.class = cl.df.clean[idx,  "class_label"]
      myanno.subclass = cl.df.clean[idx,  "subclass_label"]
      myanno.cl = cl.df.clean[idx,  "cl"]
      
      my.mean.class = t(apply(mydat, 1, 
                           function(x) {y=aggregate(x, list(myanno.class), mean); 
                                        yy=y$x ; names(yy)=y$Group.1 ; return(yy)}))
      my.mean.subclass = t(apply(mydat, 1, 
                           function(x) {y=aggregate(x, list(myanno.subclass), mean); 
                                        yy=y$x ; names(yy)=y$Group.1 ; return(yy)}))
      my.mean.cluster  = t(apply(mydat, 1, 
                           function(x) {y=aggregate(x, list(myanno.cl), mean); 
                                        yy=y$x ; names(yy)=y$Group.1 ; return(yy)}))
      train.mmean.dat <- list()
      train.mmean.dat[["cluster" ]] = my.mean.cluster
      train.mmean.dat[["subclass"]] = my.mean.subclass
      train.mmean.dat[["class"]]    = my.mean.class

      save(train.mmean.dat, file=paste0("train.mmean.list.", dd, ".clean.rda"))
      print(paste0("train.mmean.list.", dd, ".clean.rda"))
   }
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
replace_cl_by_label <- function (pred, cl.df) {
   mypred = pred
   cl.df = as.data.frame(cl.df)
   cl.number.idx = which(mypred %in% cl.df$cl)
   if (length(cl.number.idx) > 0) {
      mypred.num = mypred[cl.number.idx]
      mypred[cl.number.idx] = cl.df[match(mypred.num, cl.df$cl), "cluster_label"]
   }
   mypred
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
replace_label_by_cl <- function (pred, cl.df) {
   mypred = pred
   cl.df = as.data.frame(cl.df)
   cl.label.idx = which(mypred %in% cl.df$cluster_label)
   if (length(cl.label.idx) > 0) {
      mypred.lbl = mypred[cl.label.idx]
      mypred[cl.label.idx] = cl.df[match(mypred.lbl, cl.df$cluster_label), "cl"]
   }
   mypred
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
replace_subclass_by_cl <- function (pred, cl.df) {
   mypred = pred
   cl.df = as.data.frame(cl.df)
   if ("subclass_label" %in% colnames(cl.df)) {
      subclass.label.idx = which(mypred %in% cl.df$subclass_label)
      if (length(subclass.label.idx) > 0) {
         mypred.lbl = mypred[subclass.label.idx]
         mypred[subclass.label.idx] = cl.df[match(mypred.lbl, cl.df$subclass_label), "cl"]
      }
   }
   mypred
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
replace_class_by_cl <- function (pred, cl.df) {
   mypred = pred
   cl.df = as.data.frame(cl.df)
   if ("class_label" %in% colnames(cl.df)) {
      class.label.idx = which(mypred %in% cl.df$class_label)
      if (length(class.label.idx) > 0) {
         mypred.lbl = mypred[class.label.idx]
         mypred[class.label.idx] = cl.df[match(mypred.lbl, cl.df$class_label), "cl"]
      }
   }
   mypred
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
add_labels <- function(in.df, cl.df) {
   pred.df = as.data.frame(in.df)
   cl.df = as.data.frame(cl.df)
   idx = match(pred.df$cl, cl.df$cl)
   if ("cluster_label" %in% colnames(cl.df)) pred.df$cluster_label = cl.df[idx, "cluster_label"]
   if ("subclass_label" %in% colnames(cl.df)) pred.df$subclass_label = cl.df[idx, "subclass_label"]
   if ("neighborhood" %in% colnames(cl.df)) pred.df$neighborhood = cl.df[idx, "neighborhood"]
   if ("class_label" %in% colnames(cl.df)) pred.df$class_label = cl.df[idx, "class_label"]
   return(pred.df)
}

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
get_last_field <- function(x, key) {
   unlist(lapply(strsplit(x, key), function(y) {n=length(y); y[n]}))
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
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
#' @param dat  to_be_added
#' @param ref.dat  to_be_added
#' @param k  to_be_added
#' @param method  to_be_added
#' @param dim  to_be_added
#'
#' @return to_be_added
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
      #library(RANN)
      knn.result = RANN::nn2(ref.dat, dat, k=k)
    }
    else if(method %in% c("Annoy.Euclidean", "Annoy.Cosine","cor")){
      #library(BiocNeighbors)      
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
#' @param dat  to_be_added
#' @param ref.dat  to_be_added
#' @param k  to_be_added
#' @param method  to_be_added
#' @param dim  to_be_added
#' @param batch.size  to_be_added
#' @param mc.cores  to_be_added
#'
#' @return to_be_added
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
#' @param x  to_be_added
#' @param batch.size  to_be_added
#' @param FUN  to_be_added
#' @param mc.cores  to_be_added
#' @param .combine  to_be_added
#' @param ... 
#'
#' @return to_be_added
#'
#' @keywords internal
batch_process <- function(x, batch.size, FUN, mc.cores=1, .combine="c",...)
  {
    #require(foreach)
    #require(doMC)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_index <- function(cl.dat, method= c("Annoy.Cosine","cor","Annoy.Euclidean"),fn=tempfile(fileext=".idx"))
  {
    #library(BiocNeighbors)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_index_bs <- function(cl.dat, method= c("Annoy.Cosine","cor","Annoy.Euclidean"),sample.markers.prop=0.8, iter=100, mc.cores=10,fn=tempfile(fileext=".idx"))
  {
    #library(BiocNeighbors)
    #require(doMC)
    #require(foreach)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
map_cells_knn_bs <- function(topk=1, test.dat, iter=100,cl.dat=NULL,train.index.bs=NULL, method = c("Annoy.Cosine","cor"), mc.cores=20, ...)
  {
    #require(doMC)
    #require(foreach)
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
    #library(data.table)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
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
    #require(doMC)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
update_gene_score_ds <- function(gene.score, ds, to.remove, cl.bin, de=NULL, max.num=1000,mc.cores=20)
  {
    #require(doMC)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
check_pairs_ds <- function(de.dir, to.add, genes,cl.bin, de=NULL, mc.cores=10,max.num=1000)
  {
    #require(doMC)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
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
    #library(parallel)    
    #require(doMC)
    #require(foreach)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
select_markers_pair_group_top_ds <- function(g1,g2,ds, genes, cl.bin, select.sign="up", n.markers=20,mc.cores=1, ...)
{
  #require(matrixStats)
  #require(data.table)
  #require(arrow)
  #require(dplyr)

  #require(doMC)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
select_pos_markers_ds <- function(de.dir, cl, select.cl, genes, cl.bin, n.markers=1,  mc.cores=1,out.dir="cl.markers",...)
  {
    #library(parallel)    
    #require(doMC)
    #require(foreach)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
select_top_pos_markers_ds<- function(ds, cl, select.cl, genes, cl.bin, n.markers=3, mc.cores=10,...)
  {
    #library(parallel)    
    #require(doMC)
    #require(foreach)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
select_markers_groups_top_ds <- function(ds, cl.group, select.groups=names(cl.group), n.markers=3,mc.cores=1,...)
  {

    #library(parallel)
    #library(parallel)    
    #require(doMC)
    #require(foreach)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
select_markers_groups <- function(de.dir, cl.group, genes, cl.bin, select.groups= unique(cl.group$group), n.markers=20,mc.cores=1,byg=10, ...)
  {
    ds = open_dataset(de.dir)
    cl.group$cl = as.character(cl.group$cl)
    group_pair=create_pairs(unique(cl.group$group))
    #library(parallel)
    #require(doMC)
    #require(foreach)
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
#' @param in.df to_be_added
#' @param cl.df to_be_added
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



#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal

prep_parquet_de_all_pairs <- function (norm.dat, cl, cl.bin = NULL, cl.bin.size = 100, de.param = de_param(),
    method = "fast_limma", mc.cores = 1, pairs.fn = "pairs.parquet", cl.bin.fn="cl.bin.rda", ...)
{
    cn <- as.character(sort(unique(cl)))
    pairs = create_pairs(cn)
    pairs = as.data.frame(pairs)
    pairs$pair = row.names(pairs)
    pairs$pair_id = 1:nrow(pairs)
    if (is.null(cl.bin)) {
       cl.bin.size = min(100, length(cn)/mc.cores)
       cl.bin = data.frame(cl = cn, bin = ceiling((1:length(cn)/cl.bin.size)))
    }
    #library(arrow)
    write_parquet(pairs, sink = pairs.fn)
    save(cl.bin, file=cl.bin.fn) 
    de.result = prep_parquet_de_selected_pairs(norm.dat, cl = cl, pairs = pairs,
        cl.bin = cl.bin, de.param = de.param, method = method,
        mc.cores = mc.cores, ...)
    #de.result = de_selected_pairs(norm.dat, cl = cl, pairs = pairs,
    #    de.param = de.param, method = method, mc.cores = mc.cores, ...)
    return(de.result$de.genes)
}

#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @keywords internal
prep_parquet_de_selected_pairs <- function (norm.dat, cl, pairs, cl.bin = NULL, cl.size = NULL,
    de.param = de_parm(), method = "fast_limma", cl.means = NULL,
    cl.present = NULL, cl.sqr.means = NULL, use.voom = FALSE,
    counts = NULL, mc.cores = 1, out.dir = NULL, summary.dir = NULL,
    top.n = 500, overwrite = FALSE, return.df = FALSE, return.summary = !is.null(summary.dir))
{
    #library(arrow)
    method <- match.arg(method, choices = c("fast_limma", "limma",
        "chisq", "t.test"))
    #require(parallel)
    if (use.voom & is.null(counts)) {
        stop("The use.voom = TRUE parameter requires a raw count matrix via the counts parameter.")
    }
    if (is.null(cl.size)) {
        cl.size <- table(cl)
        cl.size = setNames(as.integer(cl.size), names(cl.size))
    }
    pairs.fn = NULL
    if (length(pairs) == 1) {
        pairs.fn = pairs
        pairs = open_dataset(pairs.fn)
    }
    else {
        pairs = as.data.frame(pairs)
        if (is.null(pairs$pair)) {
            pairs$pair = row.names(pairs)
        }
        if (is.null(pairs$pair_id)) {
            pairs$pair_id = 1:nrow(pairs)
        }
    }
    select.cl <- unique(c(pairs %>% pull(P1), pairs %>% pull(P2)))
    select.cl <- intersect(select.cl, names(cl.size)[cl.size >=
        de.param$min.cells])
    cl <- cl[cl %in% select.cl]
    if (is.factor(cl)) {
        cl = droplevels(cl)
    }
    if (is.null(pairs$bin.x)) {
        if (is.null(cl.bin)) {
            cl.bin.size = min(100, length(select.cl)/mc.cores)
            cl.bin = data.frame(cl = select.cl, bin = ceiling((1:length(select.cl)/cl.bin.size)))
        }
        pairs = pairs %>% left_join(cl.bin, by = c(P1 = "cl")) %>%
            left_join(cl.bin, by = c(P2 = "cl"))
    }
    pairs = pairs %>% filter(P1 %in% select.cl & P2 %in% select.cl) %>%
        collect()
    cl.size = cl.size[select.cl]
    if (is.null(cl.means)) {
        cl.means <- as.data.frame(get_cl_means(norm.dat, cl))
    }
    else {
        cl.means <- as.data.frame(cl.means)
    }
    if (is.null(cl.present)) {
        cl.present <- as.data.frame(get_cl_present(norm.dat,
            cl, de.param$low.th))
    }
    else {
        cl.present <- as.data.frame(cl.present)
    }
    if (is.null(cl.sqr.means)) {
        cl.sqr.means <- as.data.frame(get_cl_sqr_means(norm.dat,
            cl))
    }
    else {
        cl.sqr.means <- as.data.frame(cl.sqr.means)
    }
    if (method == "limma") {
        #require("limma")
        norm.dat <- as.matrix(norm.dat[, names(cl)])
        cl <- setNames(as.factor(paste0("cl", cl)), names(cl))
        design <- model.matrix(~0 + cl)
        colnames(design) <- levels(as.factor(cl))
        if (use.voom & !is.null(counts)) {
            v <- limma::voom(counts = as.matrix(counts[row.names(norm.dat),
                names(cl)]), design = design)
            fit <- limma::lmFit(object = v, design = design)
        }
        else {
            fit <- limma::lmFit(object = norm.dat[, names(cl)],
                design = design)
        }
    }
    else if (method == "fast_limma") {
        fit = simple_lmFit(norm.dat, cl = cl, cl.means = cl.means,
            cl.sqr.means = cl.sqr.means)
    }
    else if (method == "t.test") {
        cl.vars <- as.data.frame(get_cl_vars(norm.dat, cl, cl.means = cl.means))
    }
    #require(doMC)
    #require(foreach)
    mc.cores = min(mc.cores, ceiling(nrow(pairs)/5000))
    registerDoMC(cores = mc.cores)
    de_combine <- function(result.1, result.2) {
        #library(data.table)
        de.genes = c(result.1$de.genes, result.2$de.genes)
        if (!is.null(result.1$de.summary)) {
            de.summary = rbindlist(result.1$de.summary, result.2$de.summary)
            return(list(de.genes = de.genes, de.summary = de.summary))
        }
        else {
            return(list(de.genes = de.genes))
        }
    }
    if (!is.null(out.dir)) {
        if (!dir.exists(out.dir)) {
            dir.create(out.dir)
        }
    }
    if (!is.null(summary.dir)) {
        if (!dir.exists(summary.dir)) {
            dir.create(summary.dir)
        }
    }
    all.bins = sort(unique(c(pairs$bin.x, pairs$bin.y)))
    de_list = foreach(bin1 = 1:length(all.bins), .combine = "c") %:%
        foreach(bin2 = bin1:length(all.bins), .combine = "c") %dopar%
       {
            x = all.bins[bin1]
            y = all.bins[bin2]
            if (!overwrite & !is.null(out.dir)) {
                if (file.exists(file.path(out.dir, paste0("bin.x=",
                  x), paste0("bin.y=", y)))) {
                  return(list(de.genes = NULL, de.summary = NULL))
                }
            }
            #library(dplyr)
            #library(arrow)
            #library(data.table)
            tmp.pairs = pairs %>% filter(bin.x == x & bin.y ==
                y | bin.x == y & bin.y == x)
            if (is.null(tmp.pairs) | nrow(tmp.pairs) == 0) {
                return(list(de.genes = NULL, de.summary = NULL))
            }
            de.genes = sapply(1:nrow(tmp.pairs), function(i) {
                pair = unlist(tmp.pairs[i, c("P1", "P2")])
                if (method == "limma") {
                  #require("limma")
                  df = de_pair_limma(pair = pair, cl.present = cl.present,
                    cl.means = cl.means, design = design, fit = fit)
                }
                else if (method == "fast_limma") {
                  df = de_pair_fast_limma(pair = pair, cl.present = cl.present,
                    cl.means = cl.means, fit = fit)
                }
                else if (method == "t.test") {
                  df = de_pair_t.test(pair = pair, cl.present = cl.present,
                    cl.means = cl.means, cl.vars = cl.vars, cl.size = cl.size)
                }
                else if (method == "chisq") {
                  df = de_pair_chisq(pair = pair, cl.present = cl.present,
                    cl.means = cl.means, cl.size = cl.size)
                }
                if (!is.null(de.param$min.cells)) {
                  cl.size1 <- cl.size[as.character(pair[1])]
                  cl.size2 <- cl.size[as.character(pair[2])]
                }
                else {
                 cl.size1 <- NULL
                 cl.size2 <- NULL
                }
                stats = de_stats_pair(df, de.param = de.param,
                  cl.size1, cl.size2, return.df = return.df)
            }, simplify = F)
            pair = tmp.pairs %>% pull(pair)
            names(de.genes) = pair
            if (return.summary) {
                de.summary = scrattch.bigcat::de_pair_summary(de.genes, cl.bin=cl.bin, out.dir = summary.dir,
                  pairs = tmp.pairs, return.df = is.null(summary.dir))
            }
            if (!is.null(out.dir)) {
                cat("Export", bin1, bin2, "\n")
                result = scrattch.bigcat::export_de_genes(de.genes, cl.means,
                  out.dir = out.dir, pairs = tmp.pairs,
                  mc.cores = 1, top.n = top.n)
                cat("Finish Export", x, y, "\n")
            }
            out = list(de.genes = de.genes, de.summary = de.summary)
          
            de.summary = NULL
            de.genes = NULL
            out
        }

    de.genes = do.call("c", de_list[names(de_list) == "de.genes"])
    names(de.genes) = gsub("de.genes.", "", names(de.genes))
    de.summary = do.call("c", de_list[names(de_list) == "de.summary"])
    if (is.null(de.summary)) {
        return(de.genes)
    }
    return(list(de.genes, de.summary))
}


