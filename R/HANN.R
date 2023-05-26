#' INFO -- PLEASE ADD --
#'
#' @param in.df to_be_added
#' @param cl.df to_be_added
#'
#' @return ???
#'
#' @export
run_mapping_on_taxonomy <- function(query.dat, 
                             Taxonomy="AIT12.0_mouse", prefix="", TaxFN=NA, 
                             prebuild=FALSE, newbuild=FALSE, 
                             mapping.method=c('flat', 'hierarchy'),
                             iter=100, mc.cores=7, blocksize=50000, dist.method="cor", topk=1, 
                             subsample_pct=0.9, top.n.genes=15, rm.clusters=NA, 
                             flag.serial=TRUE, flag.parallel.tmp=FALSE, flag.fuzzy=FALSE) 
{
   print("### Training Templates for Query Data")
   train.list = build_train_list_on_taxonomy ( TaxFN, Taxonomy, pre.train.list, 
                                               query.genes=rownames(query.dat),
                                               prefix=prefix, mapping.method=mapping.method, 
                                               prebuild=prebuild, newbuild=newbuild,
                                               mc.cores=mc.cores, subsample_pct=subsample_pct )

   
   print("### Mapping on Training Templates")
   query.dat.mapped = mapping_on_taxonomy (query.dat, train.list, blocksize=blocksize, mc.cores=mc.cores, 
                                method=dist.method, iter=iter, topk=topk, 
                                flag.serial=flag.serial, flag.parallel.tmp, flag.fuzzy=flag.fuzzy )

   print("### Adding cluster related info")
   query.dat.mapped$cl.df = train.list$cl.df

   return(query.dat.mapped)
}

#' Starting point for optimized tree mapping
#'
#' @param TaxFN to_be_added
#' @param Taxonomy to_be_added
#' @param pre.train.list to_be_added
#' @param query.genes to_be_added
#' @param prefix to_be_added
#' @param mapping.method to_be_added
#' @param prebuild to_be_added
#' @param newbuild to_be_added
#' @param mc.cores to_be_added
#' @param div_thr to_be_added
#' @param subsample_pct to_be_added
#' @param top.n.genes to_be_added
#' @param n.group.genes to_be_added
#' @param rm.cl to_be_added
#'
#' @return Mapping results
#'
#' @export
build_train_list_on_taxonomy <- function ( TaxFN=NA, Taxonomy, pre.train.list=NA, 
                                           query.genes=NA, prefix="", mapping.method=c('flat', 'hierarchy'),
                                           prebuild=FALSE, newbuild=FALSE,
                                           mc.cores=10, div_thr=3, subsample_pct= 0.9, top.n.genes=15, 
                                           n.group.genes=3000, rm.cl=c() )
{

   if (is.na(Taxonomy)) {
      if (is.na(TaxFN)) {
         print("error : Either Taxonomy or TaxFN needs to be specified")
         break
      } else {
         print(paste("loading taxonomy train templates", TaxFN))
         load(TaxFN)
      }
   } else {
      if (Taxonomy=="AIT20.0_macaque") { 
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/Macaquet_MTG"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT20.0_macaque"
      }
      if (Taxonomy=="AIT17.BG.1_mouse") { #WB_TH
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT17.BG.1_mouse"
      }
      if (Taxonomy=="AIT17.FB.1_mouse") { #WB_TH
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT17.FB.1_mouse"
         TaxDir = file.path(TaxDir , "Templates")
      }
      if (Taxonomy=="AIT17.TR_mouse") { #WB_TH
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT17.TR_mouse"
         TaxDir = file.path(TaxDir , "Templates")
      }
      if (Taxonomy %in% c("AIT17.2_mouse", "AIT17.TH.1_mouse")) { #WB_TH
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT17.2_mouse"
         TaxDir = file.path(TaxDir , "Templates")
      }
      if (Taxonomy=="AIT17.0_mouse_5200") { #WB_2022_Yao
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT17.0_mouse_5200"
         TaxDir = file.path(TaxDir , "Templates")
      }
      if (Taxonomy=="AIT17.0_mouse") { #WB_2022_Yao
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT17.0_mouse"
         TaxDir = file.path(TaxDir , "Templates")
      }
      if (Taxonomy=="AIT16.2_mouse") { #WB_TH
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT16.2_mouse"
         TaxDir = file.path(TaxDir , "Templates")
      }
      if (Taxonomy=="AIT16.1_mouse") { #WB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT16.1_mouse"
      }
      if (Taxonomy=="AIT16.0_mouse") { #WB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb/v3"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT16.0_mouse"
         TaxDir = file.path(TaxDir , "Templates")
      }
      if (Taxonomy=="AIT14.0_mouse") { #WB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb_nuclei"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT13.1_mouse"
      }
      if (Taxonomy %in% c("AIT13.2_mouse","AIT13.TH.1_mouse")) { #WB_TH
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT13.2_mouse"
      }
      if (Taxonomy=="AIT13.0_mouse") { #WB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb/v3"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT13.0_mouse"
      }
      if (Taxonomy=="AIT12.1_mouse") { #WB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb/v2"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT12.1_mouse"
      }
      if (Taxonomy=="AIT12.0_mouse") { #WB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/wb/v2"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT12.0_mouse"
      }
      if (Taxonomy %in% c("AIT11.0_mouse", "AIT9.4_mouse", "AIT9.BG.1_mouse")) { #BG
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/Mapping_PatchSeq/BasalGanglia/FB_Joint_Taxonomy"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT11.0_mouse"
      }
      if (Taxonomy=="AIT10.0_mouse") { #WB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/10X_analysis/wholebrain_v3/v1"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT10.0_mouse"
      }
      if (Taxonomy=="AIT9.0_mouse") { #FB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/forebrain_new"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT9.0_mouse"
         TaxDir = file.path(TaxDir , "Templates")
      }
      if (Taxonomy=="AIT9.LO10_mouse") { #FB
         TrainDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/forebrain_new"
         TaxDir   = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT9.0_mouse/LO10"
      }
      #TaxDir = file.path(TaxDir , "Templates")

      # set hierarchy level for each  Taxonomy
      if (mapping.method=="flat") {
         nlevel=2
      } else {
         nlevel=4
         if (Taxonomy=="AIT20.0_macaque") { 
            nlevel = 3
         }
      }
      if (is.na(subsample_pct)) nlevel_str = paste0(nlevel, "_noboot") 
      else nlevel_str = nlevel
   
      if (is.na(TaxFN)) TaxFN = file.path(TaxDir, paste0('train.list.nlevel', nlevel_str, '_marker_index.rda'))

      if (file.exists(TaxFN)) {
         print(paste("### ", Taxonomy, " Base Training Template (markers & indices) are being read ..."))
         load(TaxFN)
      } else {
         print("### Building Base Training Teamplates ...")
         pre.train.list = NA

         if (Taxonomy=="AIT20.0_macaque") train.list = build_train_list_20( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT17.BG.1_mouse") train.list = build_train_list_WB17_BG( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT17.FB.1_mouse") train.list = build_train_list_WB17_FB( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT17.TR_mouse") train.list = build_train_list_WB17_TR( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT17.2_mouse") train.list = build_train_list_TH17_2( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT17.0_mouse") train.list = build_train_list_WB17( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT17.0_mouse_5200") train.list = build_train_list_WB17_5200( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT14.0_mouse") train.list = build_train_list_WBnuc( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT13.2_mouse") train.list = build_train_list_TH13_2(
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir,
                                                        prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT16.2_mouse") train.list = build_train_list_TH16_2( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT16.1_mouse") train.list = build_train_list_WB16nuc( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT16.0_mouse") train.list = build_train_list_WB16( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT13.0_mouse") train.list = build_train_list_WB13( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT12.0_mouse") train.list = build_train_list_WBint( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT11.0_mouse" || Taxonomy=="AIT9.4_mouse") train.list = build_train_list_BG( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT10.0_mouse") train.list = build_train_list_WB( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT9.0_mouse")  train.list = build_train_list_FB( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)
         if (Taxonomy=="AIT9.LO10_mouse")  train.list = build_train_list_FBLO10( 
                                                        pre.train.list, query.genes=NA, TrainDir, TaxDir, 
							prefix="", nlevel=nlevel, TaxFN=TaxFN)

         require(doMC)
         require(foreach)
         registerDoMC(cores=mc.cores)

         MI_str = paste0('nlevel', nlevel_str, '_marker_index')
         train.list = build_marker_index_tree_cl ( train.list, pre.train.list=pre.train.list, query.genes=NA,
                                                   outdir=file.path(TaxDir, MI_str),
                                                   div_thr=3, subsample_pct=subsample_pct,
                                                   top.n.genes=15, n.group.genes=3000, mc.cores=mc.cores)
         save(train.list, file=train.list$TaxFN)
      }

      print("### Base Training Teamplates Are Ready.\n")

      gdx = match(train.list$all.markers, query.genes)
      if (any(is.na(gdx))) {
         print("### Alert! Some marker genes are missing in your data")
         train.list$TaxFN = gsub("train.list", paste0(prefix, "_train.list"), TaxFN) 
         if (prebuild && file.exists(train.list$TaxFN)) {
            print("### Marker gene selection is already done. let's load it")
            print(paste("### !!! You chose to use the exisiting marker_index for", prefix, "data."))
            print(TaxFN)
            print("### Enter 'c' to continue. Or rerun with prebuild=FALSE")
            tmp = load(train.list$TaxFN)
         } else {
            print(paste(sum(is.na(gdx)), 'genes are missing :'))
            print(head(train.list$all.markers[is.na(gdx)]))
            print("...")
            print(tail(train.list$all.markers[is.na(gdx)]))
            if (newbuild) {
               print("### Rebuilding the marker genes and indices using only availble genes.")
               pre.train.list = NA
               #prefix = paste0(prefix, "_rebuilt")
               train.list$cl.dat = train.list$cl.dat[train.list$all.markers,]
               select.markers = intersect(train.list$all.markers, query.genes)
               train.list$select.markers = select.markers
            } else {
               print("### Using the subset of existing marker genes and indices.")
               pre.train.list = train.list
               train.list$cl.dat = pre.train.list$cl.dat[train.list$all.markers,]
               select.markers = intersect(pre.train.list$all.markers, query.genes)
               train.list$select.markers = select.markers
            }
            
            MI_str = paste0(prefix, '_nlevel', nlevel, '_marker_index')

            require(doMC)
            require(foreach)
            registerDoMC(cores=mc.cores)
            train.list = build_marker_index_tree_cl ( train.list, pre.train.list=pre.train.list, query.genes,
                                                      outdir=file.path(TaxDir, MI_str),
                                                      div_thr=3, subsample_pct=subsample_pct,
                                                      top.n.genes=15, n.group.genes=3000, mc.cores=mc.cores)
            save(train.list, file=train.list$TaxFN)
         
         }
      } else {
         print("### All marker genes are in your data set.")
         print("### Precalculated marker genes and indices for hierarchical knn mapping will be  used.")
      }
   
      print("### Training Teamplates Are Ready for Query Data.\n")
   }
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
   #print("call_ANN_cl") ; browser()
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
   require(doMC)
   require(foreach)

   iter = min(iter, length(marker_index[[1]]$index_tree))
   mc.cores = min(mc.cores, iter)
   registerDoMC(cores=mc.cores)

   library(data.table)
   map.list <- foreach(i=1:iter, .combine="c") %dopar% {
      cat('\r', 'iter', i, '/', iter)
      assigned.df = create_initial_df(colnames(query), "root") 
      assigned.df = predict_HKNN_cl( lvl=1, iter_i=i, marker_index,  
                              query, 
                              assigned.df, 
                              children = train.list$children, 
                              ancestor = train.list$ancestor, 
                              nodestr = "", topk )
      assigned.df$cl = replace_subclass_by_cl (assigned.df$cl, train.list$cl.df)

      if ("ancestor_markers" %in% names(train.list) && train.list$nlvl>2) {
         assigned.df$path.cor = cor_ancestor_markers(query, assigned.df, train.list$cl.dat, train.list$ancestor_markers)
      } else {
         assigned.df$path.cor = NA
      }
      assigned.df = add_labels (assigned.df, train.list$cl.df)
      list(assigned.df)
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
mapping_on_taxonomy <- function ( query.dat, train.list, 
                               blocksize=5000, mc.cores=5, method='cor', iter=100, topk=1, 
                               flag.serial=TRUE, flag.parallel.tmp=FALSE, flag.fuzzy=FALSE )
{
   print(paste0("Mapping by block(", blocksize, ")...."))
   qdat = query.dat[intersect(rownames(query.dat), train.list$all.markers), ]

   if (flag.serial) {
      assigned = mapping_by_block_serial ( qdat, train.list, blocksize=blocksize, mc.cores=mc.cores, 
                                           method=method, iter=iter, topk=topk, flag.fuzzy=flag.fuzzy ) 
   } else {
      if (flag.parallel.tmp ){   
         assigned = mapping_by_block_parallel_tmp ( qdat, train.list, blocksize=blocksize, mc.cores=mc.cores, 
                                                    method=method, iter=iter, topk=topk, flag.fuzzy=flag.fuzzy ) 
      } else {
         assigned = mapping_by_block_parallel ( qdat, train.list, blocksize=blocksize, mc.cores=mc.cores, 
                                                method=method, iter=iter, topk=topk, flag.fuzzy=flag.fuzzy ) 
      }
   }

   if (length(class(assigned))==2) assigned = mapping_summary(assigned)
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
mapping_by_block_serial <- function ( query.dat, train.list, blocksize=10000, mc.cores=7, method="cor", iter, topk=1, flag.fuzzy=FALSE ) 
{
 

   require(doMC)
  # require(foreach)
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

   require(doMC)
  # require(foreach)
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
   require(doMC)
  # require(foreach)
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


   library(dplyr)
   if (nlvl == 3) {
      if ("neighborhood" %in% colnames(cl.df) && "subclass_label" %in% colnames(cl.df)) {
         ancestor = cl.df %>% select(neighborhood, subclass_label, cl) %>% 
                        mutate(level1=neighborhood, level2=subclass_label, level3=cl) %>% 
                        select (level1, level2, level3)
      }
      if ("Level1_label" %in% colnames(cl.df) && "Level2_label" %in% colnames(cl.df)) {
         ancestor = cl.df %>% select(Level1_label, Level2_label, cl) %>% 
                        mutate(level1=Level1_label, level2=Level2_label, level3=cl) %>% 
                        select (level1, level2, level3)
      }
      ancestor$level0 = rep("root", nrow(ancestor))
      out = ancestor[, c(4,1,2,3)]
   }
   if (nlvl == 2) {
      ancestor = cl.df %>% select(group, cl) %>% 
                        mutate(level1=group, level2=cl) %>% 
                        select (level1, level2)
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
     if ("neighborhood" %in% colnames(cl.df) && "subclass_label" %in% colnames(cl.df)) {
       # root / neighborhood / subclass / cluster 
       nodestr = "root" 
       neighborhoods = unique(cl.df$neighborhood)
       children[[nodestr]] = neighborhoods
       for (nn in neighborhoods) {
         nn.str = paste0(nodestr, ":", nn)
         #print(nn.str)
         nn.cl.df = cl.df %>% filter(neighborhood == nn)
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
     if ("Level1_label" %in% colnames(cl.df) && "Level2_label" %in% colnames(cl.df)) {
       # root / Level1_label / subclass / cluster 
       nodestr = "root" 
       Level1_labels = unique(cl.df$Level1_label)
       children[[nodestr]] = Level1_labels
       for (nn in Level1_labels) {
         nn.str = paste0(nodestr, ":", nn)
         #print(nn.str)
         nn.cl.df = cl.df %>% filter(Level1_label == nn)
         nn.Level2 = unique(nn.cl.df$Level2_label)
         children[[nn.str]] = nn.Level2
   
         for (ss in nn.Level2) {
            ss.str = paste0(nn.str, ":", ss)
            #print(ss.str)
            ss.cl.df = nn.cl.df %>% filter(Level2_label == ss)
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
      # root / subclass(group) / cluster(cl)
      nodestr = "root" 
      groups = unique(cl.df$group)
      children[[nodestr]] = groups

      for (ss in groups) {
         ss.str = paste0(nodestr, ":", ss)
         #print(ss.str)
         ss.cl.df = cl.df %>% filter(group == ss)
         ss.cls = unique(ss.cl.df$cl)
         children[[ss.str]] =ss.cls
      
         for (cc in ss.cls) {
            cc.str = paste0(ss.str, ":", cc)
            children[[cc.str]] = NA
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

   # train.list$nlvl=4 : root / neighborhood / subclass / cluster 
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

   library(dplyr)

   # comb.dat  & imputed data
   tmp = load(file.path(anal_dir, "comb.dat.rda"))
   tmp = load(file.path(anal_dir, "impute.dat.refine.rda"))
   tmp = load(file.path(anal_dir, "cl.clean.rda")) # cl.df.clean & cl.clean
   genename.sel = rownames(impute.dat)

   # for each dat.list ("Smartseq_cells", "10X_cells_v3", "10X_cells")
   require(doMC)
   require(foreach)
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
      myanno.neighborhood = cl.df.clean[idx,  "neighborhood"]
      myanno.subclass = cl.df.clean[idx,  "subclass_label"]
      myanno.cl = cl.df.clean[idx,  "cl"]
      
      my.mean.neighbor = t(apply(mydat, 1, 
                           function(x) {y=aggregate(x, list(myanno.neighborhood), mean); 
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
      train.mmean.dat[["neighbor"]] = my.mean.neighbor

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
   if ("Level2_label" %in% colnames(cl.df)) {
      Level2_label.idx = which(mypred %in% cl.df$Level2_label)
      if (length(Level2_label.idx) > 0) {
         mypred.lbl = mypred[Level2_label.idx]
         mypred[Level2_label.idx] = cl.df[match(mypred.lbl, cl.df$Level2_label), "cl"]
      }
   }
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
#' @param index.bs to_be_added
#' @param qdat to_be_added
#' @param prev to_be_added
#' @param children  to_be_added
#' @param gg  to_be_added
#'
#' @return ???
#'
#' @keywords internal
call_ANN_0 <- function(index.bs, qdat, prev, children, gg) {
   out = cl.clean[colnames(qdat)]
   if (gg == "node") out.nnn = as.character(cl.df.clean[match(out, cl.df.clean$cl), "neighborhood"])
   if (gg %in% neighborhoods) out.nnn = as.character(cl.df.clean[match(out, cl.df.clean$cl), "subclass_label"])
   if (gg %in% subclasses) out.nnn = as.character(cl.df.clean[match(out, cl.df.clean$cl), "cl"])
   names(out.nnn) = colnames(qdat)
   return(out.nnn)
}