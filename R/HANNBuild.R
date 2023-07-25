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
#' @import doMC
#' @import foreach
#'
#' @return Mapping results
#'
#' @export
build_train_list_on_taxonomy <- function ( TaxFN=NA, Taxonomy, 
                                           TaxHome='/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/',
                                           pre.train.list=NA, 
                                           query.genes=NA, 
                                           prefix="", 
                                           mapping.method=c('flat', 'one-step', 'hierarchy'),
                                           nlevel=2, 
                                           prebuild=FALSE, 
                                           newbuild=FALSE,
                                           mc.cores=10, 
                                           div_thr=3, 
                                           subsample_pct= 0.9, 
                                           top.n.genes=15, 
                                           n.group.genes=3000, 
                                           rm.cl=c()){

   ##
   if (is.na(Taxonomy)) {
      if (is.na(TaxFN)) {
         print("error : Either Taxonomy or TaxFN needs to be specified")
         break
      } else {
         print(paste("loading taxonomy train templates", TaxFN))
         load(TaxFN)
      }
   } 
   else {
      ##
      TrainDir = file.path(TaxHome, Taxonomy)
      TaxDir   = file.path(TaxHome, Taxonomy)

      # set hierarchy level for each  Taxonomy
      if (mapping.method %in% c('flat', 'one-step')) {
         nlevel=2
      }
      
      ##
      if(is.na(subsample_pct)) nlevel_str = paste0(nlevel, "_noboot") else nlevel_str = nlevel
   
      ##
      if(is.na(TaxFN)) TaxFN = file.path(TaxDir, paste0('train.list.nlevel', nlevel_str, '_marker_index.rda'))

      ##
      if(file.exists(TaxFN)) {
         print(paste("### ", Taxonomy, " Base Training Template (markers & indices) are being read ..."))
         load(TaxFN)
      }else{
         print("### Building Base Training Teamplates ...")
         pre.train.list = NA
         ##
         train.list = build_train_list_default(pre.train.list, 
                                                query.genes=NA, 
                                                TrainDir, 
                                                TaxDir,
                                                prefix="", 
                                                nlevel=nlevel, 
                                                TaxFN=TaxFN)
         ##
         registerDoMC(cores=mc.cores)

         ##
         MI_str = paste0('nlevel', nlevel_str, '_marker_index')
         train.list = build_marker_index_tree_cl (train.list, 
                                                   pre.train.list=pre.train.list, 
                                                   query.genes=NA,
                                                   outdir=file.path(TaxDir, MI_str),
                                                   div_thr=3, 
                                                   subsample_pct=subsample_pct,
                                                   top.n.genes=15, 
                                                   n.group.genes=3000, 
                                                   mc.cores=mc.cores)
         save(train.list, file=train.list$TaxFN)
      }

      print("### Base Training Teamplates Are Ready.\n")

      ##
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
            ##
            registerDoMC(cores=mc.cores)
            train.list = build_marker_index_tree_cl(train.list, 
                                                      pre.train.list=pre.train.list, 
                                                      query.genes,
                                                      outdir=file.path(TaxDir, MI_str),
                                                      div_thr=3, 
                                                      subsample_pct=subsample_pct,
                                                      top.n.genes=15, 
                                                      n.group.genes=3000, 
                                                      mc.cores=mc.cores)
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

#' INFO -- PLEASE ADD -- @CK please document this function
#'
#' @param x  to_be_added
#' @param key  to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_default <- function(pre.train.list=NA, 
                                       query.genes=NA, 
                                       TrainDir, 
                                       TaxDir, 
                                       prefix="", 
                                       nlevel=4, 
                                       TaxFN, 
                                       rm.cl=c())
{
   tmp = load(file.path(TrainDir, "cl.df.rda"))
   #cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)
   if("class_label" %in% colnames(cl.df)){
      cl.df$class_label = gsub("/","_",cl.df$class_label)
   }
   if("neighborhood_label" %in% colnames(cl.df)){
      cl.df$neighborhood_label = gsub("/", "_", cl.df$neighborhood_label)
   }
   if("subclass_label" %in% colnames(cl.df)){
      cl.df$subclass_label = gsub("/","_",cl.df$subclass_label)
   }
   tt21 = table(cl.df$subclass_label, cl.df$class_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$class_label   = tt21.major[cl.df$subclass_label]# cl.df$Level1_label

   tmp=load(file.path(TrainDir, "cl.means.rda"))
   train.cl.dat = cl.means 
   colnames(train.cl.dat) = colnames(cl.means)
   common.cluster = intersect(colnames(train.cl.dat), cl.df$cl)
   if (length(rm.cl)>0) {
      common.cluster = setdiff(common.cluster, rm.cl)
   }

   train.cl.dat = train.cl.dat[, common.cluster]
   cl.df = cl.df %>% filter(cl %in% common.cluster)

   if (!is.list(pre.train.list)) {
      print("=====================================")
      print("  marker/index generation begins...  ")
      print("  it will take a day  :)             ")
      print("=====================================")
   }
   # WholeBrain Parquet
   dsFN = file.path(TrainDir, "de_parquet")
   pairsFN = file.path(TrainDir, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(TrainDir, "cl.bin.rda"))

   select.markers.short = rownames(train.cl.dat)
   if (is.na(query.genes)) query.genes = select.markers.short

   if (is.list(pre.train.list)) {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   } else {
      select.markers = intersect(select.markers.short, query.genes)
   }

   train.cl.dat = train.cl.dat[select.markers,]

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df
   train.list$nlvl      = nlevel
   train.list$dsFN      = dsFN
   train.list$cl.bin    = cl.bin
   train.list$select.markers = select.markers
   train.list$TaxFN     = TaxFN
   train.list$TaxDir    = TaxDir

   return(train.list)
}
