#' build_train_list_default : building train.list
#'
#' @param pre.train.list to_be_added
#' @param query.genes to_be_added
#' @param TaxDir to_be_added
#' @param prefix to_be_added
#' @param nlevel to_be_added
#' @param TaxFN to_be_added
#' @param rm.cl to_be_added
#' @param root.markers to_be_added
#'
#' @return train.list
#'
#' @export

build_train_list_default <- function(pre.train.list=NA, query.genes=NA, TrainDir, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=NULL, root.markers=NULL)
{
   tmp = load(file.path(TrainDir, "cl.df.rda"))
   cl.df = as.data.frame(cl.df)
   cl.df$cl = as.character(cl.df$cl)
   if ("class_label" %in% colnames(cl.df)) {
      cl.df$class_label = gsub("/","_",cl.df$class_label)
   }
   if ("neighborhood_label" %in% colnames(cl.df)) {
      cl.df$neighborhood_label = gsub("/", "_", cl.df$neighborhood_label)
   }
   if ("subclass_label" %in% colnames(cl.df)) {
      cl.df$subclass_label = gsub("/","_",cl.df$subclass_label)
   }
   if ("class_label" %in% colnames(cl.df)) {
      tt21 = table(cl.df$subclass_label, cl.df$class_label)
      tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
      names(tt21.major) = rownames(tt21)
      cl.df$class_label   = tt21.major[cl.df$subclass_label]# cl.df$Level1_label
   }

   tmp=load(file.path(TrainDir, "cl.means.rda"))
   train.cl.dat = cl.means 

   # clusters
   colnames(train.cl.dat) = colnames(cl.means)
   common.cluster = intersect(colnames(train.cl.dat), cl.df$cl)
   if (!is.null(rm.cl)) {
      common.cluster = setdiff(common.cluster, rm.cl)
   }
   cl.df = cl.df %>% filter(cl %in% common.cluster)

   # genes
   select.markers = rownames(train.cl.dat)
   if (is.na(query.genes)) query.genes = select.markers
   if (is.list(pre.train.list)) {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   } else {
      select.markers = intersect(select.markers, query.genes)
   }

   # update train.cl.dat
   train.cl.dat = train.cl.dat[select.markers, common.cluster]

   if (!is.list(pre.train.list)) {
      print("=====================================")
      print("  marker/index generation begins...  ")
      print("  it will take a day  :)             ")
      print("=====================================")
   }
   # WholeBrain Parquet
   dsFN = file.path(TrainDir, "de_parquet")
   load(file.path(TrainDir, "cl.bin.rda"))

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df
   train.list$nlvl      = nlevel
   train.list$dsFN      = dsFN
   train.list$cl.bin    = cl.bin
   train.list$TaxFN     = TaxFN
   train.list$TaxDir    = TaxDir
   train.list$select.markers = select.markers
   train.list$root.markers   = root.markers

   return(train.list)
}
