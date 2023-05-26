#' INFO -- PLEASE ADD --
#'
#' @param x  to_be_added
#' @param key  to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB17_FB <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   tmp = load(file.path(TaxDir, "cl.clean.rda"))
   #cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "cl.means.clean.rda"))
   train.cl.dat = cl.means 
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(TaxDir, "select.markers.short.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x  to_be_added
#' @param key  to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB17_TR <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   tmp = load(file.path(TaxDir, "cl.clean.rda"))
   cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "cl.means.clean.rda"))
   train.cl.dat = cl.means.clean 
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(TaxDir, "select.markers.short.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_20 <- function(pre.train.list=NA, query.genes=NA, TrainDir, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   tmp = load(file.path(TrainDir, "cl.df.rda"))
   #cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)
   cl.df$subclass_label = gsub("/","_",cl.df$subclass_label)
   cl.df$neighborhood_label = gsub("/", "_", cl.df$neighborhood_label)
   cl.df$Level2_label = cl.df$subclass_label
   cl.df$Level1_label = cl.df$neighborhood_label
   tt21 = table(cl.df$subclass_label, cl.df$neighborhood_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$neighborhood   = tt21.major[cl.df$subclass_label]# cl.df$Level1_label

   tmp=load(file.path(TrainDir, "cl.means.rda"))
   train.cl.dat = cl.means 
   colnames(train.cl.dat) = gsub("\\.0$", "", colnames(cl.means))
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB17_BG <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   tmp = load(file.path(TaxDir, "BG.cl.df.rda"))
   cl.df = as.data.frame(BG.cl.df)
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "BG.cl.means.rda"))
   train.cl.dat = BG.cl.means 
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))

   load(file.path(TaxDir, "select.markers.short.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB17_FB <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   tmp = load(file.path(TaxDir, "cl.clean.rda"))
   #cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "cl.means.clean.rda"))
   train.cl.dat = cl.means 
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(TaxDir, "select.markers.short.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB17_TR <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   tmp = load(file.path(TaxDir, "cl.clean.rda"))
   cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "cl.means.clean.rda"))
   train.cl.dat = cl.means.clean 
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(TaxDir, "select.markers.short.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_TH17_2 <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   tmp = load(file.path(TaxDir, "cl.df.rda"))
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "cl.means.rda"))
   train.cl.dat = cl.means
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   #ds = open_dataset(dsFN)
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))

   load(file.path(TaxDir, "select.markers.TH.rda")) # pooled
   if (is.na(query.genes)) query.genes = select.markers.TH

   if (is.list(pre.train.list)) {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   } else {
      select.markers = intersect(select.markers.TH, query.genes)
   }

   train.cl.dat = train.cl.dat[select.markers,]

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df
   train.list$nlvl      = nlevel
   train.list$dsFN      = dsFN
   #train.list$ds        = ds
   train.list$cl.bin    = cl.bin
   train.list$select.markers = select.markers
   train.list$all.markers = select.markers
   train.list$TaxFN     = TaxFN
   train.list$TaxDir    = TaxDir

   return(train.list)
}

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB17 <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   tmp = load(file.path(TaxDir, "cl.clean.rda"))
   cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "cl.means.clean.rda"))
   train.cl.dat = cl.means.clean 
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(TaxDir, "select.markers.short.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WBnuc <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   #tmp = load(file.path(WBDIR, "anno.df.rda"))
   tmp = load(file.path(TaxDir, "cl.clean.rda"))
   cl.df = cl.df.clean

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "train.cl.means.rda"))
   train.cl.dat = cl.means 
   common.cluster = intersect(colnames(train.cl.dat),cl.df$cl)
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(WBDIR, "conserved.markers.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB16nuc <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   #tmp = load(file.path(WBDIR, "anno.df.rda"))
   tmp = load(file.path(WBDIR, "cl.clean.rda"))
   cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(WBDIR, "nuclei.impute.cl.means.rda"))
   train.cl.dat = cl.means 
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(WBDIR, "select.markers.short.rda")) # pooled
   select.markers.short  = intersect(rownames(cl.means), select.markers.short)

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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB16 <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   #tmp = load(file.path(WBDIR, "anno.df.rda"))
   tmp = load(file.path(WBDIR, "cl.clean.rda"))
   cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(WBDIR, "cl.means.clean.rda"))
   train.cl.dat = cl.means.clean 
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(WBDIR, "select.markers.short.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB13 <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   #tmp = load(file.path(WBDIR, "anno.df.rda"))
   tmp = load(file.path(WBDIR, "cl.clean.rda"))
   cl.df = cl.df.clean
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "train.cl.means.rda"))
   train.cl.dat = cl.means 
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
   CLSDIR = gsub("/v3", "",WBDIR)
   dsFN = file.path(CLSDIR, "comb_de_parquet")
   pairsFN = file.path(CLSDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(CLSDIR, "cl.bin.rda"))
   load(file.path(WBDIR, "select.markers.short.rda")) # pooled
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WBint <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   CLSDIR = WBDIR
   tmp=load(file.path(WBDIR, "anno.df.rda"))
   cl.df = read.csv(file.path(WBDIR, "short.cl.df.csv"))
   cl.df$cl = as.character(cl.df$cluster_id)
   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = cl.df$Level1_label
   tmp=load(file.path(TaxDir, "train.mmean.list.10Xv3.rda"))
   train.cl.dat = train.mmean.dat[["cluster"]]
   common.cluster = intersect(colnames(train.cl.dat),as.character(anno.df$cl))
   if (length(rm.cl)>0) {
      common.cluster = setdiff(common.cluster, rm.cl)
   }

   train.cl.dat = train.cl.dat[, common.cluster]
   rownames(train.cl.dat) = train.mmean.dat[["genename"]]

   cl.df = cl.df %>% filter(cl %in% common.cluster)

   if (!is.list(pre.train.list)) {
      print("=====================================")
      print("  marker/index generation begins...  ")
      print("  it will take a day  :)             ")
      print("=====================================")
   }
   # WholeBrain Parquet
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))
   load(file.path(WBDIR, "select.markers.rda")) # pooled
   if (is.na(query.genes)) query.genes = select.markers

   if (is.list(pre.train.list)) {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   } else {
      select.markers = intersect(select.markers, query.genes)
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

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_WB <- function( pre.train.list=NA, query.genes, CLSDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   tmp=load(file.path(CLSDIR, "cl.HQ.subdiv.rda"))
   cl.df.HQ = cl.df.HQ.subdiv

   # WB train.mmeta
   tmp=load(file.path(CLSDIR, "train.mmean.list.HQ.rda"))
   train.cl.dat = train.mmean.dat[["cluster"]]
   common.cluster = intersect(colnames(train.cl.dat),as.character(cl.df.HQ$cl))

   if (length(rm.cl)>0) {
      common.cluster = setdiff(common.cluster, rm.cl)
   }
   train.cl.dat = train.cl.dat[, common.cluster]
   rownames(train.cl.dat) = train.mmean.dat[["genename"]]

   cl.df.HQ = cl.df.HQ %>% filter(cl %in% common.cluster)

   if (!is.list(pre.train.list)) {
      print("=====================================")
      print("  marker/index generation begins...  ")
      print("  it will take a day  :)             ")
      print("=====================================")
   }
   # WholeBrain Parque
   dsFN = file.path(CLSDIR, "de_parquet")
   pairsFN = file.path(CLSDIR, "all.pairs.parqeut")
   load(file.path(CLSDIR, "select.markers.rda")) # pooled
   all.pairs = read_parquet(pairsFN)

   if (!is.list(pre.train.list)) {
      select.markers = intersect(select.markers, query.genes)
   } else {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   }

   train.cl.dat = train.cl.dat[select.markers,]

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df.HQ
   train.list$nlvl      = nlevel
   train.list$dsFN      = dsFN
   train.list$all.pairs    = all.pairs
   train.list$select.markers = select.markers
   train.list$TaxFN  = TaxFN
   train.list$TaxDir = TaxDir

   return(train.list)
}

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_TH13_2 <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   tmp = load(file.path(TaxDir, "cl.df.rda"))
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "cl.means.rda"))
   train.cl.dat = cl.means
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   #ds = open_dataset(dsFN)
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))

   load(file.path(TaxDir, "select.markers_TH.rda")) # pooled
   if (is.na(query.genes)) query.genes = select.markers

   if (is.list(pre.train.list)) {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   } else {
      select.markers = intersect(select.markers, query.genes)
   }

   train.cl.dat = train.cl.dat[select.markers,]

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df
   train.list$nlvl      = nlevel
   train.list$dsFN      = dsFN
   #train.list$ds        = ds
   train.list$cl.bin    = cl.bin
   train.list$select.markers = select.markers
   train.list$all.markers = select.markers
   train.list$TaxFN     = TaxFN
   train.list$TaxDir    = TaxDir

   return(train.list)
}

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_TH16_2 <- function(pre.train.list=NA, query.genes=NA, WBDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   tmp = load(file.path(TaxDir, "cl.df.rda"))
   cl.df$cl = as.character(cl.df$cl)

   tt21 = table(cl.df$Level2_label, cl.df$Level1_label)
   tt21.major = colnames(tt21)[apply(tt21, 1, function(x){which.max(x)})]
   names(tt21.major) = rownames(tt21)

   cl.df$subclass_label = cl.df$Level2_label
   cl.df$neighborhood   = tt21.major[cl.df$Level2_label]# cl.df$Level1_label

   tmp=load(file.path(TaxDir, "cl.means.rda"))
   train.cl.dat = cl.means
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
   dsFN = file.path(WBDIR, "comb_de_parquet")
   pairsFN = file.path(WBDIR, "pairs.parquet")
   #ds = open_dataset(dsFN)
   all.pairs = read_parquet(pairsFN)
   load(file.path(WBDIR, "cl.bin.rda"))

   load(file.path(TaxDir, "select.markers_TH.rda")) # pooled
   if (is.na(query.genes)) query.genes = select.markers

   if (is.list(pre.train.list)) {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   } else {
      select.markers = intersect(select.markers, query.genes)
   }

   train.cl.dat = train.cl.dat[select.markers,]

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df
   train.list$nlvl      = nlevel
   train.list$dsFN      = dsFN
   #train.list$ds        = ds
   train.list$cl.bin    = cl.bin
   train.list$select.markers = select.markers
   train.list$all.markers = select.markers
   train.list$TaxFN     = TaxFN
   train.list$TaxDir    = TaxDir

   return(train.list)
}

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_FBLO10 <- function( pre.train.list=NA, query.genes, CLSDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   tmp=load(file.path(CLSDIR, "cl.clean.rda"))
   cl.df.clean = cl.df.neighborhood

   # FB train.mmeta
   tmp=load(file.path(CLSDIR, "train.list.cl.mean.90.rda"))
   train.cl.dat = FB.90.cl.mean
   common.cluster = intersect(colnames(train.cl.dat),as.character(cl.df.clean$cl))
   train.cl.dat = train.cl.dat[, common.cluster]

   cl.df.clean = cl.df.clean %>% filter(cl %in% common.cluster)

   if (!is.list(pre.train.list)) {
      print("=====================================")
      print("  marker/index generation begins...  ")
      print("  it will take a day  :)             ")
      print("=====================================")
   }
   # WholeBrain Parquet
   dsFN = file.path(CLSDIR, "de_parquet")
   pairsFN = file.path(CLSDIR, "pair.parquet")
   load(file.path(CLSDIR, "select.markers.rda")) # pooled
   all.pairs = read_parquet(pairsFN)

   if (!is.list(pre.train.list)) {
      select.markers = intersect(select.markers, query.genes)
   } else {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   }

   train.cl.dat = train.cl.dat[select.markers,]

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df.clean
   train.list$nlvl      = nlevel 
   train.list$dsFN      = dsFN
   train.list$all.pairs   = all.pairs
   train.list$select.markers = select.markers
   train.list$TaxFN     = TaxFN
   train.list$TaxDir = TaxDir

   return(train.list)
}

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_FB <- function( pre.train.list=NA, query.genes, CLSDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   tmp=load(file.path(CLSDIR, "cl.clean.rda"))
   cl.df.clean = cl.df.neighborhood

   # FB train.mmeta
   tmp=load(file.path(CLSDIR, "train.mmean.list.10X_cells_v3.clean.rda"))
   train.cl.dat = train.mmean.dat[["cluster"]]
   common.cluster = intersect(colnames(train.cl.dat),as.character(cl.df.clean$cl))
   train.cl.dat = train.cl.dat[, common.cluster]
   rownames(train.cl.dat) = train.mmean.dat[["genename"]]

   cl.df.clean = cl.df.clean %>% filter(cl %in% common.cluster)

   if (!is.list(pre.train.list)) {
      print("=====================================")
      print("  marker/index generation begins...  ")
      print("  it will take a day  :)             ")
      print("=====================================")
   }
   # WholeBrain Parquet
   dsFN = file.path(CLSDIR, "de_parquet")
   pairsFN = file.path(CLSDIR, "pair.parquet")
   load(file.path(CLSDIR, "select.markers.rda")) # pooled
   all.pairs = read_parquet(pairsFN)

   if (!is.list(pre.train.list)) {
      select.markers = intersect(select.markers, query.genes)
   } else {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   }

   train.cl.dat = train.cl.dat[select.markers,]

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df.clean
   train.list$nlvl      = nlevel 
   train.list$dsFN      = dsFN
   train.list$all.pairs   = all.pairs
   train.list$select.markers = select.markers
   train.list$TaxFN     = TaxFN
   train.list$TaxDir = TaxDir

   return(train.list)
}

#' INFO -- PLEASE ADD --
#'
#' @param x to_be_added
#' @param key to_be_added
#'
#' @return ???
#'
#' @keywords internal
build_train_list_BG <- function( pre.train.list=NA, query.genes, BGDIR, TaxDir, prefix="", nlevel=4, TaxFN, rm.cl=c())
{
   BGDIR = "/home/changkyul/CK/Mapping_PatchSeq/BasalGanglia/FB_Joint_Taxonomy"
   load(file.path(BGDIR, "BG.cl.df.clean.rda"))
   cl.df.clean = BG.cl.df.clean

   #BG train.mmeta
   tmp=load(file.path(BGDIR, "train.mmean.list.FB.10X_cells_v3.clean.BG.rda"))
   train.cl.dat = train.mmean.dat[["cluster"]]
   common.cluster = intersect(colnames(train.cl.dat),as.character(cl.df.clean$cl))
   train.cl.dat = train.cl.dat[, common.cluster]

   cl.df.clean = cl.df.clean %>% filter(cl %in% common.cluster)

   # ForeBrain Parquet
   FBDIR = "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/forebrain_new"

   dsFN = file.path(FBDIR, "de_parquet")
   pairsFN = file.path(FBDIR, "pair.parquet")
   load(file.path(FBDIR, "select.markers.rda")) # pooled
   all.pairs = read_parquet(pairsFN)

   if (!is.list(pre.train.list)) {
      print("=====================================")
      print("  marker/index generation begins...  ")
      print("  it will take a day  :)             ")
      print("=====================================")
   }

   if (!is.list(pre.train.list)) {
      all.markers = NA
      select.markers = intersect(select.markers, query.genes)
   } else {
      select.markers = intersect(pre.train.list$all.markers, query.genes)
   }

   select.markers = intersect(select.markers, query.genes)
   train.cl.dat = train.cl.dat[select.markers,]

   train.list <-list()
   train.list$cl.dat    = train.cl.dat
   train.list$cl.df     = cl.df.clean
   train.list$nlvl      = 4 
   train.list$dsFN      = dsFN
   train.list$all.pairs   = all.pairs
   train.list$all.markers = all.markers
   train.list$select.markers = select.markers
   train.list$TaxFN  = TaxFN
   train.list$TaxDir = TaxDir

   return(train.list)
}

