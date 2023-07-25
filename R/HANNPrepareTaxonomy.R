#' Prepare taxonomy for optimized tree mapping
#'
#' This function writes an anndata object in the correct format for all downstream mapping algoriths.
#'
#' @param AIT.anndata scrattch.mapping .h5ad object from `loadTaxonomy()`
#' @param big.dat pre-computed big.data object, @CK can you provide more info on what this is?
#'
#' @import patchseqtools
#' @import scrattch.hicat
#'
#' @return A copy of the anndata object that is also written to disk
#'
#' @export
prepareHANNTaxonomy = function(AIT.anndata, 
                               big.dat=NULL){

   ##                              
   count = t(AIT.anndata$X)
   rownames(count) = AIT.anndata$var_names
   colnames(count) = AIT.anndata$obs_names
   cl = AIT.anndata$uns[['HANN']]$cl
   names(cl) = AIT.anndata$obs_names

   ## Define taxonomy locations
   HANN.dir = file.path(AIT.anndata$uns$taxonomyDir, paste0(AIT.anndata$uns$taxonomyName, "_HANN"))
   de.dir  = file.path(HANN.dir, "de_parquet")
   sum.dir = file.path(HANN.dir, "de_summary")
   dat.dir = file.path(HANN.dir, "data_parquet")
   pairs.FN = file.path(HANN.dir, "pairs.parquet")
   cl.bin.FN = file.path(HANN.dir, "cl.bin.rda")

   ##
   if (is.null(big.dat)) {
      ## convert norm.data matrix to big.dat in parquet
      count = count[, sample_cells[cl,200]]
      cl = cl[colnames(count)]
      count.sparse = as(count, 'sparseMatrix')
      big.dat=convert_big.dat_parquet(count.sparse,  dir=dat.dir)
      save(big.dat, file=file.path(HANN.dir, "big.dat.rda"))
   } else {
      save(big.dat, file=file.path(HANN.dir, "big.dat.rda"))
   } 

   ## cluster stat
   cl.stats = get_cl_stats_big(big.dat, cl=cl, stats=c("means","present","sqr_means"), mc.cores=15)
   save(cl.stats, file=file.path(HANN.dir, "cl.stats.big.rda"))
   cl.means = cl.stats$means
   save(cl.means, file=file.path(HANN.dir, "cl.means.rda"))
   cl.df = rearrange_cl.df(cl.df.users, cl.hierarchy) 
   save(cl.df, file=file.path(HANN.dir, "cl.df.rda"))

   ## calculate all-pairwise DE genes 
   ## I would prefer to remove these lines @CK.
   if(file.exists(de.dir)) system(paste0("rm -r ", de.dir))
   if(file.exists(sum.dir)) system(paste0("rm -r ", sum.dir))

   ##
   de.result = prep_parquet_de_all_pairs(norm.dat=NULL, 
                                          cl=cl, 
                                          cl.bin=NULL, 
                                          mc.cores=30, 
                                          pairs.fn=pairs.FN,  
                                          cl.bin.fn=cl.bin.FN, 
                                          cl.means=cl.stats$means, 
                                          cl.present=cl.stats$present, 
                                          cl.sqr.means=cl.stats$sqr_means, 
                                          out.dir=de.dir, 
                                          summary.dir=sum.dir)
   return(HANN.dir)
}

#' Rearrange cl.df to internal format labels
#'
#' This function writes an anndata object in the correct format for all downstream mapping algoriths.
#'
#' @param counts count[gene x cell]
#' @param cl assigned cluster
#' @return A copy of the anndata object that is also written to disk
#'
#' @keywords internal
rearrange_cl.df = function (cl.df, cl.hierarchy) {
   ##
   nlevel=length(cl.hierarchy)
   if (nlevel==2) {
      hier_label = c('root', 'cl')
      cl         = cl.df[, cl.hierarchy[2]]
      cluster_id = cl
      df = data.frame(cl, cluster_id)
   }
   if (nlevel==3) {
      hier_label  = c('root', 'class_label', 'cl')
      class_label = cl.df[, cl.hierarchy[2]]
      cl          = cl.df[, cl.hierarchy[3]]
      cluster_id  = cl
      df = data.frame(cl, cluster_id, subclass_label)
   }
   if (nlevel==4) {
      hier_label = c('root', 'class_label', 'subclass_label', 'cl')
      class_label    = cl.df[, cl.hierarchy[2]]
      subclass_label = cl.df[, cl.hierarchy[3]]
      cl             = cl.df[, cl.hierarchy[4]]
      cluster_id = cl
      df = data.frame(cl, cluster_id, subclass_label, class_label)
   }
   library(dplyr)
   rownames(df) = df$cl
   return(df)
}

