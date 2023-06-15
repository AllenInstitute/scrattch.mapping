#' run_prepareTaxonomy : call prepareTaxonomy() from h5ad file
#'
#' @param h5adFN filename of anndata
#' @export
run_prepareTaxonomy <- function( h5adFN, big.dat.FN=NULL ) {

   adata = read_h5ad(h5adFN)
   norm.dat = t(adata$X)
   rownames(norm.dat) = adata$var_names
   colnames(norm.dat) = adata$obs_names
   cl = adata$uns[['HANN']]$cl
   names(cl) = adata$obs_names

   if (is.null(big.dat.FN)) big.dat=NULL
   else load(big.dat.FN)

   AIT.dir = prepareTaxonomy ( count = norm.dat,
                     cl    = cl,
                     cl.df = adata$uns[['HANN']]$cl.df,
                     cl.hierarchy = adata$uns[['HANN']]$cl.hierarchy,
                     AIT.str      = adata$uns[['taxonomyName']],
                     taxonomy.dir = adata$uns[['taxonomyDir']],
                     big.dat = big.dat )
   print(paste0("Taxonomy is ready in", AIT.dir))
}
#' Prepare taxonomy for optimized tree mapping
#'
#' This function writes an anndata object in the correct format for all downstream mapping algoriths.
#'
#' @param counts count[gene x cell]
#' @param cl assigned cluster
#' @param cl.df.users cluster annotation 
#' @param cl.hierarchy : hierarhcy in clusters : cluster(cl)/subclass_label/neighborhood/root
#' @param AIT.str taxonomy id
#' @param taxonomy.dir Output directly to write h5ad file
#'
#' @import patchseqtools
#' @import scrattch.hicat
#'
#' @return A copy of the anndata object that is also written to disk
#'
#' @export
prepareTaxonomy <- function (
   count, cl, cl.df.users, cl.hierarchy, AIT.str,
   lognormal = NULL,
   taxonomy.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/",
   big.dat=NULL,
   h5ad_out = FALSE ) {

   ## Define taxonomy locations
   AIT.dir = file.path(taxonomy.dir, AIT.str)
   if (!file.exists(AIT.dir)) dir.create(AIT.dir)
   de.dir  = file.path(AIT.dir, "de_parquet")
   sum.dir = file.path(AIT.dir, "de_summary")
   dat.dir = file.path(AIT.dir, "data_parquet")
   h5ad.FN = file.path(AIT.dir, paste0(AIT.str, ".h5ad"))
   pairs.FN = file.path(AIT.dir, "pairs.parquet")
   cl.bin.FN = file.path(AIT.dir, "cl.bin.rda")

   if (is.null(big.dat)) {
      ## convert norm.data matrix to big.dat in parquet
      count.sparse = as(count, 'sparseMatrix')
      big.dat=convert_big.dat_parquet(count.sparse,  dir=dat.dir)
      save(big.dat, file=file.path(AIT.dir, "big.dat.rda"))
   } else {
      save(big.dat, file=file.path(AIT.dir, "big.dat.rda"))
   } 

   ## cluster stat
   cl.stats = get_cl_stats_big(big.dat, cl=cl, stats=c("means","present","sqr_means"), mc.cores=15)
   save(cl.stats, file=file.path(AIT.dir, "cl.stats.big.rda"))
   cl.means = cl.stats$means
   save(cl.means, file=file.path(AIT.dir, "cl.means.rda"))

   cl.df = rearrange_cl.df (cl.df.usrs, cl.hierarchy) 
   save(cl.df, file=file.path(AIT.dir, "cl.df.rda"))

   ## calculate all-pairwise DE genes 
   if (file.exists(de.dir)) system(paste0("rm -r ", de.dir))
   if (file.exists(sum.dir)) system(paste0("rm -r ", sum.dir))

   de.result = prep_parquet_de_all_pairs(norm.dat=NULL, cl=cl, cl.bin=NULL, mc.cores=30, 
                  pairs.fn=pairs.FN,  cl.bin.fn=cl.bin.FN, cl.means=cl.stats$means, 
                  cl.present=cl.stats$present, cl.sqr.means=cl.stats$sqr_means, 
                  out.dir=de.dir, summary.dir=sum.dir)
   return(AIT.dir) 
}


#' Rearrange cl.df to internal format labels
#'
#' This function writes an anndata object in the correct format for all downstream mapping algoriths.
#'
#' @param counts count[gene x cell]
#' @param cl assigned cluster
#' @return A copy of the anndata object that is also written to disk
#'
#' @export
rearrange_cl.df <- function (cl.df, cl.hierarchy) {

   nlevel=1 + length(cl.hierarhcy)
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

