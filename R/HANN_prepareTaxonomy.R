#' run_prepareTaxonomy : call prepareTaxonomy() from h5ad file
#'
#' @param h5adFN     filename of anndata with uns[['HANN']]
#' @export
run_prepareTaxonomy <- function( h5adFN ) {
   library(anndata)
   library(reticulate)
   library(Matrix)  
   adata = read_h5ad(h5adFN)
   bdata = adata$T
   rm(adata)
   ngenes_1 = length(bdata$obs_names)-1
   norm.dat = bdata$chunk_X(0:ngenes_1)
   rownames(norm.dat) = bdata$obs_names
   colnames(norm.dat) = bdata$var_names

   cl           = bdata$uns[['HANN']]$cl
   names(cl)    = bdata$var_names
   if (sum(names(cl) != bdata$var_names) > 0) {
      print("cl and bdata$obs_names do not match")
   }
   cl.df        = bdata$uns[['HANN']]$cl.df
   cl.df        = py_to_r(cl.df)

   cl.hierarchy = bdata$uns[['HANN']]$cl.hierarchy
   AIT.str      = bdata$uns[['HANN']]$taxonomyName
   taxonomy.dir = bdata$uns[['HANN']]$taxonomyDir

   AIT.dir = prepareTaxonomy ( count        = norm.dat,
                               cl           = cl,
                               cl.df        = cl.df,
                               cl.hierarchy = cl.hierarchy,
                               AIT.str      = AIT.str,
                               taxonomy.dir = taxonomy.dir )
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
   count, 
   cl, cl.df.users, cl.hierarchy, AIT.str,
   taxonomy.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/",
   h5ad_out = FALSE ) {
   library(scrattch.bigcat)
   library(dplyr)

   ## Define reference-taxonomy folder and filenames
   AIT.dir = file.path(taxonomy.dir, AIT.str)
   if (!file.exists(AIT.dir)) dir.create(AIT.dir)
   de.dir  = file.path(AIT.dir, "de_parquet")
   sum.dir = file.path(AIT.dir, "de_summary")
   dat.dir = file.path(AIT.dir, "data_parquet")
   h5ad.FN = file.path(AIT.dir, paste0(AIT.str, ".h5ad"))
   pairs.FN = file.path(AIT.dir, "pairs.parquet")
   cl.bin.FN = file.path(AIT.dir, "cl.bin.rda")


   ## convert norm.data matrix to big.dat in parquet, if big.dat is not available (default)
   count = count[, sample_cells(cl,200)]
   cl = cl[colnames(count)]
   count.sparse = as(count, 'sparseMatrix')
   library(data.table)
   big.dat=convert_big.dat_parquet(count.sparse,  dir=dat.dir)

   print("## cluster stat")
   cl.stats = get_cl_stats_big(big.dat, cl=cl, stats=c("means","present","sqr_means"), mc.cores=15)
   save(cl.stats, file=file.path(AIT.dir, "cl.stats.big.rda"))
   cl.means = cl.stats$means
   save(cl.means, file=file.path(AIT.dir, "cl.means.rda"))

   cl.df = rearrange_cl.df (cl.df.users, cl.hierarchy) 
   save(cl.df, file=file.path(AIT.dir, "cl.df.rda"))

   print("## calculate all-pairwise DE genes ")
   if (file.exists(de.dir)) system(paste0("rm -r ", de.dir))
   if (file.exists(sum.dir)) system(paste0("rm -r ", sum.dir))

   de.result = prep_parquet_de_all_pairs(norm.dat=NULL, cl=cl, cl.bin=NULL, mc.cores=30, 
                  pairs.fn=pairs.FN,  cl.bin.fn=cl.bin.FN, cl.means=cl.stats$means, 
                  cl.present=cl.stats$present, cl.sqr.means=cl.stats$sqr_means, 
                  out.dir=de.dir, summary.dir=sum.dir, de.param=de_param())
   return(AIT.dir) 
}
# de_param default values
#     de_param(
#       low.th = 1,
#       padj.th = 0.01,
#       lfc.th = 1,
#       q1.th = 0.5,
#       q2.th = NULL,
#       q.diff.th = 0.7,
#       de.score.th = 150,
#       min.cells = 4,
#       min.genes = 5
#     )

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
   nlevel=length(cl.hierarchy)
   if (nlevel==2) {
      hier_label = c('root', 'cl')
      cl         = cl.df[, cl.hierarchy[2]]
      cluster_id = cl
      df = data.frame(cl, cluster_id)
   }
   if (nlevel==3) {
      hier_label     = c('root', 'subclass_label', 'cl')
      subclass_label = cl.df[, cl.hierarchy[2]]
      cl             = cl.df[, cl.hierarchy[3]]
      cluster_id     = cl
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

