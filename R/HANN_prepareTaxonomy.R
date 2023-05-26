#' Prepare taxonomy for optimized tree mapping
#'
#' This function writes an anndata object in the correct format for all downstream mapping algoriths.
#'
#' @param counts count[gene x cell]
#' @param cl assigned cluster
#' @param cl.df cluster anno with hierarchy : cluster(cl)/subclass_label/neighborhood/root
#' @param AIT.str taxonomy id
#' @param lognormal a logCPM cell x gene matrix.  If counts provided, logCPM is calculated.
#' @param taxonomy.dir Output directly to write h5ad file
#'
#' @import patchseqtools
#' @import scrattch.hicat
#'
#' @return A copy of the anndata object that is also written to disk
#'
#' @export
prepareTaxonomy <- function(count, cl, cl.df, AIT.str,
   lognormal = NULL,
   taxonomy.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/"){

   ## Define taxonomy locations
   AIT.dir = file.path(taxonomy.dir, AIT.str)
   if (!file.exists(AIT.dir)) dir.create(AIT.dir)
   de.dir  = file.path(AIT.dir, "de_parquet")
   sum.dir = file.path(AIT.dir, "de_summary")
   dat.dir = file.path(AIT.dir, "data_parquet")
   h5ad.FN = file.path(AIT.dir, paste0(AIT.str, ".h5ad"))
   pairs.FN = file.path(AIT.dir, "pairs.parquet")
   cl.bin.FN = file.path(AIT.dir, "cl.bin.rda")

   ## count matrix to sparse matrix
   count.sparse = as(count, 'sparseMatrix')
   if (!lognormal) {
      norm.dat.sparse = logCPM(count.sparse)
   } else {
      norm.dat.sparse = count.sparse
      raw.sparse = 2^(count.sparse) - 1
      rownames(raw.sparse) = rownames(count.sparse)
      colnames(raw.sparse) = colnames(count.sparse)
      count.sparse = raw.sparse
   }

   ## convert norm.data matrix to big.dat in parquet
   big.dat=convert_big.dat_parquet(norm.dat.sparse,  dir=dat.dir)
   save(big.dat, file="big.dat.rda")
   
   ## cluster stat
   cl.stats = get_cl_stats_big(big.dat, cl=cl, stats=c("means","present","sqr_means"), mc.cores=15)
   save(cl.stats, file=file.path(AIT.dir, "cl.stats.big.rda"))
   cl.means = cl.stats$means
   save(cl.means, file=file.path(AIT.dir, "cl.means.rda"))
   save(cl.df, file=file.path(AIT.dir, "cl.df.rda"))

   ## calculate all-pairwise DE genes 
   if (file.exists(de.dir)) system(paste0("rm -r ", de.dir))
   if (file.exists(sum.dir)) system(paste0("rm -r ", sum.dir))

   de.result = prep_parquet_de_all_pairs(norm.dat=NULL, cl=cl, mc.cores=30, 
                  pairs.fn=pairs.FN,  cl.bin.fn=cl.bin.FN,
                  cl.means=cl.stats$means, cl.present=cl.stats$present, cl.sqr.means=cl.stats$sqr_means, 
                  out.dir=de.dir, summary.dir=sum.dir)
   
   ##
   print("===============================")
   print(AIT.dir)
   print(dir(AIT.dir))
   print("===============================")

   ## prepare h5ad
   tmp.cl.df = cl.df[match(cl, cl.df$cl),]
   row.names(tmp.cl.df) = names(cl)
   adata = AnnData(X=t(norm.dat.sparse), obs=tmp.cl.df, var=colnames(norm.dat.sparse))#, raw=t(count.sparse))
   adata$obs = tmp.cl.df
   adata$raw = t(count.sparse)
   write_h5ad(adata, h5ad.FN)

   ##
   return(adata)
}
