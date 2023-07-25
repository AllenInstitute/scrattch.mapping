#' Hierarchical approximate nearest neighbor mapping 
#'
#' @param querydat : count matrix (log-normalized)  Ngene X Ncell
#' @param Taxonomy : Taxonomy id
#' @param TaxHome  : diretory path holding taxonomy refrence files
#' @param mapping.method : 'flat(one-step)'/'hierarchy'
#'             run 'flat(one-step)' first
#'            'hierarchy' mapping is recommended for mapping across platforms,
#'             yet it is slow until we come up with further parallelization.
#' @param prebuild : TRUE/FALSE
#'            if prebuild==TRUE : if the template with 'prefix' exists already, use it.
#'            if prebuild==FALSE : do not use preexisting template and build new template
#'                                 with 'prefix'
#' @param newbuild : TRUE/FALSE
#'            When any marker is missing in the query data,
#'            if newbuild==TRUE : it will find markers all over again from the gene list
#'                                in query data. (can take a whole day)
#'            if newbuild==FALSE : it will use the subset of markers which are in the query data 
#'            selected markers are saved using 'prefix' as a key
#' @param prefix : your data descriptor  and platform
#'            '10x_nuclei', 'MFISH', 'patchseq', ....
#'            matching templates are saved using 'prefix' as a key
#'            when you mapping the data from the same platform, you don't need to rebuild the
#'            template for the same platform and use the existing template with prebuild=TRUE
#' @param mc.cores: (default 20) number of cores to be used in parallel processing
#' @param iter : (default=100) number of iteration in mapping with subsampling to estimate
#'               the confidence of mapping
#' @param blocksize : (default=5000) processing by block to avoid the crash
#' @param Taxonomy :   AIT id for available Taxonomies for mapping in the table
#'
#' @return map.freq
#'            cl     : all clusters a sample is mapped in N iterations of mapping with
#'                     sub-sampled markers
#'            freq : frequencies a sample is mapped to each cluster, cl
#'            dist  : distance to cluster template centroid (mean marker gene count)
#'            path.cor : correlation of markers along the path of the hierarchy to the terminal node(cluster), in hierarchical mapping
#' @return best.map.df
#' @param      : all clusters a sample is mapped in N iterations of mapping with
#'                     sub-sampled markers
#'            freq : frequencies a sample is mapped to each cluster, cl
#'            dist  : distance to cluster template centroid (mean marker gene count)
#'            path.cor : correlation of markers along the path of the hierarchy to the terminal node(cluster), in hierarchical mapping
#' @return best.map.df
#'            best.cl  : the cluster a sample is mapped with highest freq in map.freq
#'            prob     : probablity of a sample being mapped to best.cl cluster out of  N iterations
#'            avg.dist : distance to the template cluster mean
#'            avg.path.cor : correlation of markers along the path of the hierarchy to the terminal node(cluster), in hierarchical mapping
#'            avg.cor  : correlation to template cluster mean
#'            cor.zscore : z-normalized value of over avg.cor
#' @return cl.df : taxonomy cluster annotation matched by "cl" (mapped[["best.map.df"]]$best.cl)
#'
#' @export
hannMap <- function(query.dat,
                     Taxonomy="",
                     TaxHome="",
                     prefix="",
                     TaxFN=NA,
                     prebuild=FALSE,
                     newbuild=FALSE,
                     mapping.method=c('flat', 'one-step', 'hierarchy'),
                     nlevel=4,
                     iter=100,
                     mc.cores=7,
                     blocksize=50000,
                     dist.method="cor",
                     topk=1,
                     subsample_pct=0.9,
                     top.n.genes=15,
                     rm.clusters=NA,
                     flag.serial=TRUE,
                     flag.parallel.tmp=FALSE,
                     flag.fuzzy=FALSE){

   ## Check for required files here

   ##
   print("### Training Templates for Query Data")
   train.list = build_train_list_on_taxonomy(TaxFN, 
                                             Taxonomy, 
                                             TaxHome=TaxHome,
                                             query.genes=rownames(query.dat),
                                             prefix=prefix, 
                                             mapping.method=mapping.method,
                                             nlevel=nlevel, 
                                             prebuild=prebuild, 
                                             newbuild=newbuild,
                                             mc.cores=mc.cores, 
                                             subsample_pct=subsample_pct)

   ##
   print("### Mapping on Training Templates")
   query.dat.mapped = mapping_on_taxonomy(query.dat, 
                                          train.list, 
                                          blocksize=blocksize, 
                                          mc.cores=mc.cores,
                                          method=dist.method, 
                                          iter=iter, 
                                          topk=topk,
                                          flag.serial=flag.serial, 
                                          flag.parallel.tmp, 
                                          flag.fuzzy=flag.fuzzy)

   ##
   print("### Adding cluster related info")
   query.dat.mapped$cl.df = train.list$cl.df
   query.dat.mapped$all.markers = train.list$all.markers
   return(query.dat.mapped)
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
