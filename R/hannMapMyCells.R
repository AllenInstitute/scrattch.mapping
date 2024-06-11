hannMapMyCells = function(AIT_anndata,
                          anndata_path,
                          extended_result_path,
                          chunk_size = 1000,
                          n_processors = 3,
                          normalization = "log2CPM",
                          tmp_dir = NULL){
  if(missing(anndata_path)){
    anndata_path = file.path(AIT_anndata$uns$taxonomyDir, paste0(AIT_anndata$uns$taxonomyName, ".h5ad"))
  }
  source_python("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Inkar/repositories/scrattch.mapping/R/mapMyCellsMapping.py")
  fromSpecifiedMarkersRunner(anndata_path, extended_result_path, chunk_size, n_processors, normalization, tmp_dir)

  AIT_anndata = read_h5ad(anndata_path)
  return(AIT_anndata)
}