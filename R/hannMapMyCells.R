hannMapMyCells = function(AIT_anndata,
                          anndata_path,
                          extended_result_path,
                          chunk_size = 1000,
                          n_processors = 3,
                          normalization = "log2CPM",
                          tmp_dir = NULL){
  delete_tmp_folder <- FALSE
  tryCatch(
    {
      cell_type_mapper <- import("cell_type_mapper")
  
      if(missing(anndata_path)){
        anndata_path = file.path(AIT_anndata$uns$taxonomyDir, paste0(AIT_anndata$uns$taxonomyName, ".h5ad"))
      }

      if (missing(tmp_dir) || is.null(tmp_dir) || tmp_dir == "") {
        tmp_dir <- paste0("./temp_folder_", format(Sys.time(), "%Y%m%d-%H%M%S"))
        dir.create(tmp_dir)
        delete_tmp_folder <- TRUE
      }

      extended_result_path <- paste0(".", file.path(tmp_dir, "mapmycells_results.json"))
      precomputed_stats_path <- toString(cell_type_mapper$utils$output_utils$uns_to_precomputed_stats(h5ad_path=anndata_path, 
                                         uns_key="MapMyCells_HANN_precomp_stats", tmp_dir=tmp_dir))
      print(precomputed_stats_path)

      serialized_query_markers <- cell_type_mapper$utils$anndata_utils$read_uns_from_h5ad(h5ad_path=anndata_path)$MapMyCells_query_markers
      query_markers_filename <- paste0(paste0("query_markers_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".h5")
      query_markers_output_path <- file.path(tmp_dir, query_markers_filename)

      query_markers_json = toJSON(clean_for_json(serialized_query_markers))
      write(query_markers_json, query_markers_output_path)

      config <- list(
          'precomputed_stats'= list(
              'path' = precomputed_stats_path
          ),
          'query_markers'= list(
              'serialized_lookup' = query_markers_output_path
          ),
          'query_path' = anndata_path,
          'extended_result_path' = extended_result_path,
          'type_assignment' = list(
              'normalization' = normalization,
              'chunk_size' = chunk_size,
              'n_processors' = n_processors
          )
      )
  
      runner <- cell_type_mapper$cli$from_specified_markers$FromSpecifiedMarkersRunner(
          args=c(),
          input_data=config)
      runner$run()

      AIT_anndata = read_h5ad(anndata_path)
      return(AIT_anndata)
    },
    error = function(e) {
      errorMessage <- conditionMessage(e)
      cat("Error message:", errorMessage, "\n")
    },
    finally = {
      if (delete_tmp_folder) {
        unlink(tmp_dir, recursive = TRUE)
      }
    }
  )
}