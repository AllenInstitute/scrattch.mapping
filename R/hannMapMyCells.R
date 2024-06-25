#' @export

mapHANNMapMyCells = function(AIT_anndata,
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
        tmp_dir <- paste0("temp_folder_", format(Sys.time(), "%Y%m%d-%H%M%S"))
        tmp_dir <- file.path(getwd(), tmp_dir)
        dir.create(tmp_dir)
        delete_tmp_folder <- TRUE
      }

      if (missing(extended_result_path)) {
        extended_result_path <- file.path(tmp_dir, "mapmycells_results.json")
      }

      precomp_stats_output_path <- toString(cell_type_mapper$utils$output_utils$uns_to_precomputed_stats(h5ad_path=anndata_path, 
                                         uns_key="MapMyCells_HANN_precomp_stats", tmp_dir=tmp_dir))
      print(precomp_stats_output_path)

      serialized_query_markers <- cell_type_mapper$utils$anndata_utils$read_uns_from_h5ad(h5ad_path=anndata_path)$MapMyCells_HANN_query_markers
      query_markers_filename <- paste0(paste0("query_markers_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".json")
      query_markers_output_path <- file.path(tmp_dir, query_markers_filename)

      query_markers_json = toJSON(
                              cell_type_mapper$utils$utils$clean_for_uns_deserialization(
                              cell_type_mapper$utils$utils$clean_for_json(serialized_query_markers)),
                            pretty = TRUE)
      write(query_markers_json, query_markers_output_path)

      config <- list(
          'precomputed_stats'= list(
              'path' = precomp_stats_output_path
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

      ## Extract mapping results
      mapmycells_results_json = fromJSON(extended_result_path)

      ## Build mapping results dataframe
      mappingResults=list()
      for (hierarcy_level in names(mapmycells_results_json$results)) {
        if (hierarcy_level != "cell_id") {
          mappingResults[[paste0(hierarcy_level, ".assignment")]] = 
              mapmycells_results_json$results[[hierarcy_level]]$assignment
          mappingResults[[paste0(hierarcy_level, ".bootstrapping_probability")]] = 
              mapmycells_results_json$results[[hierarcy_level]]$bootstrapping_probability
          mappingResults[[paste0(hierarcy_level, ".avg_correlation")]] = 
              mapmycells_results_json$results[[hierarcy_level]]$avg_correlation
        }
      }
      
      ## Combine results into dataframe
      mappingAnno = Reduce(cbind, mappingResults)
      rownames(mappingAnno) = mapmycells_results_json$results$cell_id
      colnames(mappingAnno) = names(mappingResults)

      ## Build mapping class object
      resultAnno <- mappingClass(annotations = as.data.frame(mappingAnno),
                                 detailed_results = list("MapMyCellsHANN" = mapmycells_results_json))
      
      ## Return annotations and detailed model results
      return(resultAnno)
    },
    error = function(e) {
      errorMessage <- conditionMessage(e)
      cat("Error message:", errorMessage, "\n")
    },
    finally = {
      if (delete_tmp_folder) {
        print(paste("Deleting temp folder", tmp_dir))
        unlink(tmp_dir, recursive = TRUE)
      }
    }
  )
}