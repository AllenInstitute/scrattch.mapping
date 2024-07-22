#' @export

mapHANNMapMyCells = function(AIT_anndata,
                          anndata_path=NULL,
                          chunk_size = 1000,
                          n_processors = 3,
                          normalization = "log2CPM",
                          tmp_dir = NULL,
                          user_extended_result_path=NULL,
                          user_precomp_stats_path=NULL,
                          user_query_markers_path=NULL){
  tryCatch(
    {
      temp_folder = tmp_dir
      extended_result_path = user_extended_result_path

      if (is.null(AIT_anndata$uns$hierarchical[[AIT.anndata$uns$mode]])){
        stop(paste("ERROR. Provided mode doesn't exists for hierarchical:", AIT.anndata$uns$mode))
      }

      if(is.null(anndata_path)){
        anndata_path = file.path(AIT_anndata$uns$taxonomyDir, paste0(AIT_anndata$uns$taxonomyName, ".h5ad"))
      }

      if (is.null(temp_folder) || temp_folder == "") {
        temp_folder <- paste0("temp_folder_", format(Sys.time(), "%Y%m%d-%H%M%S"))
        temp_folder <- file.path(getwd(), temp_folder)
        dir.create(temp_folder)
      }

      if (is.null(extended_result_path)) {
        extended_result_path <- file.path(temp_folder, "mapmycells_results.json")
      }
      
      if (file.exists(extended_result_path)) {
        stop(paste("ERROR. Extended result path already exists:", extended_result_path))
      }
      precomp_stats_output_path = user_precomp_stats_path
      precomp_stats_output_path = get_precomp_stats(AIT_anndata, anndata_path, 
                                                    precomp_stats_output_path, temp_folder)

      query_markers_output_path = user_query_markers_path
      query_markers_output_path = get_marker_genes(AIT_anndata, query_markers_output_path, temp_folder)

      get_mapmycells_results(anndata_path, extended_result_path, 
                             precomp_stats_output_path, query_markers_output_path, 
                             normalization, chunk_size, n_processors)

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
      # remove the temp folder is it was code generated
      if (is.null(tmp_dir) || tmp_dir == "") {
        print(paste("Deleting temp folder", temp_folder))
        unlink(temp_folder, recursive = TRUE)
      }
      else {
        # remove the files is they were code generated
        if(is.null(user_precomp_stats_path)) {
          file.remove(precomp_stats_output_path)
        }
        if(is.null(user_query_markers_path)) {
          file.remove(query_markers_output_path)
        }
        if(is.null(user_extended_result_path)) {
          file.remove(extended_result_path)
        }
      }
    }
  )
}

get_precomp_stats = function(AIT_anndata, anndata_path, precomp_stats_output_path, temp_folder) {
  # retrieve precompute stats from anndata
  if(is.null(precomp_stats_output_path)) {
    uns_keys = list("hierarchical", AIT_anndata$uns$mode, "precomp_stats")
    precomp_stats_output_path <- toString(cell_type_mapper$utils$output_utils$uns_to_precomputed_stats(h5ad_path=anndata_path, 
                                      uns_keys_list=uns_keys, tmp_dir=temp_folder))
  }
  return(precomp_stats_output_path)
}

get_marker_genes = function(AIT_anndata, query_markers_output_path, temp_folder) {
  # retrieve query markers from anndata
  if(is.null(query_markers_output_path)) {
    serialized_query_markers <- AIT_anndata$uns$hierarchical[[AIT.anndata$uns$mode]][["query_markers"]]
    query_markers_filename <- paste0(paste0("query_markers_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".json")
    query_markers_output_path <- file.path(temp_folder, query_markers_filename)

    query_markers_json = toJSON(
                            cell_type_mapper$utils$utils$clean_for_uns_deserialization(
                            cell_type_mapper$utils$utils$clean_for_json(serialized_query_markers)),
                          pretty = TRUE)
    write(query_markers_json, query_markers_output_path)
  }
  return(query_markers_output_path)
}

get_mapmycells_results = function(anndata_path, extended_result_path, 
                precomp_stats_output_path, query_markers_output_path, 
                normalization, chunk_size, n_processors) {
  config <- list(
    'query_path' = anndata_path,
    'extended_result_path' = extended_result_path,
    'precomputed_stats'= list(
        'path' = precomp_stats_output_path
    ),
    'query_markers'= list(
        'serialized_lookup' = query_markers_output_path
    ),
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
}