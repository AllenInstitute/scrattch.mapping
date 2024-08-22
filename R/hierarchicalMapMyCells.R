#' @export

hierarchicalMapMyCells = function(AIT_anndata,
                            query_data,
                            chunk_size = 1000,
                            n_processors = 3,
                            normalization = "log2CPM",
                            tmp_dir = NULL,
                            user_extended_result_path=NULL,
                            user_precomp_stats_path=NULL,
                            user_query_markers_path=NULL){
  print("Hierarchical mapping")
  tryCatch(
    {
      temp_folder = tmp_dir
      extended_result_path = user_extended_result_path

      if (is.null(AIT_anndata$uns$hierarchical[[AIT.anndata$uns$mode]])){
        stop(paste("ERROR. Provided mode doesn't exists for hierarchical:", AIT.anndata$uns$mode))
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

      taxonomy_anndata_path = file.path(AIT_anndata$uns$taxonomyDir, paste0(AIT_anndata$uns$taxonomyName, ".h5ad"))
      anndata_path = get_anndata_path(taxonomy_anndata_path)

      precomp_stats_output_path = get_precomp_stats(AIT_anndata, anndata_path, 
                                                    user_precomp_stats_path, temp_folder)

      query_markers_output_path = get_marker_genes(AIT_anndata, user_query_markers_path, temp_folder)

      query_data_output_path = get_query_data_path(query_data, temp_folder)

      get_mapmycells_results(query_data_output_path, extended_result_path, 
                             precomp_stats_output_path, query_markers_output_path, 
                             normalization, chunk_size, n_processors)

      ## Extract mapping results
      mapmycells_results_json = fromJSON(extended_result_path)

      ## Build mapping results dataframe
      mappingResults=list()
      for (hierarcy_level in names(mapmycells_results_json$results)) {
        if (hierarcy_level != "cell_id") {
          mappingResults[[paste0(hierarcy_level, "_Hierarchical")]] = 
              mapmycells_results_json$results[[hierarcy_level]]$assignment
          mappingResults[[paste0("bootstrapping_probability.Hierarchical.", hierarcy_level)]] = 
              mapmycells_results_json$results[[hierarcy_level]]$bootstrapping_probability
          mappingResults[[paste0("avg_correlation.Hierarchical.", hierarcy_level)]] = 
              mapmycells_results_json$results[[hierarcy_level]]$avg_correlation
        }
      }

      ## Return annotations and detailed model results
      return(mappingResults)
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
        file.remove(query_data_output_path)

        # remove the files is they were code generated
        if(!file.exists(taxonomy_anndata_path)) {
          file.remove(anndata_path)
        }
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

get_anndata_path = function(anndata_path) {
  # Check if anndata path exists; if does not, write it out to temp - show WARNING
  if (!file.exists(anndata_path)) {
    print(paste0(paste("WARNING: INVALID FILE PATH, ERROR in AIT.anndata$uns taxonomyDir and taxonomyName:", anndata_path),
          ". Writing the AIT.anndata to temperary location, SAVE anndata or FIX path to OPTIMIZE this step."))
    anndata_filename <- paste0(paste0("temp_anndata_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".h5ad")
    anndata_path <- file.path(temp_folder, anndata_filename)
    AIT_anndata$write_h5ad(anndata_path)
  }
  return(anndata_path)
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

get_query_data_path = function(query_data, temp_folder) {
  query_data_filename = paste0(paste0("query_data_", format(Sys.time(), "%Y%m%d-%H%M%S")), ".h5ad")
  query_data_output_path = file.path(temp_folder, query_data_filename)
  query_mat = Matrix::t(query_data)
  adata_query = AnnData(X = query_mat,
                        obs = data.frame(row.names = c(rownames(query_mat))),
                        var = data.frame(row.names = c(colnames(query_mat)))
                        )
  adata_query$write_h5ad(query_data_output_path)
  return(query_data_output_path)
}

get_mapmycells_results = function(query_data_output_path, extended_result_path, 
                precomp_stats_output_path, query_markers_output_path, 
                normalization, chunk_size, n_processors) {
  config <- list(
    'query_path' = query_data_output_path,
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