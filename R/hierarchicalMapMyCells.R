#' This function maps query data against the AIT (Shiny) reference taxonomy using hierarchical mapping.
#'
#' @param AIT_anndata A reference taxonomy anndata object.
#' @param query_data A logCPM normalized matrix to be annotated.
#' @param chunk_size Number of rows each worker process should load at a time from the query dataset.
#' @param n_processors Number of independent worker processes to spin up.
#' @param normalization Normalization of the h5ad files; must be either 'raw' or 'log2CPM'.
#' @param tmp_dir Temporary directory for writing out the hierarchical files.
#' @param user_extended_result_path Full file path and name where the original mapping results will be saved.
#' @param user_precomp_stats_path Alternative path to the user provided precompute stats HDF5 file. Will be generated, if not provided.
#' @param user_query_markers_path Alternative path to the user provided query markers JSON file. Will be generated, if not provided.
#' 
#' Note: this hierarchical mapping is a wrapper around cell_type_mapper,   
#'       and call's it's functions to generate needed files needed for mapping.
#' 
#' @import anndata
#'
#' @return List of mapping results with labels and scores.
#'
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

      taxonomy_anndata_path = file.path(AIT_anndata$uns$taxonomyDir, paste0(AIT_anndata$uns$title, ".h5ad"))
      anndata_path = get_anndata_path(taxonomy_anndata_path, temp_folder)

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

#' This function saves the AIT reference taxonomy to a temp folder as h5ad, if the provided file path is invalid.
#' @param anndata_path Local file path of the AIT reference taxonomy (h5ad file).
#' @return Local file path to the AIT reference taxonomy h5ad file.
#'
#' @keywords internal
get_anndata_path = function(anndata_path, temp_folder) {
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

#' Retrieves the precompute stats data from the AIT taxonomy and saves it in H5AD format. precomputed_stats.h5 file contains a serialization of the cell type taxonomy tree, 
#' as well as statistics about the reference cells that have been assigned to the cell types in the taxonomy. 
#'
#' @param AIT_anndata A reference taxonomy anndata object.
#' @param anndata_path Local file path of the AIT reference taxonomy (h5ad file).
#' @param precomp_stats_output_path Alternative path to the user provided precompute stats HDF5 file. Will be generated, if not provided.
#' @param temp_folder Temporary directory for writing out the hierarchical files (the code will clean these up after itself).
#'
#' @return File path to the precompute stats file.
#'
#' @keywords internal
get_precomp_stats = function(AIT_anndata, anndata_path, precomp_stats_output_path, temp_folder) {
  # retrieve precompute stats from anndata
  if(is.null(precomp_stats_output_path)) {
    uns_keys = list("hierarchical", AIT_anndata$uns$mode, "precomp_stats")
    precomp_stats_output_path <- toString(cell_type_mapper$utils$output_utils$uns_to_precomputed_stats(h5ad_path=anndata_path, 
                                      uns_keys_list=uns_keys, tmp_dir=temp_folder))
  }
  return(precomp_stats_output_path)
}

#' Retrieves the query markers from the AIT taxonomy and saves it in JSON format. The markers are, for every pair of leaf nodes in the cell type taxonomy tree, 
#' every gene that could conceivably be a marker gene for discriminating between those two cell types. 
#'
#' @param AIT_anndata A reference taxonomy anndata object.
#' @param query_markers_output_path Alternative path to the user provided query markers JSON file. Will be generated, if not provided.
#' @param temp_folder Temporary directory for writing out the hierarchical files (the code will clean these up after itself).
#'
#' @return File path of the query markers file.
#'
#' @keywords internal
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

#' Saves the query data count matrix as AnnData H5AD file. This file path is later passed into the hierarchical mapping algorithm.
#'
#' @param query_data A logCPM normalized matrix to be annotated.
#' @param temp_folder Temporary directory for writing out the hierarchical files (the code will clean these up after itself).
#'
#' @return File path of the query anndata file.
#'
#' @keywords internal
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

#' This function maps query data against the AIT (Shiny) reference taxonomy using hierarchical mapping.
#'
#' @param query_data_output_path Full file path to the query data, a logCPM normalized matrix to be annotated that's saved as AnnData H5AD file.
#' @param extended_result_path Full file path and name where the original mapping results will be saved.
#' @param precomp_stats_output_path Full file path to the precompute stats HDF5 file.
#' @param query_markers_output_path Full file path to the query markers JSON file.
#' @param chunk_size Number of rows each worker process should load at a time from the query dataset.
#' @param n_processors Number of independent worker processes to spin up.
#' @param normalization Normalization of the h5ad files; must be either 'raw' or 'log2CPM'.
#' 
#' @return 
#'
#' @keywords internal
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

#' Lists all the function parameters and their descriptions. 
#' @keywords internal
list_function_params = function() {
  command <- "python -m cell_type_mapper.cli.from_specified_markers --help"
  output <- system(command, intern = TRUE) 
  cat(output, sep = "\n") 
}