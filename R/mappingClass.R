#' mappingClass S4 class definition
#'
#' Define the mappingClass that stores mapping results and detailed method outputs from scrattch.mapping
#'
#' @export
setClass(
  "mappingClass",
  representation(
    annotations = "data.frame",
    detailed_results = "list"
  )
)

#' Constructor function for mappingClass S4 class
#'
#' This function instantiates a mappingClass S4 class object.
#'
#' @param annotations A reference taxonomy anndata object.
#' @param detailed_results The number of cells to retain per cluster (default = 100).
#'
#' @examples
#' resultAnno <- mappingClass(
#'   annotations = data.frame(
#'     map.Corr = c("Exc", "Inh", "Inh"),
#'     score.Corr = c(0.9, 0.9, 0.9)
#'   ),
#'   detailed_results = list(
#'     corr = NA,
#'     tree = membership, ## For patchseq directory
#'     seurat = NA
#'   )
#' )
#'
#' @return Instance of mappingClass S4 class.
#'
#' @export
mappingClass <- function(annotations, detailed_results) {
  new("mappingClass", 
        annotations = annotations, 
        detailed_results = detailed_results)
}

#' Get cell type annotations
#'
#' Extract cell type annotations from mappingClass S4 class
#'
#' @return mapping results as a data.frame with labels in map.Method and scores in score.Method 
#'
#' @export
setGeneric("getMappingResults", 
  function(x) standardGeneric("getMappingResults")
)

#' Get cell type annotations
#'
#' Extract cell type annotations from mappingClass S4 class
#'
#' @return mapping results as a data.frame with labels in map.Method and scores in score.Method 
#'
#' @export
setMethod(
  "getMappingResults",
  signature = "mappingClass",
  function(object) {
    return(object@annotations)
  }
)
