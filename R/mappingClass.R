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
#' @param detailed_results a method-specific set of additional mapping output
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
  function(object, scores = TRUE) standardGeneric("getMappingResults")
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
  signature(object="mappingClass", scores="logical"),
  definition = function(object, scores = TRUE) {
    mapping.anno = object@annotations
    if (!scores) {
      score.cols = sapply(mapping.anno, is.numeric)
      mapping.anno = mapping.anno[, !score.cols]
    }
    return(mapping.anno)
  }
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
  signature(object="mappingClass", scores="missing"),
  definition = function(object) {
    # If scores not provided, set as TRUE by default.  Unclear why we need a separate solution for this.
    scores = TRUE
    mapping.anno = object@annotations
    if (!scores) {
      score.cols = sapply(mapping.anno, is.numeric)
      mapping.anno = mapping.anno[, !score.cols]
    }
    return(mapping.anno)
  }
)
