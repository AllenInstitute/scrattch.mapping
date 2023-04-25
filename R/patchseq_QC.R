#' Save marker genes for patchSeqQC
#'
#' This function write a file called `QC_markers.RData` that contains all the variables required for applying the patchseq QC algorithm `pathseqtools` (which is an more flexible version of the `patchSeqQC` algorithm).  This is only used for patch-seq analysis.
# ----- Subclass calls for each cell
# ----- Broad class class calls for each cell
# ----- Distinction of neuron (e.g., mappable type) vs. non-neuron (e.g., contamination type)
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param subsample The number of cells to retain per cluster (default = 100).
#' @param subclass.column Column name corresponding to the moderate-resolution cell types used for the cell types of interest (default = "subclass_label").
#' @param class.column Column name corresponding to the low-resolution cell types used for the off-target cell types (default = "class_label").
#' @param off.target.types A character vector of off-target (also known as 'contamination') cell types.  This must include at least one of the cell types found in "class.column" and/or "subclass.column" (both columns are checked)
#' @param num.markers The maximum number of markers to calculate per node per direction (default = 50)
#' @param shinyFolder = The location to save shiny output (default = current working directory).
#' 
#' Nothing is returned; however an R data object called "QC_markers.RData" is returned with the following variables
#' markers, 
#' countsQC, 
#' cpmQC, 
#' classBr, 
#' subclassF, 
#' allMarkers, 
#' 
#' @import patchseqtools
#' @import scrattch.hicat
#'
#' @return AIT.anndata
#'
#' @export
buildPatchseqTaxonomy = function(AIT.anndata,
                                   subsample = 100,
                                   subclass.column = "subclass_label",
                                   class.column = "class_label",
                                   off.target.types = c("Glia","glia","non-neuronal","Non-neuronal"),
                                   num.markers = 50,
                                   shinyFolder = paste0(getwd(),"/")
){
  
  ## Required variables
  counts = Matrix::t(AIT.anndata$layers["counts"])
  metadata = AIT.anndata$obs

  ## Cell type labels
  if(!is.element(subclass.column, colnames(metadata))){stop(paste(subclass.column,"is not a column in the metadata data frame."))}
  if(!is.element(class.column, colnames(metadata))){stop(paste(class.column,"is not a column in the metadata data frame."))}
  if(!dir.exists(shinyFolder)){"Specified taxonomy folder does not exist."}

  ## Ensure variable naming scheme matches assumptions
  metadata$subclass_label = metadata[,subclass.column]  # For compatibility with existing code.
  metadata$class_label = metadata[,class.column]  # For compatibility with existing code.
  
  ## Subsample and filter metadata and data
  kpSamp2  = subsampleCells(metadata$subclass_label, subsample)
  goodSamp = !is.na(metadata$class_label)  # For back-compatibility; usually not used
  kpSamp2  = kpSamp2&goodSamp              # For back-compatibility; usually not used
  annoQC   = metadata[kpSamp2,]
  annoQC$subclass_label = make.names(annoQC$subclass_label)
  datQC    = as.matrix(counts[,kpSamp2])

  ## Define class and subclass labels
  ## --- We wrap on-target types by class but retain off-target types by subclass
  offTarget = is.element(annoQC$class_label, off.target.types) | is.element(annoQC$subclass_label, off.target.types)
  if(sum(offTarget)==0){stop("No valid off-target classes or subclasses are provided. Please update off.target.types accordingly.")}
  
  classBr   = annoQC$subclass_label
  classBr[!offTarget] = annoQC$class_label[!offTarget]
  classBr   = factor(classBr)
  subclassF = factor(annoQC$subclass_label)
  
  print("Define and output marker genes for each broad class and off-target subclass.") 
  ## -- These are selected using some reasonable approach that could probably be improved, if needed.    
  markers    = defineClassMarkers(datQC, subclassF, classBr, numMarkers = 50)
  allMarkers = unique(unlist(markers))
  rownames(datQC) = make.names(rownames(datQC))
  countsQC   = datQC[allMarkers,]
  cpmQC      = cpm(datQC)[allMarkers,]  ## Only use of scrattch.hicat in this function

  ##
  save(markers, countsQC, cpmQC, classBr, subclassF, allMarkers, file=file.path(shinyFolder,"QC_markers.RData"))

  ## Identify offtarget cells to filter out.
  AIT.anndata$uns$offTarget.patchseq = is.element(metadata$class_label, off.target.types) | is.element(metadata$subclass_label, off.target.types)
  AIT.anndata$uns$QC_markers = file.path(shinyFolder,"QC_markers.RData")
  AIT.anndata$write_h5ad(file.path(shinyFolder,"AI_taxonomy_patchseq.h5ad"))

  ##
  return(AIT.anndata)
}
