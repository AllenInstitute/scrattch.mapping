#' Save marker genes for patchSeqQC
#'
#' This function saves all the variables required for applying the patchseq QC algorithm `pathseqtools` (which is an more flexible version of the `patchSeqQC` algorithm) to AIT.anndata$uns. This is only used for patch-seq analysis.  Requirements for input include:
# ----- Subclass calls for each cell
# ----- Broad class class calls for each cell
# ----- Distinction of neuron (e.g., mappable type) vs. non-neuron (e.g., contamination type)
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param subsample The number of cells to retain per cluster (default = 100).
#' @param subclass.column Column name corresponding to the moderate-resolution cell types used for the cell types of interest (default = "subclass_label").
#' @param class.column Column name corresponding to the low-resolution cell types used for the off-target cell types (default = "class_label").
#' @param off.target.types A character vector of off-target (also known as 'contamination') cell types.  This must include at least one of the cell types found in "class.column" and/or "subclass.column" (both columns are checked)
#' @param mode.name A name to identify the new taxonomy version.
#' @param num.markers The maximum number of markers to calculate per node per direction (default = 50)
#' @param taxonomyDir = The location to save shiny output (default = current working directory).
#' 
#' The following variables are added to AIT.anndata$uns
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
#' @return AIT.anndata An updated AIT.anndata variable with the above content added to AIT.anndata$uns for the relevant mode.name.
#'
#' @export
buildPatchseqTaxonomy = function(AIT.anndata,
                                 mode.name = "patchseq", ## "Inhibitory"
                                 subsample = 100,
                                 subclass.column = "subclass_label",
                                 class.column = "class_label",
                                 off.target.types = c("Glia","glia","non-neuronal","Non-neuronal"), ## "Gluta", "NN"
                                 num.markers = 50,
                                 taxonomyDir = file.path(AIT.anndata$uns$taxonomyDir)
){

  ## Ensure filtering mode doesn't already exist
  if(mode.name %in% names(AIT.anndata$uns$filter)){ print(paste0("Print ", mode.name, " already in Taxonomy, you will be overwriting the previous mode files.")) }

  ## Create the required files for patchSeqQC and determine offtarget cells
  if(!is.element("counts", names(AIT.anndata$layers))){stop("`counts` must exist in AIT.anndata$layers, check taxonomy.")}
  if(!is.element(subclass.column, colnames(AIT.anndata$obs))){stop(paste(subclass.column,"is not a column in the metadata data frame."))}
  if(!is.element(class.column, colnames(AIT.anndata$obs))){stop(paste(class.column,"is not a column in the metadata data frame."))}
  if(!dir.exists(file.path(taxonomyDir))){"Specified taxonomy folder does not exist."}

  ## Determine taxonomy mode directory (Move to utilty function)
  if(mode.name == "standard"){ taxonomyModeDir = file.path(taxonomyDir) } else { taxonomyModeDir = file.path(file.path(taxonomyDir), mode.name) }
  if(!dir.exists(taxonomyModeDir)){  dir.create(taxonomyModeDir, showWarnings = FALSE) }

  ## Copy metadata
  metadata = AIT.anndata$obs

  ## Ensure variable naming scheme matches assumptions
  metadata$subclass_label = AIT.anndata$obs[,subclass.column]  # For compatibility with existing code.
  metadata$class_label = AIT.anndata$obs[,class.column]  # For compatibility with existing code.
  
  ## Subsample and filter metadata and data
  kpSamp2  = subsampleCells(metadata$subclass_label, subsample)
  goodSamp = !is.na(metadata$class_label)  # For back-compatibility; usually not used
  kpSamp2  = kpSamp2 & goodSamp            # For back-compatibility; usually not used
  annoQC   = metadata[kpSamp2,]
  annoQC$subclass_label = make.names(annoQC$subclass_label)
  datQC    = as.matrix(Matrix::t(AIT.anndata$layers["counts"])[,kpSamp2])

  ## Define class and subclass labels
  ## --- We wrap on-target types by class but retain off-target types by subclass
  offTarget = is.element(annoQC$class_label, off.target.types) | is.element(annoQC$subclass_label, off.target.types)
  if(sum(offTarget)==0){stop("No valid off-target classes or subclasses are provided. Please update off.target.types accordingly.")}
  
  ## 
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

  ## Identify off target cells to filter out.
  AIT.anndata$uns$filter[[mode.name]] = is.element(metadata$class_label, off.target.types) | is.element(metadata$subclass_label, off.target.types)
  
  ## Save patchseqQC information to uns
  AIT.anndata$uns$QC_markers[[mode.name]]$allMarkers = allMarkers
  AIT.anndata$uns$QC_markers[[mode.name]]$markers    = markers
  AIT.anndata$uns$QC_markers[[mode.name]]$countsQC   = countsQC
  AIT.anndata$uns$QC_markers[[mode.name]]$cpmQC      = cpmQC
  AIT.anndata$uns$QC_markers[[mode.name]]$classBr    = classBr
  AIT.anndata$uns$QC_markers[[mode.name]]$subclassF  = subclassF
  AIT.anndata$uns$QC_markers[[mode.name]]$qc_samples = colnames(countsQC) # since colnames are lost
  AIT.anndata$uns$QC_markers[[mode.name]]$qc_genes   = rownames(countsQC) # since rownames are lost
  
  ##################
  ## ------- Modify the dendrogram and save
  ##

  ## Load the complete dendrogram
  dend = readRDS(file.path(AIT.anndata$uns$dend[["standard"]]))

  ## Prune dendrogram to remove off.target types
  dend = prune(dend, setdiff(labels(dend), unique(AIT.anndata$obs$cluster_label[!AIT.anndata$uns$filter[[mode.name]]])))

  ## Save dendrogram
  saveRDS(dend, file.path(taxonomyModeDir,"dend.RData"))

  ## Store the pruned dendrogram, in dend list under "patchseq" mode.name
  AIT.anndata$uns$dend[[mode.name]] = file.path(taxonomyModeDir,"dend.RData")

  ## Save patch-seq taxonomy anndata
  AIT.anndata$write_h5ad(file.path(taxonomyDir, "AI_taxonomy.h5ad"))

  ## Update markers after pruning
  dend = addDendrogramMarkers(AIT.anndata, mode=mode.name)

  ##
  return(AIT.anndata)
}