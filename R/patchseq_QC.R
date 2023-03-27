#' Save marker genes for patchSeqQC
#'
#' This function write a file called `QC_markers.RData` that contains all the variables required for applying the patchseq QC algorithm `pathseqtools` (which is an more flexible version of the `patchSeqQC` algorithm).  This is only used for patch-seq analysis.
# ----- Subclass calls for each cell
# ----- Broad class class calls for each cell
# ----- Distinction of neuron (e.g., mappable type) vs. non-neuron (e.g., contamination type)
#'
#' @param counts A matrix of counts for the reference data, or character string with a file location of "counts.feather".  If it appears cpm or logCPM matrix is provided, an warning will be thrown.
#' @param metadata Data frame of metadata with rows corresponding to cells/nuclei, and either row names or a column called "sample_id" corresponding to cell names. This matrix must include entries for all cells in norm.data.  Could also be a file.  Columns can be numeric, categorical, or factors and must include "subclass.column" and "class.column"
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
#'
#' @export
writePatchseqQCmarkers = function(counts,
                                   metadata,
                                   subsample = 100,
                                   subclass.column = "subclass_label",
                                   class.column = "class_label",
                                   off.target.types = c("Glia","glia","non-neuronal","Non-neuronal"),
                                   num.markers = 50,
                                   shinyFolder = paste0(getwd(),"/")
){
  
  ## Checks and data formatting
  #counts
  if(class(counts)=="character"){
    if(file.exists(counts)){
      data_t <- try({feather(counts)})
      if(class(data_t)=="try-error"){stop("counts must be in feather format with a gene column + data matrix included.")}
      counts <- as.matrix(data_t[,colnames(data_t)!="gene"])
      if(!is.element(class(counts[,1]),c("numeric","integer"))){stop("counts doesn't have numeric entries")}
      genes <- try(data_t$gene)
      if(class(genes)=="try-error"){stop("counts must have a column called gene with gene names.")}
      rownames(counts) <- genes
    } else {
      stop(paste(counts,"is not a valid filename."))
    }
  } else {
    if(!grepl("atrix",as.character(class(counts)))){stop("counts must be some kind of matrix or a file name.")}
  }
  # Check if values are numeric and likely to be counts
  if(!is.element(class(counts[,1]),c("numeric","integer"))){stop("norm.data doesn't have numeric entries")}
  if(max(counts)<50) warning("Maximum value of counts matrix is <50.  Please check counts matrix is not log normalized.")
  if(abs(sum(counts[,1])-1000000)<1) warning("Maximum value of counts sums to 1e6.  Please check counts and not CPM was input.")
  
  #metadata 
  if(class(metadata)=="character"){
    if(file.exists(metadata)){
      metadata <- try({feather(metadata)})
      if(class(data_t)=="try-error"){stop("metadata file must be in feather format.")}
      metadata <- as.data.frame(metadata)
    } else {
      stop(paste(metadata,"is not a valid filename."))
    }
  } else {
    if(class(metadata)!="data.frame"){stop("metadata must be a data frame or a file name.")}
  }
  if(is.element("sample_id",colnames(metadata))) rownames(metadata) <- metadata$sample_id
  if(length(intersect(rownames(metadata),colnames(counts)))==0){stop("Metadata sample ids (or row names) do not match sample ids in the data (column names).")}
  if(length(setdiff(colnames(counts),rownames(metadata)))>0){stop("Some samples are missing in the metadata file based on matching sample ids between counts columns and metadata rows.")}
  metadata <- metadata[colnames(counts),]
  
  #celltype labels
  if(!is.element(subclass.column,colnames(metadata))){stop(paste(subclass.column,"is not a column in the metadata data frame."))}
  if(!is.element(class.column,colnames(metadata))){stop(paste(class.column,"is not a column in the metadata data frame."))}
  metadata$subclass_label = metadata[,subclass.column]  # For compatibility with existing code.
  metadata$class_label = metadata[,class.column]  # For compatibility with existing code.
  
  # Ensure directory exists, if not create it
  dir.create(shinyFolder, showWarnings = FALSE)
  
  
  # Subsample and filter metadata and data
  kpSamp2  <- subsampleCells(metadata$subclass_label, subsample)
  goodSamp <- !is.na(metadata$class_label)  # For back-compatibility; usually not used
  kpSamp2  <- kpSamp2&goodSamp              # For back-compatibility; usually not used
  annoQC   <- metadata[kpSamp2,]
  annoQC$subclass_label = make.names(annoQC$subclass_label)
  datQC    <- as.matrix(counts[,kpSamp2])

  # Define class and subclass labels
  # --- We wrap on-target types by class but retain off-target types by subclass
  offTarget <- is.element(annoQC$class_label,off.target.types)|is.element(annoQC$subclass_label,off.target.types)
  if(sum(offTarget)==0){stop("No valid off-target classes or subclasses are provided. Please update off.target.types accordingly.")}
  
  classBr   <- annoQC$subclass_label
  classBr[!offTarget] = annoQC$class_label[!offTarget]
  classBr   <- factor(classBr)
  subclassF <- factor(annoQC$subclass_label)


  ## This is the main code here
  
  print("Define and output marker genes for each broad class and off-target subclass.") 
  # -- These are selected using some reasonable approach that could probably be improved, if needed.    
  markers    <- defineClassMarkers(datQC,subclassF,classBr,numMarkers = 50)
  allMarkers <- unique(unlist(markers))
  rownames(datQC) <- make.names(rownames(datQC))
  countsQC   <- datQC[allMarkers,]
  cpmQC      <- cpm(datQC)[allMarkers,]  # Only use of scrattch.hicat in this function

  save(markers, countsQC, cpmQC, classBr, subclassF, allMarkers, file=paste0(shinyFolder,"QC_markers.RData"))
  print("Complete!")

}
