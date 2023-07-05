#' This function applies patchseqQC, given a taxonomy and query data
#'
#' @param AIT.anndata A reference taxonomy object.
#' @param query.data A count or CPM matrix for the query data.
#' @param query.metadata A data frame of metadata for the query data. 
#' @param verbose Should status be printed to the screen? 
#' 
#' @return A new query.metadata file with appended QC columns
#'
#' @export
applyPatchseqQC = function(AIT.anndata, 
                           query.data,
                           query.metadata,
                           verbose=FALSE){
  
  ## Checks
  if(!all(colnames(query.data) == rownames(query.metadata))){stop("Colnames of `query.data` and rownames of `meta.data` do not match.")}
  if(is.null(AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]])){stop(paste("QC_markers file must be provided for mode",AIT.anndata$uns$mode," in AIT.anndata in advance."))}
  
  ## Convert query.data to CPM
  if(is.element("data.frame",class(query.data))){stop("`query.data` should be a matrix or a sparse matrix, not a data.frame.")}
  if(max(query.data)<20){
    warning("`query.data` should not be log2-normalized. Converting back to linear space.")
    if (is.matrix(query.data)) {
      query.data <- 2^query.data - 1
    }
    else {
      query.data@x <- 2^query.data@x - 1
    }
  }
  query.cpm <- cpm(query.data)
  
  ## Load the reference files
  # ---- This includes markers, countsQC, cpmQC, classBr, subclassF, and allMarkers
  #load(AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]])
  allMarkers = AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]]$allMarkers 
  markers    = AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]]$markers
  countsQC   = AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]]$countsQC
  cpmQC      = AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]]$cpmQC
  classBr    = AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]]$classBr
  subclassF  = AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]]$subclassF
  rownames(cpmQC) <- rownames(countsQC) <- AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]]$qc_genes
  colnames(cpmQC) <- colnames(countsQC) <- AIT.anndata$uns$QC_markers[[AIT.anndata$uns$mode]]$qc_samples
  
  if(verbose) print("Format the reference and patch-seq data")
  ## -- NOTE: relevant reference data and type assignments are stored in refMarkerFile
  tmp                 <- cpmQC # countsQC
  rownames(tmp)       <- make.names(rownames(tmp))
  facs_df             <- as.data.frame(t(tmp[allMarkers,])+1)
  facs_df$sample_id   <- rownames(facs_df)
  facs_df$major_type  <- make.names(classBr)    
  facs_df$contam_type <- make.names(subclassF)  
  
  tmp              <- cpm(query.cpm) #
  rownames(tmp)    <- make.names(rownames(tmp))
  allMarkers       <- intersect(allMarkers,rownames(tmp))  # To account for differences in transcriptome
  pat_df           <- as.data.frame(t(tmp[allMarkers,rownames(query.metadata)])+1)
  pat_df$sample_id <- rownames(pat_df)
  
  ##
  if(verbose) print("Define which type each patch-seq cell is assigned to, based on maximal marker expression.")  
  nm          <- names(markers) <- make.names(names(markers))
  isOn        <- substr(nm,nchar(nm)-2,nchar(nm))=="_on"
  useThese    <- nm[isOn&(!is.element(nm,paste0(nm,"_on")))]
  subclassDat <- calcContamAllTypes(pat_df, markers[useThese])  # Identify subclass based on marker gene expression
  subclass    <- colnames(subclassDat)[subclassDat %>% apply(1,which.max)]
  subclass    <- gsub("_on","",subclass)
  
  pat_df$contam_type <- subclass  
  pat_df$major_type  <- make.names(classBr)[match(pat_df$contam_type,make.names(subclassF))]
  pat_df$contam_type <- paste0(pat_df$contam_type,"_on")
  
  ##
  if(verbose) print("Remove genes not in patch-seq transcriptome (ideally this wouldn't do anything)")
  for (i in 1:length(markers))
    markers[[i]] <- intersect(markers[[i]],allMarkers)
  
  if(verbose) print("Calculate the QC metrics")
  qcMetrics <- calculatePatchSeqQCMetrics2(pat_df,facs_df,markers)
  
  ##
  if(verbose) print("Set NMS>0.4 flag and determine most contaminated type")
  qcMetrics$Norm_Marker_Sum.0.4 <- c(TRUE,FALSE)[(qcMetrics$marker_sum_norm<0.40)+1]
  cls               <- make.names(sort(intersect(classBr,subclassF)))
  contaminationType <- cls[apply(qcMetrics[,cls],1,which.max)]
  qcMetrics$contaminationType   <- contaminationType
  
  ##
  if(verbose) print("Update annotations, including colors and ids for each new metadata.")
  cn      <- c("quality_score","marker_sum_norm","Norm_Marker_Sum.0.4","contaminationType","contam_sum")
  annoNew <- query.metadata
  for (i in 1:length(cn))      # Remove duplicate column names, if any
    annoNew <- annoNew[,!grepl(cn[i],colnames(annoNew))]
  annoNew <- cbind(annoNew,qcMetrics[,cn])
  annoNew$quality_score = pmax(annoNew$quality_score,0)
  annoNew$marker_sum_norm = pmax(annoNew$marker_sum_norm,0)
  annoNew$contam_sum      = pmax(annoNew$contam_sum,0)
  annoNew <- annoNew %>% annotate_num("quality_score", scale = "linear")    
  annoNew <- annoNew %>% annotate_num("marker_sum_norm", scale = "linear")  
  annoNew <- annoNew %>% annotate_num("contam_sum", scale = "linear") 
  annoNew <- annoNew %>% annotate_cat("Norm_Marker_Sum.0.4")
  annoNew <- annoNew %>% annotate_cat("contaminationType")
  
  if(verbose) print("Return the updated annotations")
  return(annoNew)
}