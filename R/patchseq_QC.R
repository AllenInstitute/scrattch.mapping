# #' Patchseq QC
# #'
# #' (previously: p3_patchseq_data_QC.R)
# #'
# #' @param refFolder Directory containing the Shiny taxonomy (full or patch-seq).
# #' @param patchseq_type The data storage method for patchseq: "Robj" or "feather"
# #' @param patchseq_dir
# #' @param patchseq_batch
# #' @param patchseq_robj_dir
# #' @param refMarkerFile QC_markers.RData
# #'
# #' @return Organized reference object ready for mapping against.
# #'
# #' @export
# patchseq_QC = function(refFolder, patchseq_type, patchseq_dir, patchseq_batch, patchseq_robj_dir, refMarkerFile){

#     ############
#     ## Read in the patch-seq data/metadata and perform the PatchseqQC
#     print("Read in and filter the patchseq data")
#     if (patchseq_type == "feather") {
#         ## Load in patchseq data and anno
#         query.meta   <- read_feather(file.path(patchseq_dir,"anno.feather")) 
#         query.data   <- feather(file.path(patchseq_dir,"data.feather"))
#         ## Ensure data and anno are in the same order
#         query.meta   <- query.meta[match(query.data$sample_id,query.meta$sample_id),] 
#         annoPat_all <- query.meta
#         datPat_all  <- as.matrix(query.data[,names(query.data)!="sample_id"])
#         rownames(datPat_all) <- annoPat_all$sample_id
#         datPat_all  <- t(datPat_all)
#     }else if(patchseq_type == "Robj"){
#         # Loading the query data, cpm & samp.dat
#         load(file.path(patchseq_robj_dir, paste0(patchseq_batch, "_mouse_patchseq_star2.0_samp.dat.Rdata")))
#         load(file.path(patchseq_robj_dir, paste0(patchseq_batch, "_mouse_patchseq_star2.0_cpm.Rdata")))

#         annoPat_all = samp.dat
#         annoPat_all$sample_id = as.character(samp.dat$patched_cell_container)
#         datPat_all = cpmR[,as.character(samp.dat$exp_component_name)]
#         colnames(datPat_all)=as.character(samp.dat$patched_cell_container)
#     }else{
#         stop("Unrecogognized patchseq data type argument, please provide either: 'Robj' or 'feather'")
#     }

#     print("Format the reference and patch-seq data")
#     ## -- NOTE: relevant reference data and type assignments are stored in refMarkerFile
#     tmp                 <- cpmQC # countsQC
#     rownames(tmp)       <- make.names(rownames(tmp))
#     facs_df             <- as.data.frame(t(tmp[allMarkers,])+1)
#     facs_df$sample_id   <- rownames(facs_df)
#     facs_df$major_type  <- as.character(classBr)     # Defined in previous section
#     facs_df$contam_type <- as.character(subclassF)   # Defined in previous section
#     facs_df             <- facs_df %>% filter(major_type != "Low Quality")

#     tmp              <- cpm(datPat_all) #
#     rownames(tmp)    <- make.names(rownames(tmp))
#     allMarkers       <- intersect(allMarkers,rownames(tmp))  # To account for differences in transcriptome
#     pat_df           <- as.data.frame(t(tmp[allMarkers,annoPat_all$sample_id])+1)
#     pat_df$sample_id <- rownames(pat_df)

# }

# #' Patchseq QC
# #'
# #' (previously: setup_patchseq_reference.R)
# #'
# #' @param GEXRef A reference taxonomy object.
# #' @param refFolder_pseq The data storage method for patchseq: "Robj" or "feather"
# #' @param qcCellsPerCluster
# #'
# #' @return Organized reference object ready for mapping against.
# #'
# #' @export

# c('annoReference', 'exprReference', 'datReference', 'dend', 'clustersUse', 'clusterInfo', 'kpSamp', 'referenceData', 'varFeatures')

# library(feather)
# library(dplyr)
# library(patchseqtools)  # devtools::install_github("AllenInstitute/patchseqtools")
# library(patchSeqQC)     # devtools::install_github('PavlidisLab/patchSeqQC')
# library(scrattch.io)    # devtools::install_github("AllenInstitute/scrattch"); scrattch::install_scrattch_deps(); scrattch::install_scrattch()
# library(scrattch.hicat)
# options(stringsAsFactors = FALSE)



# patchseq_qc_markers = function(GEXRef, refFolder_pseq){

#     print("Running marker gene setup for patchseqQC")
#     ## Update labeling for patchseqQC, for class we wrap neurons by class but retain non-neurons by subclass
#     classBr <- GEXRef$annoReference$subclass_label
#     classBr[GEXRef$annoReference$class_label!="Non-Neuronal"] = GEXRef$annoReference$class_label[GEXRef$annoReference$class_label!="Non-Neuronal"]
#     classBr   <- factor(classBr)
#     subclassF <- factor(GEXRef$annoReference$subclass_label)

#     ##
#     print("Define and output marker genes for each broad class and contamination class (use 50 for now)") 
#     if(file.exists(file.path(refFolder_pseq, "QC_markers.RData"))){
#       load(file.path(refFolder_pseq, "QC_markers.RData"))
#     }else{
#       print("Marker file does not exist for patchseq reference, generating now.")
#       markers    <- defineClassMarkers(GEXRef$datReference, subclassF, classBr, numMarkers = 50)
#       allMarkers <- unique(unlist(markers))
#     }

#     print("Format the reference and patch-seq data")
#     ## -- NOTE: relevant reference data and type assignments are stored in refMarkerFile
#     tmp                 <- cpmQC # countsQC
#     rownames(tmp)       <- make.names(rownames(tmp))
#     facs_df             <- as.data.frame(t(tmp[allMarkers,])+1)
#     facs_df$sample_id   <- rownames(facs_df)
#     facs_df$major_type  <- as.character(classBr)     # Defined in previous section
#     facs_df$contam_type <- as.character(subclassF)   # Defined in previous section
#     facs_df             <- facs_df %>% filter(major_type != "Low Quality")

#     tmp              <- cpm(datPat_all) #
#     rownames(tmp)    <- make.names(rownames(tmp))
#     allMarkers       <- intersect(allMarkers,rownames(tmp))  # To account for differences in transcriptome
#     pat_df           <- as.data.frame(t(tmp[allMarkers,annoPat_all$sample_id])+1)
#     pat_df$sample_id <- rownames(pat_df)


#     print("Define which type each patch-seq cell is assigned to, based on maximal marker expression.")  
#     nm          <- names(markers)
#     isOn        <- substr(nm, nchar(nm)-2, nchar(nm))=="_on"
#     useThese    <- nm[isOn&(!is.element(nm,paste0(nm,"_on")))]
#     subclassDat <- calcContamAllTypes(pat_df, markers[useThese])  # Identify subclass based on marker gene expression
#     subclass    <- colnames(subclassDat)[subclassDat %>% apply(1,which.max)]
#     subclass    <- gsub("_on","",subclass)

#     pat_df$contam_type <- subclass  
#     pat_df$major_type  <- as.character(classBr)[match(pat_df$contam_type,as.character(subclassF))]
#     pat_df$contam_type <- paste0(pat_df$contam_type,"_on")
#     pat_df             <- pat_df %>% filter(major_type != "Low Quality")


#     print("Remove genes not patch-seq transcriptome (ideally this wouldn't do anything)")
#     for (i in 1:length(markers))
#       markers[[i]] <- intersect(markers[[i]], allMarkers)
#     if ("Low Quality" %in% names(markers)) markers = markers[-match("Low Quality", names(markers))]  
#     print("Calculate the QC metrics")
#     qcMetrics <- calculatePatchSeqQCMetrics2(pat_df,facs_df,markers)


#     print("Set NMS>0.4 flag and determine most contaminated type")
#     qcMetrics$Norm_Marker_Sum.0.4 <- c(TRUE,FALSE)[(qcMetrics$marker_sum_norm<0.40)+1]
#     cls               <- sort(setdiff(classBr,c("GABAergic","Glutamatergic", "Low Quality")))
#     contaminationType <- cls[apply(qcMetrics[,cls],1,which.max)]
#     qcMetrics$contaminationType   <- contaminationType

#     save(markers, annoQC, countsQC, cpmQC, normQC, classBr, subclassF, allMarkers, file=file.path(refFolder_pseq, "QC_markers.RData"))

# }

#    refFolder_FACS    = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"
#    refFolder_pseq    = "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/Mapping_PatchSeq/PatchSeq_test/mouse_patchseq_VISp_20210412_collapsed40_cpm_test"
#    file.type         = "feather"
#    toneFN            = ""
#    cellsPerCluster   = 1000
#    completeDendFile  = "dend.RData"
#    qcCellsPerCluster = 100
#    refDendFile       = "dend.reference.RData"


# dir.create(patchseq_dir)
# load(file.path(refFolder_pseq, refMarkerFile))

# # Load libraries
# library(feather)
# library(dplyr)
# library(patchseqtools)  # devtools::install_github("AllenInstitute/patchseqtools")
# library(patchSeqQC)     # devtools::install_github('PavlidisLab/patchSeqQC')
# library(scrattch.io)    # devtools::install_github("AllenInstitute/scrattch"); scrattch::install_scrattch_deps(); scrattch::install_scrattch()
# library(scrattch.hicat)
# options(stringsAsFactors = FALSE)


# print("Format the reference and patch-seq data")
# ## -- NOTE: relevant reference data and type assignments are stored in refMarkerFile
# tmp                 <- cpmQC # countsQC
# rownames(tmp)       <- make.names(rownames(tmp))
# facs_df             <- as.data.frame(t(tmp[allMarkers,])+1)
# facs_df$sample_id   <- rownames(facs_df)
# facs_df$major_type  <- as.character(classBr)     # Defined in previous section
# facs_df$contam_type <- as.character(subclassF)   # Defined in previous section
# facs_df             <- facs_df %>% filter(major_type != "Low Quality")

# tmp              <- cpm(datPat_all) #
# rownames(tmp)    <- make.names(rownames(tmp))
# allMarkers       <- intersect(allMarkers,rownames(tmp))  # To account for differences in transcriptome
# pat_df           <- as.data.frame(t(tmp[allMarkers,annoPat_all$sample_id])+1)
# pat_df$sample_id <- rownames(pat_df)


# print("Define which type each patch-seq cell is assigned to, based on maximal marker expression.")  
# nm          <- names(markers)
# isOn        <- substr(nm,nchar(nm)-2,nchar(nm))=="_on"
# useThese    <- nm[isOn&(!is.element(nm,paste0(nm,"_on")))]
# subclassDat <- calcContamAllTypes(pat_df, markers[useThese])  # Identify subclass based on marker gene expression
# subclass    <- colnames(subclassDat)[subclassDat %>% apply(1,which.max)]
# subclass    <- gsub("_on","",subclass)

# pat_df$contam_type <- subclass  
# pat_df$major_type  <- as.character(classBr)[match(pat_df$contam_type,as.character(subclassF))]
# pat_df$contam_type <- paste0(pat_df$contam_type,"_on")
# pat_df             <- pat_df %>% filter(major_type != "Low Quality")


# print("Remove genes not patch-seq transcriptome (ideally this wouldn't do anything)")
# for (i in 1:length(markers))
#   markers[[i]] <- intersect(markers[[i]],allMarkers)
# if ("Low Quality" %in% names(markers)) markers = markers[-match("Low Quality", names(markers))]  
# print("Calculate the QC metrics")
# qcMetrics <- calculatePatchSeqQCMetrics2(pat_df,facs_df,markers)


# print("Set NMS>0.4 flag and determine most contaminated type")
# qcMetrics$Norm_Marker_Sum.0.4 <- c(TRUE,FALSE)[(qcMetrics$marker_sum_norm<0.40)+1]
# cls               <- sort(setdiff(classBr,c("GABAergic","Glutamatergic", "Low Quality")))
# contaminationType <- cls[apply(qcMetrics[,cls],1,which.max)]
# qcMetrics$contaminationType   <- contaminationType


# print("Update annotations")
# cn      <- c("quality_score","marker_sum_norm","Norm_Marker_Sum.0.4","contaminationType","contam_sum")
# annoNew <- annoPat_all[match(pat_df$sample_id, annoPat_all$exp_component_name),]
# for (i in 1:length(cn))      # Remove duplicate column names, if any
#   annoNew <- annoNew[,!grepl(cn[i],colnames(annoNew))]

# annoNew <- cbind(annoNew,qcMetrics[,cn])
# annoNew$quality_score = pmax(annoNew$quality_score,0)
# annoNew$marker_sum_norm = pmax(annoNew$marker_sum_norm,0)
# annoNew$contam_sum      = pmax(annoNew$contam_sum,0)
# annoNew <- annoNew %>% annotate_num("quality_score", scale = "linear")    
# annoNew <- annoNew %>% annotate_num("marker_sum_norm", scale = "linear")  
# annoNew <- annoNew %>% annotate_num("contam_sum", scale = "linear") 
# annoNew <- annoNew %>% annotate_cat("Norm_Marker_Sum.0.4")
# annoNew <- annoNew %>% annotate_cat("contaminationType")

# print("Output an updated anno.feather file. ")
# if (file.exists(file.path(patchseq_dir,"anno.feather"))){ 
#   file.remove(file.path(patchseq_dir,"anno.feather"))
# }
# write_feather(annoNew,file.path(patchseq_dir,"anno.feather")) 




# # input arguments
# args <- commandArgs(TRUE)
# #if (0) {
# #   args[1]= "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/Mapping_PatchSeq/Ref_test/mouse_patchseq_VISp_20210412_collapsed40_cpm"
# #   args[2]= "Robj"
# #   args[3]= "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/Mapping_PatchSeq/PatchSeq_test/mouse_patchseq_VISp_20210412_collapsed40_cpm"
# #   args[4]= "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
# #   args[5]= "20210331_BT014-RSC-273"
# #} else {
# #   args[1]="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_MOp_miniatlas_newtest_hpc/"
# #   args[2]=patchseq_type="Robj"
# #   args[3]=patchseq_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/Mapping_PatchSeq/Test_qc/"
# #   args[4]=patchseq_robj_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
# #   args[5]=patchseq_batch="20210324_BT014-RSC-272"
# #}
# refFolder_pseq = args[1]
# patchseq_type  = args[2] # "Robj", "feather"
# if (patchseq_type == "feather") {
#    patchseq_dir   = args[3]
#    patchseq_batch = args[4]
# } else {
#    patchseq_dir        = args[3]
#    patchseq_robj_dir   = args[4]
#    patchseq_batch      = args[5]
#    dir.create(patchseq_dir)
# }
