#' Seurat based mapping
#'
#' @param GEXRef A reference taxonomy object.
#' @param query.data A logCPM normalized matrix to be annotated.
#'
#' @return Seurat mapping results as a data.frame.
#'
#' @export
seuratMap = function(GEXRef, query.data, dims=30, k.weight=5){
    print("Seurat-based mapping")
    ## Build Query Seruat object
    query.seurat = CreateSeuratObject(query.data[GEXRef$varFeatures,])
    query.seurat = SetAssayData(query.seurat, slot = "data", new.data = query.data[GEXRef$varFeatures,], assay = "RNA")
    ## Attempt Seurat mapping
    mappingTarget = tryCatch(
        expr = {
            ## Create a data list for label transfer
            seurat.list <- list(GEXRef$referenceData[GEXRef$varFeatures,], query.seurat)
            names(seurat.list) <- c("Reference", "Query")

            ## Compute variable features for each object
            for (i in 1:length(x = seurat.list)) VariableFeatures(seurat.list[[i]]) <- GEXRef$varFeatures

            ## Seurat label transfer (celltype)
            Target.anchors <- FindTransferAnchors(reference = seurat.list[["Reference"]], query = seurat.list[["Query"]], 
                                                  dims = 1:dims, verbose=FALSE, npcs=dims)
            predictions    <- TransferData(anchorset = Target.anchors, refdata = seurat.list[["Reference"]]$cluster_label, 
                                                  dims = 1:dims, verbose=FALSE, k.weight=k.weight)
            ## Create results data.frame
            mappingTarget = data.frame(map.Tree=as.character(predictions$predicted.id), 
                                       score.Tree=predictions$prediction.score.max)
            mappingTarget
        },
        error = function(e){ 
            print("Error caught for Seurat mapping.")
            print(e)
            return(NULL)
        },
        warning = function(w){
        },
        finally = {
        }
    )
    return(mappingTarget)
}