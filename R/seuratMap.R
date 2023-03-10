#' Seurat based mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param query.data A logCPM normalized matrix to be annotated.
#'
#' @return Seurat mapping results as a data.frame.
#'
#' @export
seuratMap = function(AIT.anndata, query.data, dims=30, k.weight=5){
    print("Seurat-based mapping")

    ## Attempt Seurat mapping
    mappingTarget = tryCatch(
        expr = {

            ## Build Query Seruat object
            query.seurat = CreateSeuratObject(query.data[AIT.anndata$var_names[AIT.anndata$var$highly_variable_genes],])
            query.seurat = SetAssayData(query.seurat, slot = "data", new.data = query.data[AIT.anndata$var_names[AIT.anndata$var$highly_variable_genes],], assay = "RNA")
            
            ## Build Ref Seurat object
            ref.seurat = suppressWarnings(CreateSeuratObject(t(AIT.anndata$X[,AIT.anndata$var$highly_variable_genes]), meta.data=as.data.frame(AIT.anndata$obs)));
            ref.seurat = SetAssayData(ref.seurat, slot = "data", new.data = t(AIT.anndata$X[,AIT.anndata$var$highly_variable_genes]), assay = "RNA")
            
            ## Create a data list for label transfer
            seurat.list <- list(ref.seurat, query.seurat)
            names(seurat.list) <- c("Reference", "Query")

            ## Compute variable features for each object
            for (i in 1:length(x = seurat.list)) VariableFeatures(seurat.list[[i]]) <- AIT.anndata$var_names[AIT.anndata$var$highly_variable_genes]

            ## Seurat label transfer (celltype)
            Target.anchors <- FindTransferAnchors(reference = seurat.list[["Reference"]], query = seurat.list[["Query"]], 
                                                  dims = 1:dims, verbose=FALSE, npcs=dims)
            predictions    <- TransferData(anchorset = Target.anchors, refdata = seurat.list[["Reference"]]$cluster_label, 
                                                  dims = 1:dims, verbose=FALSE, k.weight=k.weight)
            ## Create results data.frame
            mappingTarget = data.frame(map.Seurat=as.character(predictions$predicted.id), 
                                       score.Seurat=predictions$prediction.score.max)
            mappingTarget
        },
        error = function(e){ 
            print("Error caught for Seurat mapping.")
            print(e)
            return(NULL)
        },
        finally = {
        }
    )
    return(mappingTarget)
}