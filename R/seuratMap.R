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
            use.genes    = intersect(rownames(query.data),AIT.anndata$var_names[AIT.anndata$var$highly_variable_genes])
            query.seurat = suppressWarnings(CreateSeuratObject(query.data[use.genes,]))
            query.seurat = suppressWarnings(SetAssayData(query.seurat, slot = "data", new.data = query.data[use.genes,], assay = "RNA"))
            VariableFeatures(query.seurat) <- use.genes
            
            ## Build Ref Seurat object
            ref.data   = as.matrix(BiocGenerics::t(AIT.anndata$X[,use.genes]))
            ref.seurat = suppressWarnings(CreateSeuratObject(ref.data, meta.data=as.data.frame(AIT.anndata$obs)));
            ref.seurat = suppressWarnings(SetAssayData(ref.seurat, slot = "data", new.data = ref.data, assay = "RNA"))
            VariableFeatures(ref.seurat) <- use.genes

            ## Seurat label transfer (celltype)
            Target.anchors <- suppressWarnings(FindTransferAnchors(reference = ref.seurat, query = query.seurat, 
                                                                   dims = 1:dims, verbose=FALSE, npcs=dims))
            predictions    <- suppressWarnings(TransferData(anchorset = Target.anchors, refdata = ref.seurat$cluster_label, 
                                                  dims = 1:dims, verbose=FALSE, k.weight=k.weight))
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