#' Seurat based mapping
#'
#' @param AIT.anndata A reference taxonomy anndata object.
#' @param query.data A logCPM normalized matrix to be annotated.
#' @param dims Number of principle component dimensions to use for FindTransferAnchors and TransferData (default = 30)
#' @param k.weight k.weight parameter for TransferData Seurat function (default = 5)
#' @param genes.to.use The set of genes to use for Seurat variable genes (default is the highly_variable_genes associated with the current mode). Can either be a character vector of gene names or a TRUE/FALSE (logical) vector of which genes to include.
#' @param cells.to.use The set of cells to include in Seurat integration (default is the filtered cell set associated with the current mode). Can either be a character vector of cell names or a TRUE/FALSE (logical) vector of which cells to include.
#'
#' @return Seurat mapping results as a data.frame.
#'
#' @export
seuratMap = function(AIT.anndata, 
                     query.data, 
                     dims=30, 
                     k.weight=5,
                     genes.to.use=NULL,
                     cells.to.use=NULL){
    print("Seurat-based mapping")

    ## Attempt Seurat mapping
    mappingTarget = tryCatch(
        expr = {
            ## Create the vector of genes to use.
            if(is.null(genes.to.use)){
              # If null, default to correct set of highly variable genes
              genes.to.use.vector <- AIT.anndata$var[,paste0("highly_variable_genes_",AIT.anndata$uns$mode)]
            } else if (is.logical(genes.to.use)) {
              if (length(genes.to.use)!=dim(AIT.anndata)[2]) stop("If genes.to.use is logical it must be the same length as the total number of genes in AIT.anndata.")
              genes.to.use.vector = genes.to.use
            } else if (is.character(genes.to.use)){
              genes.to.use <- intersect(genes.to.use,AIT.anndata$var_names)
              if (length(genes.to.use)==0) stop("No valid gene names provided in genes.to.use.")
              if (length(genes.to.use)<=2) stop("More than 2 valid gene names must be provided in genes.to.use.")
              genes.to.use.vector <- is.element(AIT.anndata$var_names,genes.to.use)
            } else {
              stop("genes.to.use must be a character or logical vector.")
            }
          
            ## Create the vector of cells to use.
            if(is.null(cells.to.use)){
              # If null, default to correct set of highly variable genes
              cells.to.use.vector <- !AIT.anndata$uns$filter[[AIT.anndata$uns$mode]]
            } else if (is.logical(cells.to.use)) {
              if (length(cells.to.use)!=dim(AIT.anndata)[1]) stop("If cells.to.use is logical it must be the same length as the total number of cells in AIT.anndata.")
              cells.to.use.vector = cells.to.use
            } else if (is.character(cells.to.use)){
              cells.to.use <- intersect(cells.to.use,AIT.anndata$obs_names)
              if (length(cells.to.use)==0) stop("No valid cell names provided in cells.to.use")
              if (length(cells.to.use)<=2) stop("More than 2 valid cell names must be provided in cells.to.use.")
              cells.to.use.vector <- is.element(AIT.anndata$obs_names,cells.to.use)
            } else {
              stop("cells.to.use must be a character or logical vector.")
            }

            ## Build Query Seurat object
            use.genes    = intersect(rownames(query.data),AIT.anndata$var_names[genes.to.use.vector])
            query.seurat = suppressWarnings(CreateSeuratObject(query.data[use.genes,]))
            query.seurat = suppressWarnings(SetAssayData(query.seurat, slot = "data", new.data = query.data[use.genes,], assay = "RNA"))
            VariableFeatures(query.seurat) <- use.genes
            
            ## Build Ref Seurat object
            ref.data   = as.matrix(BiocGenerics::t(AIT.anndata$X[,use.genes]))
            ref.data   = ref.data[,cells.to.use.vector]# UPDATE to filter cells
            ref.seurat = suppressWarnings(CreateSeuratObject(ref.data, meta.data=as.data.frame(AIT.anndata$obs)[cells.to.use.vector,]));
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