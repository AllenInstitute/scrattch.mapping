library(scrattch.mapping)
library(Seurat)

## Load in the data to be annotated
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/pipeline/1001_C06_MTX-2046/1001_C06.RData")

##
counts = GetAssayData(rnaseq.data, "counts")
metadata = rnaseq.data@meta.data

##
seurat.obj = createSeuratObjectForReferenceFolder(counts, 
                                     metadata)

## Additional tests required for various parameter options
# seurat.obj = createSeuratObjectForReferenceFolder(counts, 
#                                      metadata)

##
buildReferenceFolder(seurat.obj,
                    "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/shinymap_test",
                    subsample=2000,
                    feature.set=NULL,
                    save.normalized.data=FALSE)