library(scrattch.mapping)
library(Seurat)

##
source("/home/nelson.johansen/scripts/R/scrattch-map_dev/scrattch-mapping/R/buildShiny.R")

## Load in the data to be annotated
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/pipeline/1001_C06_MTX-2046/1001_C06.RData")

## Counts input
counts = GetAssayData(rnaseq.data, "counts")

## meta.data input
rnaseq.data@meta.data$cluster = rnaseq.data$seurat_clusters
meta.data = rnaseq.data@meta.data
meta.data$sample_id = NULL

## feature.set input
feature.set = VariableFeatures(rnaseq.data)

## umap.coords input
umap.coords = FetchData(rnaseq.data, vars=c("UMAP_1","UMAP_2"))

##
buildReference(counts,
                meta.data,
                feature.set,
                umap.coords,
                shinyFolder="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/shinymap_test",
                subsample=2000)