## Test that should be run after "build_test.R"

library(scrattch.mapping)
library(Seurat)

## Load in the data to be annotated
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/pipeline/1001_C06_MTX-2046/1001_C06.RData")
query.data = as.matrix(GetAssayData(rnaseq.data, "data"))

## Specify which taxonomies to map against.
taxonomies = c("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/shinymap_test")

## Map data against all taxonomies
mapping.anno = list()
for(taxonomy in taxonomies){
    ## Load the shiny taxonomy into a standard object for mapping.
    GEXRef = loadGEXRef(refFolder = taxonomy)

    ## Map! Returns a data.frame with mapping results.
    mapping.anno[[taxonomy]] = taxonomy_mapping(GEXRef=GEXRef, query.data=query.data, label.cols = c("cluster_label"))
}

##
save(mapping.anno, file="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/pipeline/1001_C06_MTX-2046/mapping_anno.RData")