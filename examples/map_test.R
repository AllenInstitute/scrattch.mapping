## Test that should be run after "build_test.R"

library(scrattch.mapping)
library(Seurat)
library(anndata)

## Load in the data to be annotated
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/pipeline/1001_C06_MTX-2046/1001_C06.RData")
query.data = as.matrix(GetAssayData(rnaseq.data, "data"))

## Specify which taxonomies to map against.
taxonomies = c("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/shinymap_test")

taxonomy = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/human_BGplus_20211020"

## Map data against all taxonomies
mapping.anno = list()
for(taxonomy in taxonomies){
    ## Load the shiny taxonomy into a standard object for mapping.
    AIT.anndata = loadTaxonomy(refFolder = taxonomy)

    ## Map! Returns a data.frame with mapping results.
    mapping.anno[[taxonomy]] = taxonomy_mapping(AIT.anndata=AIT.anndata, seurat.map=FALSE, query.data=query.data)
}

##
save(mapping.anno, file="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/pipeline/1001_C06_MTX-2046/mapping_anno.RData")