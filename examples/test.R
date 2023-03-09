## Load scrattch.mapping
library(scrattch.mapping)

## Load in example count data
library(tasic2016data)

## Compute log CPM
query.data = tasic_2016_counts
query.data = logCPM(query.data)

## Specify which taxonomies to map against.
taxonomies = c("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/human_BGplus_20211020")

## Map data against all taxonomies
mapping.anno = list()
for(taxonomy in taxonomies){
    ## Load the shiny taxonomy into a standard object for mapping.
    AIT.anndata = loadTaxonomy(refFolder = taxonomy)
    ## Map! Returns a data.frame with mapping results.
    mapping.anno[[taxonomy]] = taxonomy_mapping(AIT.anndata=AIT.anndata, 
                                                query.data=query.data, 
                                                label.cols="cluster_label",
                                                corr.map=TRUE, 
                                                tree.map=TRUE, 
                                                seurat.map=TRUE)
}