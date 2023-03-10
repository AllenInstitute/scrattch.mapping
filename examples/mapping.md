# Tutorial: Standard taxonomy mapping using flat, tree and Seurat.

In this tutorial we demonstrate how to run standard mapping algorithms using scrattch.mapping on the Tasic et al. 2016 study. Available taxonomies can be found under `/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/` on hpc.

```R
## Load scrattch.mapping
library(scrattch.mapping)

## Load in example count data
library(tasic2016data)

## Compute log CPM
query.data = tasic_2016_counts
query.data = logCPM(query.data)

## Specify which taxonomies to map against.
taxonomies = c("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016")

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
```