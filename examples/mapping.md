# Tutorial: Standard taxonomy mapping using flat, tree and Seurat.

In this tutorial we demonstrate how to run standard mapping algorithms using scrattch.mapping on the Tasic et al. 2016 study. Available taxonomies can be found under `/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/` on hpc.

```R
## Load scrattch.mapping
library(scrattch.mapping)
library(scrattch.taxonomy)

## Load in example count data
library(tasic2016data)

## Compute log CPM
query.data = tasic_2016_counts
query.data = logCPM(query.data)

## Specify which taxonomies to map against.
taxonomyDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016"

## Load the shiny taxonomy into a standard object for mapping.
AIT.anndata = loadTaxonomy(taxonomyDir = taxonomyDir, anndata_file="Tasic2016.h5ad")

## Map! Returns an S4 class with mapping results.
mapping.anno = taxonomy_mapping(AIT.anndata=AIT.anndata,
                                query.data=query.data,
                                label.cols="cluster_label", ## Which obs in AIT.anndata contain annotations to map. E.g. "class", "subclass", etc.
                                corr.map=TRUE,
                                tree.map=TRUE,
                                seurat.map=FALSE)

## Extract mapping results from S4 mappingClass
mapping.results = getMappingResults(mapping.anno)

## Extract tree mapping bootstraping table (We will improve this in the near future.)
tree.bootstraps = mapping.anno@detailed_results[["tree"]]

## Save
save(mapping.results, file="mapping_results.rda")
```
