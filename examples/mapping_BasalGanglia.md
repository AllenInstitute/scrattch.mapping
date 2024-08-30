# Tutorial: Standard taxonomy mapping using flat and hierarchical.

In this tutorial we demonstrate how to run multiple mapping algorithms using scrattch.mapping on the HMBA Macaque consensus basal ganglia taxonomy (applies to Human as well).

To run this code please use the docker provided at: `docker://njjai/scrattch_mapping:0.6.3`:

```R
library(scrattch.taxonomy)
library(scrattch.mapping)
library(reticulate)
cell_type_mapper <- import("cell_type_mapper") ## For now this has to be defined for hierarchical.map to work.

## ------------------------------------------------------
## Build Shiny taxonomy
AIT.anndata = loadTaxonomy(taxonomyDir="/PATH/TO/AIT",
                            anndata_file="HMBA_Macaque_BG_082024_AIT.h5ad")

## Annotate! In this case we are mapping the taxonomy against itself. 
## At this time only correlation (corr.map) and hierarchical (MapMyCells) are supported.
## Returns an S4 class with mapping results.
mapping.anno = taxonomy_mapping(AIT.anndata=AIT.anndata , ## Allen Institute Taxonomy loaded via `loadTaxonomy()`
                                query.data=Matrix::t(AIT.anndata$X[1:100,]), ## Gene x Cell Matrix
                                label.cols=c("Neighborhood_label", "Class_label", "Subclass_label", "Group_label"),
                                corr.map=TRUE,
                                hierarchical.map=TRUE,
                                tree.map=FALSE,
                                seurat.map=FALSE)

## Extract mapping results from S4 mappingClass
mapping.results = getMappingResults(mapping.anno, scores = FALSE)
```