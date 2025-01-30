# Tutorial: Standard taxonomy mapping using flat and hierarchical.

In this tutorial we demonstrate how to run multiple mapping algorithms using scrattch.mapping on the HMBA Macaque consensus basal ganglia taxonomy (applies to Human as well).

To run this code please use the docker provided at: `docker://njjai/scrattch_mapping:0.6.3`:

```R
library(scrattch.taxonomy)
library(scrattch.mapping)
library(reticulate)
cell_type_mapper <- import("cell_type_mapper") ## For now this has to be defined for mapmycells.hierarchical.map to work.

## ------------------------------------------------------
## Build Shiny taxonomy
AIT.anndata = loadTaxonomy(taxonomyDir="/PATH/TO/AIT",
                            anndata_file="HMBA_Macaque_BG_082024_AIT.h5ad")

## Annotate! In this case we are mapping the taxonomy against itself. 
## At this time four mapping algorithms are supported:
## -- (1) Simple correlation based mapping (corr.map)
## -- (2) Hierarchical mapping (mapmycells.hierarchical.map) - this is the method used in MapMyCells as default for all taxonomies
## -- (2) Flat mapping (mapmycells.flat.map) - this is the method used in MapMyCells for all taxonomies
## -- (3) Tree based mapping (tree.map) - This method requires a dendrogram and is the method used for several Patch-seq studies. NOT RECOMMENDED in most situations.
## -- (4) Seurat based mapping (seurat.map) - Mapping using TransferData from Seurat v4.4 with largely default parameters
## Returns an S4 class with mapping results.
hierarchy <- c("Neighborhood_label", "Class_label", "Subclass_label", "Group_label")
mapping.anno = taxonomy_mapping(AIT.anndata=AIT.anndata , ## Allen Institute Taxonomy loaded via `loadTaxonomy()`
                                query.data=Matrix::t(AIT.anndata$X[1:100,]), ## Gene x Cell Matrix
                                label.cols=hierarchy,
                                corr.map=TRUE,
                                mapmycells.hierarchical.map=TRUE,
                                mapmycells.flat.map=TRUE,
                                tree.map=TRUE,
                                seurat.map=TRUE,
                                mapmycells_params_list = list())

                                #Run list_hierarchical_params() to get a full list of parameters.
                                #Ex: mapmycells_params_list = list('type_assignment' = 
                                #                             list('bootstrap_iteration' = 100, 
                                #                             'bootstrap_factor' = 0.9))

## Extract mapping results and associated scores from S4 mappingClass
mapping.results = getMappingResults(mapping.anno, scores = TRUE)
```