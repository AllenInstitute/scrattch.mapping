# Tutorial: Standard taxonomy mapping to public AIT file.

In this tutorial we demonstrate how to run standard mapping algorithms using scrattch.mapping on query data from [Gouwens, Sorensen, et al 2020 study](https://doi.org/10.1016/j.cell.2020.09.057) using the [Tasic et al. 2018 study](https://www.nature.com/articles/s41586-018-0654-5) reference taxonomy available [in the table of available AIT taxonomies](https://github.com/AllenInstitute/AllenInstituteTaxonomy/blob/main/taxonomies.md). 

*Note: This example was built using `docker://alleninst/scrattch:1.2`.*

```R
## Load scrattch.mapping
library(reticulate)
library(scrattch.mapping)
library(scrattch.taxonomy)
cell_type_mapper <- import("cell_type_mapper")

## Set the working directory
taxonomyDir = getwd() 

## Download, load, and normalize some example patch-seq data to map
## -- These data are from Gouwens et al 2020 and would be replaced by your query data
patchFolder  <- "https://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/"
counts_url   <- "20200513_Mouse_PatchSeq_Release_count.csv.tar"
download.file(paste0(patchFolder,counts_url),counts_url)
untar(counts_url)
query.counts <- as.matrix(data.table::fread("20200513_Mouse_PatchSeq_Release_count/20200513_Mouse_PatchSeq_Release_count.csv"),rownames=1)
query.data   <- logCPM(query.counts)

## Specify which reference taxonomy to map against.
AIT.anndata = loadTaxonomy("https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Mouse_VISp_ALM_SMART_seq_04042025.h5ad")

## Set up the taxonomy for tree mapping (this is not done in advance for taxonomies on AWS)
# This step is slow! Comment and set "tree.map=FALSE" below if tree mapping is not desired
AIT.anndata = addDendrogramMarkers(AIT.anndata)

## Map! Returns an S4 class with mapping results.
mapping.anno = taxonomy_mapping(AIT.anndata=AIT.anndata,
                                query.data=query.data,
                                corr.map=TRUE,
                                tree.map=TRUE,
                                mapmycells.hierarchical.map=FALSE,
                                mapmycells.flat.map=FALSE,
                                seurat.map=TRUE,
                                mapmycells_params_list = list())
                                
                                #Run list_hierarchical_params() to get a full list of parameters.
                                #Ex: mapmycells_params_list = list('type_assignment' = 
                                #                             list('bootstrap_iteration' = 100, 
                                #                             'bootstrap_factor' = 0.9))

## Extract mapping results and associated scores from S4 mappingClass
mapping.results = getMappingResults(mapping.anno, scores = TRUE)

## Extract tree mapping bootstrapping table 
tree.bootstraps = mapping.anno@detailed_results[["tree"]]

## Extract mapmycells hierachical mapping bootstrapping probabilities for the top five cell mapped hits 
hierarchical.bootstraps = mapping.anno@detailed_results[["hierarchical"]]

## Extract mapmycells flat mapping bootstrapping probabilities for the top five cell mapped hits 
flat.bootstraps = mapping.anno@detailed_results[["flat"]]

## Save
save(mapping.anno, file="mapping_results.rda")
write.csv(mapping.results,"mapping_results.csv")
```