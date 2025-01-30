# Tutorial: Standard taxonomy mapping using flat, tree, hierarchical, and Seurat.

In this tutorial we demonstrate how to run standard mapping algorithms using scrattch.mapping on query data from [Gouwens, Sorensen, et al 2020 study](https://doi.org/10.1016/j.cell.2020.09.057) using the [Tasic et al. 2016 study](https://www.nature.com/articles/nn.4216) reference taxonomy created in the [scrattch.taxonomy example](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/examples/build_taxonomy.md). 

```R
## Load scrattch.mapping
library(reticulate)
library(scrattch.mapping)
library(scrattch.taxonomy)
cell_type_mapper <- import("cell_type_mapper")

## Download some example patch-seq data to map
## -- These data are from Gouwens et al 2020 and would be replaced by your query data
patchFolder  <- "https://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/"
counts_url   <- "20200513_Mouse_PatchSeq_Release_count.csv.tar"
download.file(paste0(patchFolder,counts_url),counts_url)
untar(counts_url)
query.counts <- as.matrix(data.table::fread("20200513_Mouse_PatchSeq_Release_count/20200513_Mouse_PatchSeq_Release_count.csv"),rownames=1)

## Compute log CPM
query.data = logCPM(query.counts)

## Specify which reference taxonomy to map against.
## -- Replace folder and file name with correct locations
taxonomyDir = getwd() 
AIT.anndata = loadTaxonomy(taxonomyDir = taxonomyDir, anndata_file="Tasic2016.h5ad")

## Map! Returns an S4 class with mapping results.
mapping.anno = taxonomy_mapping(AIT.anndata=AIT.anndata,
                                query.data=query.data,
                                corr.map=TRUE,
                                tree.map=TRUE,
                                mapmycells.hierarchical.map=TRUE,
                                mapmycells.flat.map=TRUE,
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