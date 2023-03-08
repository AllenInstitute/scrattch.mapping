# scrattch-mapping

Generalized mapping scripts for RNA-seq and Patch-seq data

## Workflow

![](https://github.com/AllenInstitute/scrattch-mapping/blob/main/schematic.png)

## Installation

*Note: slight edits to installation will be needed while repo is private.*

```
# Quickly, but without the vignettes:
devtools::install_github("AllenInstitute/scrattch-mapping")

# More slowly, but with the vignettes:
devtools::install_github("AllenInstitute/scrattch-mapping",build_vignettes=TRUE, force=TRUE)
# Note that this strategy might not work outside the docker. Vignettes are also available below.
```

## Dependencies
### Remotes:
    AllenInstitute/scrattch.io,
    AllenInstitute/scrattch.hicat,
    AllenInstitute/scrattch.vis,
    AllenInstitute/mfishtools,
    AllenInstitute/patchseqtools,
    PavlidisLab/patchSeqQC,

## Library vignettes

*Note: links below will not work while repo is private.*

1. [**Build reference directory for mapping.**](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/build_reference_taxonomy.html)  This vignette provides an example of how to convert a *completed* single cell RNA-seq analysis (e.g., a counts matrix + cell type assignments) into a standard reference taxonomy. Resulting taxonomy files are used as input for various mapping techniques in this package, and are also compatible with tools for visualiation of taxonomies at the Allen Institute and Patch-seq QC and visualization. **This process must be run first.**  
2. [**Map patch-seq data and output directory.**](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/complete_patchseq_analysis.html)  This vignette goes through how to map a small data set against a reference taxonomy. Here we use a subset of tasic2016data as an example but the intention is for mapping of patch-seq data.  
3. [**Comparing multiple mapping methods.**](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/comparison_of_mapping_methods.html)  This vignette goes through how to run multiple mapping algorithms on a data set for label transfer, and then compare and contrast the results of these mapping algorithms. It also includes a cursory look into how to use ground truth (if available) to aid in selection of optimal mapping results. Here we use the entirely of tasic2016data as an example so that we can have “ground truth” clustering results to compare against, but in practice this information is typically not available.  

## TODO
 
- [x] Package up the stand-alone script.
- [x] Add stand-alone scripts for creating a reference sn/scRNA-seq taxonomy to package
- [x] Accept various inputs (Seurat, Anndata, count matrix + metadata table, reading from reference taxonomy folder)?
- [x] Write in a module for making the taxonomy compatible with tree mapping (adding markers to dendrogram). Perhaps this should be done after the taxonomy is made for generalization. (JM: I was planning to take this one)
- [x] Add stand-alone scripts for patch-seq QC to package
- [x] Output feather with gene counts along with required Shiny files for use with scVI/scANVI etc. I’ve had to manually make a few of these for the older BG taxonomies.
- [x] This might be the same thing, but also output files for viewing data in existing patch-seq shiny folder format. 
- [ ] Standardize output from multiple modes of mapping (RNA, Patch-Seq, spatial?)
- [ ] Implement new tree mapping (Code from CK)
- [ ] Dockerize.
