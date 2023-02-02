# scrattch-mapping

Genearlized mapping scripts for RNA-seq and Patch-seq data

## Workflow

![](https://github.com/AllenInstitute/scrattch-mapping/blob/main/schematic.jpg)

## Dependencies
### Remotes:
    AllenInstitute/scrattch.io,
    AllenInstitute/scrattch.hicat,
    AllenInstitute/scrattch.vis,
    AllenInstitute/mfishtools,
    AllenInstitute/patchseqtools,
    PavlidisLab/patchSeqQC,

## TODO
 
- [x] Package up the stand-alone script.
- [x] Add stand-alone scripts for creating a reference sn/scRNA-seq taxonomy to package
- [x] Accept various inputs (Seurat, Anndata, count matrix + metadata table, reading from reference taxonomy folder)?
- [x] Write in a module for making the taxonomy compatible with tree mapping (adding markers to dendrogram). Perhaps this should be done after the taxonomy is made for generalization. (JM: I was planning to take this one)
- [ ]Add stand-alone scripts for patch-seq QC to package
- [x] Output feather with gene counts along with required Shiny files for use with scVI/scANVI etc. I’ve had to manually make a few of these for the older BG taxonomies.
- [ ]This might be the same thing, but also output files for viewing data in existing patch-seq shiny folder format
- [ ]Dockerize.
