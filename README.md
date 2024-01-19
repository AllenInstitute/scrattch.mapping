# scrattch.mapping

Generalized mapping and taxonomy building scripts for RNA-seq, Patch-seq or any gene expression data.

## Workflow

![](https://github.com/AllenInstitute/scrattch-mapping/blob/main/schematic.png)

## Documentation

You can find a detail description of all scrattch.mapping functions here: ![Documentation](https://github.com/AllenInstitute/scrattch-mapping/blob/main/scrattch.mapping_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch-mapping/blob/dev_njj/VERSIONS.md)

## Docker

We have setup a docker environemnt for scattch.mapping that contains all the required dependencies and the current version of scrattch.mapping. This docker is accessible through docker hub via: `bicore/scrattch_mapping:latest`.

#### HPC usage:

##### Non-interactive
`singularity exec --cleanenv docker://njjai/scrattch_mapping:0.41 Rscript YOUR_CODE.R`

##### Interactive
`singularity shell --cleanenv docker://njjai/scrattch_mapping:0.41`


## Installation

While we advice using the provided docker, you can also install scrattch.mapping directly from github as follows:

*Note: slight edits to installation will be needed while repo is private.  Also note that `doMC` may need to be installed manually from the download at https://r-forge.r-project.org/R/?group_id=947 if you use Windows.*

```
# Quickly, but without the vignettes:
devtools::install_github("AllenInstitute/scrattch-mapping")

# More slowly, but with the vignettes:
devtools::install_github("AllenInstitute/scrattch-mapping",build_vignettes=TRUE, force=TRUE)
```

Note that this strategy might not work outside the docker due to complicated dependencies. Vignettes are provided below.

## Usage examples

1. [**Run Flat, Tree, and Seurat taxonomy mapping**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/mapping.md) This examples shows how to use scrattch.mapping for standard taxonomy mapping.

## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-mapping/issues).

## TODO

- [ ] Generalize HANN mapping for external use.
- [ ] Standardize output from multiple modes of mapping (RNA, Patch-Seq, spatial?)

## Done

- [x] Allow other values than 100 for bootstrapping
- [x] Minor vignette updates (e.g, buildReference --> buildTaxonomy)
- [x] Allow hGenes to impact corr and Seurat mapping, but have tree and HANN ignore it (it will break the node by node marker approach)
- [x] Pass a dendrogram to buildTaxonomy if one exists already (I don't remember if this is allowed, but if not, it may be more complicated than it sounds, which is probably why I haven't done it).
- [x] load the AIT.anndata directly from the .h5ad if it exists
- [x] Fixes to addDendrogramMarker based on ET feedback.
- [x] Package up the stand-alone script.
- [x] Add stand-alone scripts for creating a reference sn/scRNA-seq taxonomy to package
- [x] Accept various inputs (Seurat, Anndata, count matrix + metadata table, reading from reference taxonomy folder)?
- [x] Write in a module for making the taxonomy compatible with tree mapping (adding markers to dendrogram). Perhaps this should be done after the taxonomy is made for generalization. (JM: I was planning to take this one)
- [x] Add stand-alone scripts for patch-seq QC to package
- [x] Output feather with gene counts along with required Shiny files for use with scVI/scANVI etc. I’ve had to manually make a few of these for the older BG taxonomies.
- [x] This might be the same thing, but also output files for viewing data in existing patch-seq shiny folder format.
- [x] Implement new tree mapping (Code from CK)
- [x] Dockerize.
