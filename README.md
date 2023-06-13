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
`singularity exec --cleanenv docker://bicore/scrattch_mapping:latest Rscript YOUR_CODE.R`

##### Interactive
`singularity shell --cleanenv docker://bicore/scrattch_mapping:latest`


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

1. [**Run heirarchical approximat nearest neighbords (HANN) taxonomy mapping**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/HANN_mapping.md) This example shows how to use Zizhen/CK's HANN taxonomy mapping.

2. [**Run Flat, Tree, and Seurat taxonomy mapping**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/mapping.md) This examples shows how to use scrattch.mapping for standard taxonomy mapping.

3. [**Build a Shiny taxonomy**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/build_taxonomy.md) This examples provides the basics for creating a new Shiny taxonomy compatible with MolGen shiny and scrattch.mapping.

4. [**Build and map against a small mouse PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/build_patchseq_taxonomy.md) This examples provides the basics for updating a Shiny taxonomy to be compatible with patchseq style mapping and visualization on MolGen Shiny.

5. [**Build and map against a human MTG PatchSeq taxonomy**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/build_MTG_patchseq_taxonomy.md) This examples provides a simplified example for creating a standard taxonomy and updating it to be compatible with patchseq style mapping and visualization on MolGen Shiny.  Data is from Hodge et al for human MTG and patch-seq examples are from Berg et al (layer 2-3, excitatory neurons). 

## Library vignettes

*Note: links below will not work while repo is private.*

1. [**Build reference directory for mapping.**](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/build_reference_taxonomy.html)  This vignette provides an example of how to convert a *completed* single cell RNA-seq analysis (e.g., a counts matrix + cell type assignments) into a standard reference taxonomy. Resulting taxonomy files are used as input for various mapping techniques in this package, and are also compatible with tools for visualiation of taxonomies at the Allen Institute and Patch-seq QC and visualization. **This process must be run first.**  

2. [**Map patch-seq data and output directory.**](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/complete_patchseq_analysis.html)  This vignette goes through how to map a small data set against a reference taxonomy. Here we use a subset of tasic2016data as an example but the intention is for mapping of patch-seq data.  

3. [**Comparing multiple mapping methods.**](http://htmlpreview.github.io/?https://github.com/AllenInstitute/mfishtools/blob/master/vignettes/comparison_of_mapping_methods.html)  This vignette goes through how to run multiple mapping algorithms on a data set for label transfer, and then compare and contrast the results of these mapping algorithms. It also includes a cursory look into how to use ground truth (if available) to aid in selection of optimal mapping results. Here we use the entirely of tasic2016data as an example so that we can have “ground truth” clustering results to compare against, but in practice this information is typically not available.  

## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-mapping/issues).

## TODO

- [ ] Generalize HANN mapping for external use.
- [ ] Standardize output from multiple modes of mapping (RNA, Patch-Seq, spatial?)
- [ ] Allow other values than 100 for bootstrapping
- [ ] Minor vignette updates (e.g, buildReference --> buildTaxonomy)
- [x] Allow hGenes to impact corr and Seurat mapping, but have tree and HANN ignore it (it will break the node by node marker approach)
- [ ] Pass a dendrogram to buildTaxonomy if one exists already (I don't remember if this is allowed, but if not, it may be more complicated than it sounds, which is probably why I haven't done it).
- [x] load the AIT.anndata directly from the .h5ad if it exists
- [x] Fixes to addDendrogramMarker based on ET feedback.


## Done

- [x] Package up the stand-alone script.
- [x] Add stand-alone scripts for creating a reference sn/scRNA-seq taxonomy to package
- [x] Accept various inputs (Seurat, Anndata, count matrix + metadata table, reading from reference taxonomy folder)?
- [x] Write in a module for making the taxonomy compatible with tree mapping (adding markers to dendrogram). Perhaps this should be done after the taxonomy is made for generalization. (JM: I was planning to take this one)
- [x] Add stand-alone scripts for patch-seq QC to package
- [x] Output feather with gene counts along with required Shiny files for use with scVI/scANVI etc. I’ve had to manually make a few of these for the older BG taxonomies.
- [x] This might be the same thing, but also output files for viewing data in existing patch-seq shiny folder format.
- [x] Implement new tree mapping (Code from CK)
- [x] Dockerize.
