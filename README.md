# scrattch.mapping

Generalized mapping scripts for RNA-seq, Patch-seq or any gene expression data. Assumes that you have built an Allen Institute Taxonomy (AIT) object with [scrattch.taxonomy](https://github.com/AllenInstitute/scrattch.taxonomy).

## Documentation

You can find a detail description of all scrattch.mapping functions here: ![Documentation](https://github.com/AllenInstitute/scrattch-mapping/blob/main/scrattch.mapping_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch-mapping/blob/dev_njj/VERSIONS.md)

## Docker

We have setup a docker environemnt for scrattch.taxonomy and scrattch.mapping that contains all the required dependencies and the current version of both scrattch.taxonomy and scrattch.mapping. This docker is accessible through docker hub via: `njjai/scrattch_mapping:0.6.3`.

#### HPC usage:

##### Non-interactive
`singularity exec --cleanenv docker://njjai/scrattch_mapping:0.6.3 Rscript YOUR_CODE.R`

##### Interactive
`singularity shell --cleanenv docker://njjai/scrattch_mapping:0.6.3`


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
2. [**Mapping to HMBA Basal Ganglia AIT**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/mapping_BasalGanglia.md) This tutorial shows how to map against the HMBA Human and Macaque Basal Ganglia consensus taxonomies.

## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-mapping/issues).

