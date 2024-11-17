# scrattch.mapping

Generalized mapping scripts for RNA-seq, Patch-seq or any gene expression data. Assumes that you have built an Allen Institute Taxonomy (AIT) object with [scrattch.taxonomy](https://github.com/AllenInstitute/scrattch.taxonomy).

## Documentation

You can find a detail description of all scrattch.mapping functions here: ![Documentation](https://github.com/AllenInstitute/scrattch-mapping/blob/main/scrattch.mapping_0.1.pdf)

Update notes are here: ![Versions](https://github.com/AllenInstitute/scrattch-mapping/blob/dev_njj/VERSIONS.md)

## Installation

### Using docker (recommended)
We have setup a docker environment for scrattch.taxonomy, scrattch.mapping, and scrattch.patchseq that contains all the required dependencies and the current version of all scrattch packages. **See [the readme](https://github.com/AllenInstitute/scrattch/blob/master/README.md#using-docker) for [the parent scrattch package](https://github.com/AllenInstitute/scrattch) for the most up-to-date docker information.**

### Directly from GitHub (strongly discouraged)

While we advise using the provided docker, you can also install scrattch.mapping directly from GitHub as follows:

```
devtools::install_github("AllenInstitute/scrattch.mapping")
```

This strategy **might not work** due to complicated dependencies. Also note that `doMC` may need to be installed manually from [HERE](https://r-forge.r-project.org/R/?group_id=947) if you use Windows. Vignettes are provided below.

## Usage examples

1. [**Run Flat, Tree, and Seurat taxonomy mapping**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/mapping.md) This examples shows how to use scrattch.mapping for standard taxonomy mapping.
2. [**Mapping to HMBA Basal Ganglia AIT**](https://github.com/AllenInstitute/scrattch-mapping/blob/main/examples/mapping_BasalGanglia.md) This tutorial shows how to map against the HMBA Human and Macaque Basal Ganglia consensus taxonomies.

## Reporting issues

If you run into any issues, please let Nelson and Jeremy know or [**create a new issue in the 'Issues' tab above**](https://github.com/AllenInstitute/scrattch-mapping/issues).

