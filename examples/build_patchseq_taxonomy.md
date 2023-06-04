# Tutorial: Building a patchseq Shiny taxonomy 

In this tutorial we demonstrate how to setup a patchseq Shiny taxonomy using scrattch.mapping for viewing on MolGen Shiny and running mapping algorithms against. 

### Required inputs:

* Standard Shiny taxonomy setup following the "build_taxonomy" tutorial.
* Query patchseq count matrix and metadata.

### Load in Tasic 2016:
```R
library(scrattch.mapping)
library(tasic2016data)

## Load in the tasic2016 data, ignore this part. Tasic specific data wrangling.
query.anno = tasic_2016_anno
query.counts = tasic_2016_counts 
query.anno = query.anno[match(colnames(query.counts),query.anno$sample_name),]
rownames(query.anno) = query.anno$sample_name  
keep = query.anno$broad_type!="Unclassified"
query.counts = query.counts[,keep]
query.logCPM = logCPM(query.counts)
query.anno = query.anno[keep,]
```

### Load and setup the base Shiny Taxonomy:
```R
## Standard shiny taxonomy
taxonomy = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016"

## Load in the taxonomy
AIT.anndata = loadTaxonomy(taxonomy)
```

### Define off target cell types.

Let's start by defining the cell classes which are off target for tasic2016 patchseq mapping and patchseqQC. Typically this is defined at the class level and are the "Non-neuronal" cell classes/types.
```R
## Identify the offtarget cell types manually.
print(unique(AIT.anndata$obs$broad_type_label))

## Add in the off.target annotation.
AIT.anndata$obs$off_target = AIT.anndata$obs$broad_type_label

## First we need to add our off.target annotation to the factor levels
levels(AIT.anndata$obs$off_target) = c(levels(AIT.anndata$obs$off_target), "Non-neuronal")

## Now lets set all Non-neuronal cells to the "Non-neuronal" off target annotation.
AIT.anndata$obs$off_target[!is.element(AIT.anndata$obs$off_target, c("GABA-ergic Neuron","Glutamatergic Neuron", "Astrocyte"))] = "Non-neuronal"
```

### Build the patchseq taxonomy:

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. **You are creating a new version of the base taxonomy which can be reused by specifying the provided `mode.name` in `scrattch.mapping::mappingMode()` as dicusssed next.**

```R
## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "patchseq", ## Give a name to off.target filterd taxonomy
                                    subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation.
                                    subclass.column = "cluster_label", ## Typically this is `subclass_label` but tasic2016 has no subclass annotation.
                                    class.column = "off_target", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal"), ## The off-target class.column labels for patchseqQC.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016")
```
The `buildPatchseqTaxonomy` function return/created the following:

* An updated AIT.anndata object for patchseq mapping and QC steps.
* Created the required marker and expression files for patchseqQC and save under 'mode.name' sub-directory

### Set scrattch.mapping mode

Now we will set scrattch.mapping to use only cells not in the off.target.types, this will filter the taxonomy and adjust the dendrogram to remove any `cluster` in `off.target.types`.

```R
AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")
```

### Map against the patchseq taxonomy:
```R
query.mapping = taxonomy_mapping(AIT.anndata= AIT.anndata,
                                  query.data = query.logCPM, 
                                  corr.map   = TRUE, # Flags for which mapping algorithms to run
                                  tree.map   = TRUE, 
                                  seurat.map = FALSE, 
                                  label.cols = c("cluster_label", "broad_type_label")) # Columns to map against from AIT.anndata$obs
```

### Determine patchseq contamination with PatchseqQC:
```R
query.mapping = applyPatchseqQC(AIT.anndata, ## A patchseq taxonomy object.
                                query.counts, ## Counts are required here.
                                query.mapping, ## Results of the previous mapping or AIT.anndata$obs, no mapping is required.
                                verbose=FALSE)
```

### Setup the patchseq Shiny taxonomy files for MolGen Shiny:
```R
buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = '/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016/patchseq_mapping',
                      query.data     = query.counts, ## Counts are required here.
                      query.metadata = query.anno,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = FALSE)  ## Set to FALSE if not needed or if buildPatchseqTaxonomy was not run.
```
