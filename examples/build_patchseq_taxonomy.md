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
query.logCPM = logCPM(query.counts[,keep])
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

Let's start by defining the cell classes which are off target for tasic2016 patchseq mapping and QC. Typically this is defined at the class level and are the "Non-neuronal" cell classes/types.

```R
## Identify the offtarget cell types manually.
print(unique(AIT.anndata$obs$broad_type_label))

## Add in the off.target annotation.
levels(AIT.anndata$obs$broad_type_label) = c(levels(AIT.anndata$obs$broad_type_label), "Non-neuronal")
AIT.anndata$obs$broad_type_label[!is.element(AIT.anndata$obs$broad_type_label, c("GABA-ergic Neuron","Glutamatergic Neuron"))] = "Non-neuronal"
```

### Build the patchseq taxonomy:

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cell classes, typically Non-neuronals, from mapping.

```R
## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation.
                                    subclass.column = "cluster_label", ## Typically this is `subclass_label` but tasic2016 has no subclass annotation.
                                    class.column = "broad_type_label", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal"), ## The off-target class.column labels for patchseqQC.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    shinyFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016")
```
The `buildPatchseqTaxonomy` function return/created the following:

* An updated AIT.anndata object for patchseq mapping and QC steps.
* Saved the patchseq taxonomy anndata to `AI_taxonomy_patchseq.h5ad`.
* Created the required marker and expression files for patchseqQC.

### Map against the patchseq taxonomy:
```R
query.mapping = taxonomy_mapping(AIT.anndata= AIT.anndata,
                                  query.data = query.logCPM, 
                                  corr.map   = TRUE, # Flags for which mapping algorithms to run
                                  tree.map   = TRUE, 
                                  seurat.map = FALSE, 
                                  label.cols = c("cluster_label", "broad_type_label"), # Columns to map against from AIT.anndata$obs
                                  patchseq = TRUE)
```

### Determine patchseq contamination with PatchseqQC:
```R
query.mapping = applyPatchseqQC = function(AIT.anndata, ## A patchseq taxonomy object.
                                            query.counts, ## Counts are required here.
                                            query.mapping, ## Results of the previous mapping or AIT.anndata$obs, no mapping is required.
                                            verbose=FALSE)
```

### Setup the patchseq Shiny taxonomy files for MolGen Shiny:
```R
buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = '/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016/patchseq_mapping',
                      query.data     = query.counts, ## Counts are required here.
                      query.metadata = query.metadata,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = FALSE)  ## Set to FALSE if not needed or if buildPatchseqTaxonomy was not run.
dir('/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016/patchseq_mapping')
```