# Tutorial: Building a patchseq Shiny taxonomy - Human MTG

In this example we demonstrate how to setup a patchseq Shiny taxonomy using scrattch.mapping for viewing on MolGen Shiny and running mapping algorithms against. This tutorial parallels the other tutorial for building a patchseq Shiny taxonomy, but using the Hodge et al taxonomy as reference and query patch-seq data from Berg et al 2020.  

### Required inputs:

* Standard Shiny taxonomy setup following the "build_taxonomy" tutorial (repeated below)
* Query patchseq count matrix and metadata (example provided below)

### Read in the REFERENCE MTG information (Hodge et al 2019)
```R
## Load the library
library(scrattch.mapping)
library(data.table) # for using fread below

## Load the complete dendrogram from this paper (saved in scrattch-mapping)
data(dend_Hodge2019) 

## Download data and metadata from the website (this is slow)
## NOT AVAILABLE YET. Instead copy Hodge2019_metadata.csv and Hodge2019_counts.csv.gz from "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.2/" to your working directory.

## Read data and metadata into R
taxonomy.metadata <- read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.2/Hodge2019_metadata.csv",row.names=1)
taxonomy.counts   <- as.matrix(fread("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.2/Hodge2019_counts.csv.gz"),rownames=1) ## This requires R.utils
colnames(taxonomy.counts) <- rownames(taxonomy.metadata) # To correct "-" to "." conversion introduced at some point. 
```

### Create the base Shiny Taxonomy for the ENTIRE Hodge et al 2019 data set
```R
## This is where our taxonomy will be created
taxonomy = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.2_test/"

## Compute top 1000 binary marker genes for clusters
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.metadata$cluster_label, 1000)

## Compute UMAP coordinates
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout

## Set rownames to your annotation and UMAP data.frames (Required!)
rownames(umap.coords) = colnames(taxonomy.counts)

## Build Shiny taxonomy 
AIT.anndata = buildTaxonomy(counts = taxonomy.counts,
              meta.data   = taxonomy.metadata,
              dend        = dend_Hodge2019,  # Can be omitted and buildTaxonomy will generate a dendrogram
              feature.set = binary.genes,
              umap.coords = umap.coords,
              taxonomyName= "MTG_Hodge2019", ## NEW!
              taxonomyDir = taxonomy,
              subsample   = 100)

## Load the taxonomy (from h5ad file name)
AIT.anndata = loadTaxonomy(taxonomy, "AI_taxonomy.h5ad")
```

### Build the patchseq taxonomy:

Now let's create a version of the taxonomy which is compatible with patchseqQC and can be filtered to remove off target cells from mapping. **You are creating a new version of the base taxonomy which can be reused by specifying the provided `mode.name` in `scrattch.mapping::mappingMode()` as dicusssed next.**

```R
## Setup the taxonomy for patchseqQC to infer off.target contamination
AIT.anndata = buildPatchseqTaxonomy(AIT.anndata,
                                    mode.name = "patchseq", ## Give a name to off.target filterd taxonomy
                                    subsample = 100, ## Subsampling is only for PatchseqQC contamination calculation.
                                    subclass.column = "subclass_label", 
                                    class.column = "class_label", ## The column by which off-target types are determined.
                                    off.target.types = c("Non-neuronal"), ## The off-target class.column labels for patchseqQC.
                                    num.markers = 50, ## Number of markers for each annotation in `class_label`
                                    taxonomyDir = taxonomy) # This will create a subfolder in the reference taxonomy directory
```
The `buildPatchseqTaxonomy` function return/created the following:

* An updated AIT.anndata object for patchseq mapping and QC steps.
* Created the required marker and expression variables for 'QC_markers' and save under 'mode.name' in the uns
* Created the required cell to cluster 'membership' variables and save under 'mode.name' in the uns
* An updated dendrogram and a saved in the 'mode.name' subdirectory

**At this point the reference taxonomy is created and ready for patch-seq mapping.**


### Read in and map patch-seq data

The rest of this example demonstrates how to read in patch-seq data and map it to the neuronal cell types from the Hodge et al 2019 taxonomy. 

### Now read in and process QUERY patch-seq data
```R
## Download data and metadata from GitHub repo for Berg et al 2022
download.file("https://github.com/AllenInstitute/patchseq_human_L23/raw/master/data/input_patchseq_data_sets.RData", "patchseq.RData", mode="wb")

## Load the data
load("patchseq.RData")

## Rename the query data and metadata for convenience
query.anno = annoPatch  # Some cell annotations for all cells from the paper 
query.logCPM = datPatch # logCPM values for all cells from the paper
```


### Set scrattch.mapping mode

Now we will set scrattch.mapping to use only cells not in the off.target.types, this will filter the taxonomy and adjust the dendrogram to remove any `cluster` in `off.target.types`.
```R
AIT.anndata = mappingMode(AIT.anndata, mode="patchseq")
```

### Map against the patchseq taxonomy

This is the step that performs the mapping.  Since we have very different query and reference data sets, we are omitting Seurat mapping in this example. 

```R
query.mapping = taxonomy_mapping(AIT.anndata= AIT.anndata,
                                 query.data = query.logCPM, 
                                 corr.map   = TRUE, # Flags for which mapping algorithms to run
                                 tree.map   = TRUE, 
                                 seurat.map = FALSE, 
                                 label.cols = c("cluster_label", "subclass_label" ,"class_label")) # Columns to map against from AIT.anndata$obs
```

### Setup the patchseq Shiny taxonomy files for human MTG:

This step outputs the files necessary for visualization of Patch-seq data with molgen-shiny tools.  *Note that this code block also generates addional QC metrics including NMS with PatchseqQC.  To perform patch-seq QC without building this directory, use the `applyPatchseqQC` function.*

```R
buildMappingDirectory(AIT.anndata    = AIT.anndata, 
                      mappingFolder  = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT15.2/TEST",
                      query.data     = query.logCPM, ## Counts or CPM are required here, but function can convert from log to linear values
                      query.metadata = query.anno,
                      query.mapping  = query.mapping,
                      doPatchseqQC   = TRUE)  ## Set to FALSE if not needed or if buildPatchseqTaxonomy was not run.
```