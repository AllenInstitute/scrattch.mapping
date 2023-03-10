# Tutorial: Building a Shiny taxonomy 

In this tutorial we demonstrate how to setup a Shiny taxonomy using scrattch.mapping for viewing on MolGen Shiny and running mapping algorithms against. 

#### Required inputs:

* Count matrix (gene x cell)
* Annotation data.frame (cell x field), with sample identifiers as rownames
* Variable and/or marker genes (vector)
* UMAP coordinates (cell x 2)

#### Build taxonomy:

```R
## Load scrattch.mapping
library(scrattch.mapping)
library(umap)

## Load in example count data and annotations
library(tasic2016data)
taxonomy.counts = tasic_2016_counts
taxonomy.anno = tasic_2016_anno

## Ensure count matrix and annotations are in the same order.
taxonomy.anno = taxonomy.anno[match(colnames(taxonomy.counts), taxonomy.anno$sample_name),]

## Ensure cluster field exists, as required by Shiny and scrattch.mapping.
taxonomy.anno$cluster = taxonomy.anno$broad_type

## Compute top 1000 binary marker genes for clusters
binary.genes = top_binary_genes(taxonomy.counts, taxonomy.anno$cluster, 1000)

## Compute UMAP coordinates
umap.coords = umap(t(taxonomy.counts))$layout

## Set rownames to your annotation and UMAP data.frames (Required!)
rownames(taxonomy.anno) = taxonomy.anno$sample_name
rownames(umap.coords) = colnames(taxonomy.counts)

## Build Shiny taxonomy 
buildTaxonomy(counts = taxonomy.counts,
                meta.data = taxonomy.anno,
                feature.set = binary.genes,
                umap.coords = umap.coords,
                shinyFolder="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016",
                subsample=2000)

## Add markers to dendrogram
shinyFolder = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016"
addDendrogramMarkers(dend      = readRDS(file.path(shinyFolder,"dend.RData")), 
                    norm.data = file.path(shinyFolder,"data_t.feather"),
                    metadata  = file.path(shinyFolder,"anno.feather"),
                    shinyFolder = shinyFolder)
```
