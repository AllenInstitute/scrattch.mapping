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
library(scrattch.mapping, lib.loc="/home/nelson.johansen/R/x86_64-pc-linux-gnu-library/4.2")
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
pcs  <- prcomp(logCPM(taxonomy.counts)[binary.genes,], scale = TRUE)$rotation
umap.coords = umap(pcs[,1:30])$layout

## Set rownames to your annotation and UMAP data.frames (Required!)
rownames(taxonomy.anno) = taxonomy.anno$sample_name
rownames(umap.coords) = colnames(taxonomy.counts)

## This is where our taxonomy will be created
taxonomyDir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/tasic_2016"

## Build Shiny taxonomy 
AIT.anndata = buildTaxonomy(counts = taxonomy.counts,
                meta.data = taxonomy.anno,
                feature.set = binary.genes,
                umap.coords = umap.coords,
                taxonomyName = "Tasic2016", ## NEW!
                taxonomyDir = taxonomyDir,
                subsample=2000)

## Add markers to dendrogram
addDendrogramMarkers(AIT.anndata = AIT.anndata)
```
