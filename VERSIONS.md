## scrattch.mapping v1.2

Improvements for mapping and patch-seq for schema.

### Major changes
* `taxonomyMapping` now works correctly for all algorithms except mapmycells.flat.map and mapmycells.hierarchical.map
* Adding html documentation for functions
* Fix for any hierarchy names (not just cluster_id or cluster_label)
* Correct treatment of high variance / marker genes

### Minor changes
* Bug fixes
* Updated documentation
* Allow query.data to be inputted as either cellxgene or genexcell for `taxonomyMapping`
* Possibly fixes an error to tree mapping that was removing many genes from the mapping <- need to confirm!

--

## scrattch.mapping v0.7.1

Updates to improve correlation mapping

### Major changes
* Updated `taxonomyMapping` wrapper to allow flexibility, while keeping defaults and back-compatibility intact
* Updated `corrMap` to allow gene list or logical input, and to default to appropriate variable gene set
* Updated `seuratMap` the same way, and also to only pull mode-specific cells and clusters
* Updated all three of the above functions to properly consider the current mode  

### Minor changes
* Synchronized versions between scrattch.taxonomy, scrattch.mapping, an scrattch.patchseq
* Error correction
* Improved documentation on mapping methods

--

## scrattch.mapping v0.55.6

Updates to streamline and fix the examples

### Major changes

### Minor changes
* Updated examples
* Related minor bug fixes and function edits
* Updated help files

--

## scrattch.mapping v0.55.5

Major updates correcting issues with hierarchical and Seurat mapping

### Major changes
* Update `seuratMap.R` to correct Seurat mapping, by rolling back to Seurat v4.4 and largely rewriting the code
* Allowing hierarchical mapping to work on modes other than 'standard'
* Returning detailed information for hierarchical mapping (e.g., top 5 matches and bootstrap probabilities) in main results

### Minor changes
* Additional minor bug fixes
* Updated help files

--


## scrattch.mapping v0.52.2

Major change to how we return mapping results.

### Major changes
`taxonomy_mapping` now returns an S4 class called `mappingClass`. We've included a single method `getMappingResults()` to extract the results table from the S4 class.

### Minor changes
Many bug fixes and path handling

--

## scrattch.mapping v0.4

Remove code associated with taxonomy building to scrattch.taxonomy.

### Major changes

### Minor changes

--

## scrattch.mapping v0.21

This version of scrattch.mapping additionally simplifies the PatchSeq workflow and file structure by moving content from additional files to the h5ad and by allowing for a precalculated dendogram to be input to `buildTaxonomy()`.

### Major changes

* `buildTaxonomy()` allows for a precalculated dendrogram to be input (back-compatible, so by default the dendrogram is generated).

* Move the patchseqQC variables from file to Anndata object (only only remaining files outside h5ad are dend.RData, de.genes.rda, and membership_information_reference.rda, which we plan to move to the h5ad in future updates).

### Minor changes

* Create a new example for generating the [Hodge et al 2019](https://www.nature.com/articles/s41586-019-1506-7) MTG taxonomy and mapping patch-seq to it.

* Allow `loadTaxonomy()` to read either a folder location or an anndata file directly.

* Error correction in addDendrogramMarkers (and other functions).
  
* Auto-correct '\' --> '/' for different platforms in file paths to allow compatibility betwee Windows and Unix, through use of a `file.path()` wrapper (beta).

---

## scrattch.mapping v0.2

Updates to this version of scrattch.mapping are primarily associated with improving and simplifying the PatchSeq workflow by allowing for Taxonomy versions which filter out target cell types and can be re-used once built. Moving towards removal of label.cols from `taxonomy_mapping()`.

Updated examples for new workflows are provided on our github page.

### Major changes

* `buildPatchseqTaxonomy()` saves off-target cells into the taxonomy anndata for filtering to specific cell types. To be generalized past patchseq in the future.

* Enabled the dynamic modification of a `base` taxonomy into a filtered taxonomy `version` for mapping using Corr, Tree and Seurat.

* Taxonomy `versions` are stored in the .h5ad file and can be reused in future mappings by an user.

### Minor changes

* `buildTaxonomy()` now also creates the anndata object, this was done only in loadTaxonomy previously.

* New parameter to `buildTaxonomy()`, `taxonomyName` allows the user to encode a taxonomy ID that is stored in the taxonomy .h5ad which is currently fixed at `AI_taxonomy.h5ad`.

* `loadTaxonomy()` will use the anndata object if it exists to load the taxonomy. If not a new anndata object will be created from the feather files and saved while loading. 

* Exposed tree mapping parameters, can now set p, low.th and bootstrap now.

* Ensure that the dendrogram labels ("cluster_label") is used for all cluster level objects instead of "cluster_1, cluster_2 ...".

* Removed all dependence on "reference.rda", which is also no longer created.

* Move multiple functions to internal only, not exposed to the user.

### Parameter renaming

(new) | (old)
taxonomyDir | shinyFolder
taxonomyDir | refFolder

### Function renaming
(new) | (old)
buildPatchseqTaxonomy | writePatchseqQCmarkers

