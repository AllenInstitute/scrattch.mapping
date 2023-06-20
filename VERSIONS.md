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

