# Tutorial: Hierarchical approximate nearest neighbors (HANN) mapping.

In this tutorial we demonstrate how to run HANN mapping using scrattch.mapping on the Tasic et al. 2016 study. For additional details on HANN mapping please refere to [confluence](http://confluence.corp.alleninstitute.org/pages/viewpage.action?pageId=117985236).

```R
## Load scrattch.mapping
library(scrattch.mapping)

## Load in example count data
library(tasic2016data)

## Compute log CPM
query.data = tasic_2016_counts
query.data = logCPM(query.data)

## Run CK's HANN mapping via scrattch.mapping
mapped = run_mapping_on_taxonomy(query.data,
                                 Taxonomy="AIT17.0_mouse",
                                 prefix="scrattch_mapping_test",               
                                 prebuild=FALSE,
                                 newbuild=TRUE,
                                 mapping.method='hierarchy',
                                 mc.cores=10,
                                 iter=10,
                                 blocksize=5000)

## Extract results
results = mapped$cl.df[match(mapped$best.map.df$best.cl, mapped$cl.df$cl),]
```