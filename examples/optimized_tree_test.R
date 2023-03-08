## Test that should be run after "build_test.R"
## export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
## sudo docker build - -t njjtemp/scrattch-mapping:2.6 -f /home/nelson.johansen/scripts/docker < Dockerfile
## sudo docker push njjtemp/scrattch-mapping:2.6
## singularity exec --cleanenv docker://njjtemp/scrattch-mapping:2.6 R CMD build --no-build-vignettes /home/nelson.johansen/scripts/R/github/scrattch-mapping/
## singularity exec --cleanenv docker://njjtemp/scrattch-mapping:2.6 Rscript /home/nelson.johansen/scripts/R/github/scrattch-mapping/examples/optimized_tree_test.R


library(scrattch.mapping)

##
# setwd("/allen/programs/celltypes/workgroups/hct/NelsonJ/glial/")
# adata = read_h5ad("smarter_annSeurat.h5ad");
# mat = t(as.matrix(adata$X)); meta.data = adata$obs;
# mat = logCPM(mat)

# ##
# # source("/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/Mapping_On_Taxonomy/hmapping/hmapping_util.R")
# mapped = run_mapping_on_taxonomy(
#                            mat,          # query data
#                            Taxonomy="AIT17.2_mouse",    #taxonomy # AIT9.0_mouse : Fore Integrated
#                            prefix="SS",               # prefix for your data platform (SS, 10X_cells_v2, 10X_cells_v3, 10X_nuclei_v3,MFISH)
#                            prebuild=FALSE,            # if the taxonomy is available for your platform, use it.
#                            newbuild=FALSE,
#                            mapping.method='hierarchy',# 'flat','hierarchy'
#                            mc.cores=10,               # lower this to 5 if the run fails due to the memory
#                            iter=10,
#                            blocksize=5000)
# results = mapped$cl.df[match(mapped$best.map.df$best.cl, mapped$cl.df$cl),]
# save(results, mapped, file="/allen/programs/celltypes/workgroups/hct/NelsonJ/glial/AIT172_mapping.rda")

# ## Load in the data to be annotated
# load("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/pipeline/1001_C06_MTX-2046/1001_C06.RData")
# mat = as.matrix(GetAssayData(rnaseq.data, "data"))

# ##
# rownames(mat) = tolower(rownames(mat))
# rownames(mat) = str_to_title(rownames(mat))

# ##
# mapped = run_mapping_on_taxonomy(
#                            mat,          # query data
#                            Taxonomy="AIT17.2_mouse",    #taxonomy # AIT9.0_mouse : Fore Integrated
#                            prefix="SS",   # prefix for your data platform (SS, 10X_cells_v2, 10X_cells_v3, 10X_nuclei_v3, MFISH)
#                            prebuild=FALSE,            # if the taxonomy is available for your platform, use it.
#                            newbuild=FALSE,
#                            mapping.method='hierarchy',# 'flat','hierarchy'
#                            mc.cores=10,               # lower this to 5 if the run fails due to the memory
#                            iter=10,
#                            blocksize=5000)
# results = mapped$cl.df[match(mapped$best.map.df$best.cl, mapped$cl.df$cl),]
# print(results)

# ##
# save(mapping.anno, file="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/NHP/pipeline/1001_C06_MTX-2046/mapping_anno.RData")