FROM rocker/tidyverse:4.2

MAINTAINER Nelson Johansen nelson.johansen@alleninstitute.org

RUN export GITHUB_PAT=1000

## Would have liked to do these next 2 RUN under FROM:python3.8 but artifact passing wasn't immediatly clear. build-essential libxml2-dev python-is-python3
RUN apt-get update
RUN apt-get install -y wget python3-dev python3-pip
RUN pip3 install anndata==0.8.0 numpy

RUN R -e 'install.packages("reticulate")'
RUN R -e 'install.packages("anndata", update=TRUE)'

RUN R -e 'install.packages("BiocManager", update=FALSE)' 
RUN R -e 'BiocManager::install(c( "AnnotationDbi", "data.table", "GO.db", \
                                  "impute", "limma", "preprocessCore", "xml2", "rols"), dependenceis=NA, update=TRUE)' 
RUN R -e 'BiocManager::install(c( "munsell", "rhdf5", "dplyr", \
                                  "optparse", "foreach", "doParallel", "futile.logger", \
                                  "ggplot2", "WGCNA"), dependenceis=NA, update=TRUE)' 
RUN R -e 'BiocManager::install(c( "randomForest", "LaplacesDemon", "reshape2", \
                                  "feather", "future", "tibble", "dendextend", \
                                  "Matrix" ), dependenceis=NA, update=TRUE)' 
RUN R -e 'BiocManager::install(c( "mgcv", "edgeR", "caret", \
                                  "ggbeeswarm", "pvclust", \
                                  "cowplot" ), dependenceis=NA, update=TRUE)'
RUN R -e 'BiocManager::install(c("bigstatsr", "umap"), dependenceis=NA, update=TRUE)' 
RUN R -e 'BiocManager::install(c("beachmat", "BiocNeighbors"), dependenceis=NA, update=TRUE)' 

RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.5.tar.gz", repos=NULL, type="source")'
RUN R -e 'install.packages("https://cloud.r-project.org/src/contrib/profmem_0.6.0.tar.gz", repos=NULL, type="source")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.8-0.tar.gz", repos=NULL, type="source")'

RUN R -e 'install.packages("reticulate")'
RUN R -e 'install.packages("arrow")'
RUN R -e 'install.packages("anndata", update=TRUE)'
RUN R -e 'install.packages("jsonlite", update=TRUE)'

## 
COPY matrixStats_1.3.0.tar.gz ./matrixStats_1.3.0.tar.gz
RUN R -e 'install.packages("matrixStats_1.3.0.tar.gz", repos=NULL, type="source")'

## Seurat setup
RUN apt-get update && \
    apt-get install -y libgeos-dev libglpk-dev
RUN R -e 'install.packages("sp")'
RUN R -e 'install.packages("spam")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz", repos=NULL, type="source",)'
RUN R -e 'install.packages("SeuratObject")'
RUN R -e 'install.packages("Seurat")'

## Remote installs
RUN R -e 'install.packages("remotes", update=TRUE)'
RUN R -e 'remotes::install_github("krlmlr/bindrcpp")'
RUN R -e 'remotes::install_github("igraph/rigraph")'
RUN R -e 'remotes::install_github("i-cyto/Rphenograph")' 
RUN R -e 'remotes::install_github("PavlidisLab/patchSeqQC")'
RUN R -e 'remotes::install_github("cysouw/qlcMatrix")'

## Allen Institute R installs
RUN R -e 'remotes::install_github("AllenInstitute/CCN")'
RUN R -e 'remotes::install_github("AllenInstitute/patchseqtools")' 
RUN R -e 'remotes::install_github("AllenInstitute/tasic2016data")'
RUN R -e 'remotes::install_github("AllenInstitute/hodge2019data")'
RUN R -e 'remotes::install_github("AllenInstitute/scrattch.io", ref="dev")' 
RUN R -e 'remotes::install_github("AllenInstitute/scrattch.vis")' 
RUN R -e 'remotes::install_github("AllenInstitute/scrattch.hicat")' 
RUN R -e 'remotes::install_github("AllenInstitute/scrattch.bigcat")'
RUN R -e 'remotes::install_github("AllenInstitute/mfishtools")'

## scrattch-taxonomy install from local source
COPY scrattch.taxonomy_0.5.11.tar.gz ./scrattch.taxonomy_0.5.11.tar.gz
RUN R -e 'install.packages("scrattch.taxonomy_0.5.11.tar.gz", repos=NULL, type="source")'

## scrattch-mapping install from local source
COPY scrattch.mapping_0.6.1.tar.gz ./scrattch.mapping_0.6.1.tar.gz
RUN R -e 'install.packages("scrattch.mapping_0.6.1.tar.gz", repos=NULL, type="source")'

# cell_type_mapper install from GitHub
RUN git clone -b update/uns/to/precomp/stats/params --single-branch https://github.com/AllenInstitute/cell_type_mapper.git
RUN pip install -r ./cell_type_mapper/requirements.txt
RUN pip install -e ./cell_type_mapper
RUN pip3 install numpy==1.26.4

## Clean up
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /tmp/downloaded_packages

## Strip binary installed libraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
RUN strip /usr/local/lib/R/site-library/*/libs/*.so

##