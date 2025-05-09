% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeMap.R
\name{treeMap}
\alias{treeMap}
\title{Tree based mapping}
\usage{
treeMap(
  AIT.anndata,
  query.data,
  p = 0.8,
  low.th = 0.1,
  bootstrap = 100,
  genes.to.use = rep(TRUE, dim(AIT.anndata)[2]),
  seed = 1
)
}
\arguments{
\item{AIT.anndata}{A reference taxonomy anndata object.}

\item{query.data}{A logCPM normalized matrix to be annotated.}

\item{p}{The proportion of marker genes to include in each iteration of the mapping algorithm.}

\item{low.th}{The minimum difference in Pearson correlation to the reference cluster mean gene expression between the top-matched cluster and others. If the difference is higher than low.th, the mapping process continues; otherwise, a random branch is chosen.}

\item{bootstrap}{Number of bootstrapping runs to calculate the membership from (default = 100)}

\item{genes.to.use}{The set of genes to use for tree mapping (default is all marker genes in the tree). Can be (1) a character vector of gene names, (2) a TRUE/FALSE (logical) vector of which genes to include, or (3) a column name in AIT.anndata$var corresponding to a logical vector of variable genes. If anything is provided the set of genes used is the intersection of all marker genes and the gene set here.}

\item{seed}{Value of the seed for reproducibility}
}
\value{
Tree mapping results as a data.frame.
}
\description{
Tree based mapping
}
