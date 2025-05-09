% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rfTreeMapping}
\alias{rfTreeMapping}
\title{Tree-based mapping}
\usage{
rfTreeMapping(
  dend,
  refDat,
  clustersF,
  mapDat = refDat,
  p = 0.7,
  low.th = 0.15,
  bootstrap = 100,
  seed = 1
)
}
\arguments{
\item{dend}{dendrogram for mapping}

\item{refDat}{normalized data of the REFERENCE data set}

\item{clustersF}{factor indicating which cluster each cell type is actually assigned to
in the reference data set}

\item{mapDat}{normalized data of the MAPPING data set.  Default is to map the data onto itself.}

\item{p}{proportion of marker genes to include in each iteration of the mapping algorithm.}

\item{low.th}{the minimum difference in Pearson correlation required to decide on which branch
to map to. otherwise, a random branch is chosen.}

\item{bootstrap}{Number of bootstrapping runs to calculate the membership from (default = 100)}

\item{seed}{added for reproducibility}
}
\value{
a matrix of confidence scores (from 0 to 100) with rows as cells and columns
as tree node/leafs.  Values indicate the fraction of permutations in which the cell
mapped to that node/leaf using the subset of cells/genes in map_dend
}
\description{
Returns the mapping membership of each cell to each node and leaf using a
tree-based method.  This is a wrapper function for map_dend.  Includes
Minor adjustments from the function of the same name in \code{mfishtools}.
}
