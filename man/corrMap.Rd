% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corrMap.R
\name{corrMap}
\alias{corrMap}
\title{Correlation based mapping}
\usage{
corrMap(
  AIT.anndata,
  query.data,
  genes.to.use = NULL,
  normalize.if.needed = TRUE
)
}
\arguments{
\item{AIT.anndata}{A reference taxonomy anndata object.}

\item{query.data}{A logCPM normalized matrix to be annotated.}

\item{genes.to.use}{The set of genes to use for correlation calculation and/or Seurat integration (default is the highly_variable_genes associated with the current mode). Can be (1) a character vector of gene names, (2) a TRUE/FALSE (logical) vector of which genes to include, or (3) a column name in AIT.anndata$var corresponding to a logical vector of variable genes.}

\item{normalize.if.needed}{Should query.data be automatically log-normalized if it contains exceedingly large values (>30). Default = TRUE.}
}
\value{
Correlation mapping results as a data.frame.
}
\description{
Correlation based mapping
}
