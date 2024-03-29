% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HANN_utils.R
\name{de_stats_pair}
\alias{de_stats_pair}
\title{Compute differential expression summary statistics based on a differential results data.frame and de_param().}
\usage{
de_stats_pair(
  df,
  de.param = de_param(),
  cl.size1 = NULL,
  cl.size2 = NULL,
  select.genes = NULL,
  return.df = FALSE
)
}
\arguments{
\item{df}{A data.frame of pairwise differential expression results (i.e. from \code{score_selected_pairs()}).}

\item{de.param}{A list of differential gene expression parameters from \code{de_param()}}

\item{cl.size1}{Optional: The number of samples in the first/high cluster}

\item{cl.size2}{Optional: The number of samples in the second/low cluster}
}
\value{
A list of filtered differential expression results containing: List

\itemize{
\item{score} The deScore value, equal to the sum of the -log10(p-values) of differentially expressed genes, with a cap of 20 per gene.
\item{up.score} The deScore value for up-regulated genes.
\item{down.score} The deScore value for down-regulated genes.
\item{num} The number of differentially expressed genes
\item{up.num} The number of up-regulated genes
\item{down.num} The number of down-regulated genes
\item{genes} Gene symbols for differentially expressed genes.
\item{up.genes} Gene symbols for up-regulated genes.
\item{down.genes} Gene symbols for down-regulated genes.
\item{de.df} The df used as input, filtered for differentially expressed genes.
}
}
\description{
Compute differential expression summary statistics based on a differential results data.frame and de_param().
}
\keyword{internal}
