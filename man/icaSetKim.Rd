\docType{data}
\name{icaSetKim}
\alias{icaSetKim}
\title{IcaSet-object containing a FastICA decomposition of gene expression microarrray-based data of bladder cancer samples.}
\description{
  Object of class \code{\link{IcaSet}} containing an ICA
  decomposition calculated by the FastICA algorithm
  (through matlab function "icasso") on bladder cancer
  expression data measured on illumina Human-6 BeadChip,
  version 2. It contains 20 independent components. The
  original expression data contain 165 tumor samples, were
  normalized by the authors of the paper with Illumina
  BeadStudio software using Quantile normalization and log2
  transformation, and are restricted to the 10000 most
  variable probe sets.
}
\author{
  Anne
}
\references{
  \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13507}
}
\keyword{datasets}

