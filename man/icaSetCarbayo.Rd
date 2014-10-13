\docType{data}
\name{icaSetCarbayo}
\alias{icaSetCarbayo}
\title{IcaSet-object containing a FastICA decomposition of gene expression microarrray-based data of bladder cancer samples.}
\description{
  Object of class \code{\link{IcaSet}} containing an ICA
  decomposition calculated by the FastICA algorithm
  (through matlab function "icasso") on bladder cancer
  expression data measured on HG-U133A Affymetrix
  microarrays. The original expression data were normalized
  with MAS5 by the authors of the paper followed by
  log2-transformation. ICA was run on the dataset
  restricted to the 10000 most variable probe sets (based
  on IQR values). 10 components were computed. Only probe
  sets/genes having an absolute projection higher than 3
  are stored in this object.
}
\author{
  Anne Biton
}
\references{
  \url{http://jco.ascopubs.org/content/24/5/778/suppl/DC1}
}
\keyword{datasets}

