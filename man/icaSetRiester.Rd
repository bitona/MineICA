\docType{data}
\name{icaSetRiester}
\alias{icaSetRiester}
\title{IcaSet-object containing a FastICA decomposition of gene expression microarrray-based data of bladder cancer samples.}
\description{
  Object of class \code{\link{IcaSet}} containing an ICA
  decomposition calculated by the FastICA algorithm
  (through matlab function "icasso") on gene expression
  data of urothelial tumors. measured on a HG-U133-plus2
  Affymetrix microarrays. It contains 20 independent
  components. The original expression data contain 93 tumor
  samples, were normalized with GCRMA with
  log2-transformation, and are restricted to the 10000 most
  variable probe sets.
}
\author{
  Anne Biton
}
\references{
  \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31684}
}
\keyword{datasets}

