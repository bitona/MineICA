\name{writeRnkFiles}
\alias{writeRnkFiles}
\title{Write rnk files containing gene projections}
\usage{
  writeRnkFiles(icaSet, abs = TRUE, path)
}
\arguments{
  \item{icaSet}{An object of class \code{IcaSet}}

  \item{abs}{If TRUE (default) the absolute projection
  values are used.}

  \item{path}{The path that will contain the rnk files.}
}
\value{
  NULL
}
\description{
  Writes the gene projection values of each component in a
  '.rnk' file for GSEA.
}
\details{
  The .rnk format requires two columns, the first
  containing the gene IDs, the second containing the
  projection values. The genes are ordered by projection
  values. The files are named "index-of-component_abs.rnk"
  if \code{abs=TRUE}, or "index-of-component.rnk" if
  \code{abs=FALSE}.
}
\author{
  Anne
}

