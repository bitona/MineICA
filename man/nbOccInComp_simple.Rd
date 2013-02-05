\name{nbOccInComp_simple}
\alias{nbOccInComp_simple}
\title{nbOccInComp_simple}
\usage{
  nbOccInComp_simple(icaSet, params,
    selectionByComp = NULL, level = c("features", "genes"))
}
\arguments{
  \item{icaSet}{An object of class \code{\link{IcaSet}}.}

  \item{params}{An object of class
  \code{\link{MineICAParams}} containing the parameters of
  the analysis. \code{cutoffSel(params)} is used as a
  threshold on the absolute projections to select the
  contributing features/genes.}

  \item{selectionByComp}{The list of components already
  restricted to the contributing features/genes (each
  element is a vector of projection values indexed by
  features or genes).}

  \item{level}{The attribute of \code{icaSet} to be used,
  the occurences of either the \code{"features"} (using
  \code{S(icaSet)}) or the \code{"genes"} (using
  \code{SByGene(icaSet)}) will be reported.}
}
\value{
  Returns a data.frame whose columns are: \code{gene} the
  feature or gene IDs, \code{nbOcc} the number of
  components the gene contributes to, \code{components} the
  indices of those components.
}
\description{
  For each feature/gene, this function returns the indices
  of the components they contribute to.
}
\examples{
data(icaSetCarbayo)
params <- buildMineICAParams(resPath="carbayo/")
nbOcc <- MineICA:::nbOccInComp_simple(icaSet=icaSetCarbayo, params=params, level="genes")
}
\author{
  Anne Biton
}
\keyword{internal}

