\name{selectWitnessGenes}
\alias{selectWitnessGenes}
\title{selectWitnessGenes}
\usage{
  selectWitnessGenes(icaSet, params,
    level = c("genes", "features"), maxNbOcc = 1,
    selectionByComp = NULL)
}
\arguments{
  \item{icaSet}{An object of class \code{\link{IcaSet}}}

  \item{params}{An object of class
  \code{\link{MineICAParams}} containing the parameters of
  the analysis, the attribute \code{cutoffSel} is used as
  the threshold.}

  \item{level}{The attribute of \code{icaSet} to be used,
  the witness elements will be either selected within the
  \code{"features"} or the \code{"genes"}}

  \item{maxNbOcc}{The maximum number of components where
  the genes can have an absolute projection value higher
  than \code{cutoffSel(params)} in order to be selected.}

  \item{selectionByComp}{The list of components already
  restricted to the contributing genes}
}
\value{
  This function returns a vector of IDs.
}
\description{
  This function selects a gene per component.
}
\details{
  Selects as feature/gene witness, for each component, the
  first gene whose absolute projection is greater than a
  given threshold in at the most \code{maxNbOcc}
  components. These witnesses can then be used as
  representatives of the expression behavior of the
  contributing genes of the components.

  When a feature/gene respecting the given constraints is
  not found, \code{maxNbOcc} is incremented of one until a
  gene is found.
}
\examples{
## load an example of IcaSet
data(icaSetCarbayo)

## define parameters: features or genes are considered to be contributor
# when their absolute projection value exceeds a threshold of 4.
params <- buildMineICAParams(resPath="carbayo/", selCutoff=4)

## selection, as gene witnesses, of the genes whose absolute projection is greater than 4
# in at the most one component. I.e, a gene is selected as a gene witness of a component
# if he has a large projection on this component only.
selectWitnessGenes(icaSet=icaSetCarbayo, params=params, level="genes", maxNbOcc=1)

## selection, as gene witnesses, of the genes whose absolute projection is greater than 4
# in at the most two components.
# I.e, a gene is selected as a gene witness of a given component if he has a large projection
# in this component and at the most another.
selectWitnessGenes(icaSet=icaSetCarbayo, params=params, level="genes", maxNbOcc=2)
}
\author{
  Anne Biton
}

