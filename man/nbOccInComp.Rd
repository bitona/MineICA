\name{nbOccInComp}
\alias{nbOccInComp}
\title{Select components the features contribute to}
\usage{
  nbOccInComp(icaSet, params, selectionByComp = NULL,
    level = c("features", "genes"), file = NULL)
}
\arguments{
  \item{icaSet}{An object of class \code{\link{IcaSet}}}

  \item{params}{An object of class
  \code{\link{MineICAParams}} containing the parameters of
  the analysis, the attribute \code{cutoffSel} is used as a
  threshold on the absolute projections to determine which
  genes contribute to the components.}

  \item{selectionByComp}{The list of components already
  restricted to the contributing genes}

  \item{level}{The attribute of \code{icaSet} to be used,
  are reported the occurences of either the
  \code{"features"} or the \code{"genes"}.}

  \item{file}{The file where the output data.frame and
  plots are written.}
}
\value{
  Returns a data.frame whose columns are: 'gene' the
  feature or gene ID, 'nbOcc' the number of components on
  which the gene contributes according to the threshold,
  'components' the indices of these components, and then
  the component indices which contain its projection
  values.
}
\description{
  For each feature/gene, this function returns the
  components they contribute to and their projection values
  across all the components.
}
\details{
  A feature/gene is considered as a contributor when its
  scaled projection value exceeds the threshold
  \code{selCutoff(icaSet)}.

  This function plots the number of times the feature/gene
  is a contributor as a function of the standard deviation
  of its expression profile.

  The created files are located in \code{genePath(params)}.
  An extensiom '.htm' and '.pdf' is respectively added to
  the \code{file} name for the data.frame and the plot
  outputs.
}
\examples{
data(icaSetCarbayo)
params <- buildMineICAParams(resPath="carbayo/")
nbOcc <- nbOccInComp(icaSet=icaSetCarbayo, params=params, level="genes", file="gene2MixingMatrix")
}
\author{
  Anne Biton
}

