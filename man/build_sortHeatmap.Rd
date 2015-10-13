\name{build_sortHeatmap}
\alias{build_sortHeatmap}
\title{Build the heatmap matrices}
\usage{
  build_sortHeatmap(icaSet, selCutoff, selectionByComp,
    level = c("features", "genes"), samplesOrder,
    featuresOrder)
}
\arguments{
  \item{icaSet}{The IcaSet object}

  \item{selCutoff}{The threshold used to select the
  contributing features/genes based on their projection
  values. Must be either of length 1 and the same treshold
  is applied to all components, or of length equal to the
  number of components and one specific threshold is used
  for each component.}

  \item{selectionByComp}{The list of gene projections per
  components already restricted to the contributing genes}

  \item{level}{A character indicating which data level is
  used to plot the heatmaps: 'features' to plot measured
  feature levels (e.g probe sets expression values),
  'genes' to plot measured gene values (e.g gene expression
  values).}

  \item{samplesOrder}{A list providing the order of the
  samples, per component, to be used in the heatmaps. If
  NULL, the contribution values of the samples are used to
  rank the columns of the heatmaps.}

  \item{featuresOrder}{A list providing the features or
  genes order, per component, to be used in the heatmaps.
  If NULL, the projection values of the genes are used to
  rank the rows of the heatmaps.}
}
\value{
  A list of matrices
}
\description{
  This function returns the matrices that will be used to
  plot the heatmaps of each component. It restricts the
  data matrix of the
  \code{\link[MineICA:IcaSet-class]{icaSet}} object to the
  contributing genes/features, and order the features/genes
  and samples.
}
\details{
  This function is called by function
  \code{\link{plot_heatmapsOnSel}} and is not likely to be
  called alone.
}
\author{
  Anne Biton
}
\keyword{internal}

