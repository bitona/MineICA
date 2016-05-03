\name{plot_heatmapsOnSel}
\alias{plot_heatmapsOnSel}
\title{Plot heatmap associated with each component}
\usage{
  plot_heatmapsOnSel(icaSet, selCutoff = 4,
    level = c("features", "genes"), samplesOrder,
    featuresOrder, selectionByComp, keepVar,
    keepComp = indComp(icaSet), doSamplesDendro = TRUE,
    doGenesDendro = TRUE,
    heatmapCol = maPalette(low = "blue", high = "red", mid = "yellow", k = 44),
    file = "", path = "", annot2col, ...)
}
\arguments{
  \item{icaSet}{The IcaSet object}

  \item{selCutoff}{A numeric threshold used to select the
  contributing genes based on their projection values. Must
  be either of length 1 and the same treshold is applied to
  all components, or of length equal to the number of
  components and one specific threshold is used for each
  component.}

  \item{samplesOrder}{A list providing the order of the
  samples, per component, to be used in the heatmaps. If
  missing, the contribution values of the samples are used
  to rank the columns of the heatmaps.}

  \item{featuresOrder}{A list providing the order of the
  genes, per component, to be used in the heatmaps. If
  missing, the projection values of the genes are used to
  rank the rows of the heatmaps.}

  \item{selectionByComp}{A list of gene projections per
  component already restricted to the contributing genes,
  if missing is computed by the function.}

  \item{level}{A character indicating which data level is
  used to plot the heatmaps: either \code{'features'} to
  represent the data at the feature levels (e.g expression
  profiles of probe sets), or \code{'genes'} to represent
  the data at the annotated-features level (e.g gene
  expression profiles).}

  \item{keepVar}{The variable labels to be considered, i.e
  a subset of the column labels of the pheno data of icaSet
  available in (\code{varLabels(icaSet)})}

  \item{keepComp}{A subset of components, must be included
  in \code{indComp(icaSet)}. By default, all components are
  used.}

  \item{doSamplesDendro}{A logical indicating whether a
  hierarchical clustering has to be performed on the data
  matrix restricted to the contributing features/genes, and
  whether the corresponding dendrogram has to be plotted,
  default is TRUE.}

  \item{doGenesDendro}{A logical indicating if the
  dendrogram of features/genes has to be plotted, default
  is FALSE.}

  \item{heatmapCol}{A list of colors used to for heatmap
  coloring (see argument \code{col} of the function
  \code{image}).}

  \item{file}{A character to add to each pdf file name.
  This function creates one file by component named
  "index-of-component_\code{file}.pdf" .}

  \item{path}{A directory for the output pdf files, must
  end with "/". Default is current directory.}

  \item{annot2col}{A vector of colours indexed by the
  levels of the variables of \code{icaSet} (i.e all the
  annotation values available in \code{pData(icaSet)}). If
  missing the colours are generated automatically using the
  function \code{annot2Color}}

  \item{...}{Additional parameters for function
  \code{heatmap.plus}}
}
\value{
  A list with one element per component, each of them being
  a list consisting of three elements:
  \describe{\item{x}{the matrix represented by the
  heatmap},\item{breaks}{the breaks used for the colours of
  the heatmap},\item{dendro}{the dendrogram}.}
}
\description{
  This function plots the heatmaps representing the
  measured values of the contributing features/genes on
  each component. It also plots the sample annotations
  above each heatmap using colours.
}
\details{
  This function restricts the data matrix of an
  \code{\link{IcaSet}} object to the contributing
  genes/features, and order features/genes and samples
  either as asked by the user or according to their values
  in the ICA decomposition.

  The heatmap is plotted using a slightly modified version
  of the function \code{heatmap.plus} from the package of
  the same name. By default in this function, the
  hierarchical clustering is calculated using the function
  \code{\link[cluster]{agnes}} with euclidean metric and
  Ward's method.
}
\examples{
\dontrun{
## load an example of IcaSet object
data(icaSetCarbayo)

## check which variables you would like to use in the heatmap
varLabels(icaSetCarbayo)
keepVar <- c("STAGE","SEX")
## Use only component 1
keepComp <- 1

## For each component, select contributing *genes* using a threshold of 2 on the absolute projection values,
## and plot heatmaps of these contributing genes by ordering genes and samples according to their contribution values
plot_heatmapsOnSel(icaSet = icaSetCarbayo, selCutoff = 2, level = "genes", keepVar = keepVar,
                   keepComp=1, doSamplesDendro = TRUE, doGenesDendro = TRUE,
                   heatmapCol = maPalette(low = "blue",high = "red", mid = "yellow", k=44),
                   file = "heatmapWithoutDendro_zval3.pdf")

## For each considered component, select contributing *features* using a threshold of 2 on the absolute projection values,
## and plot heatmaps of these contributing genes with dendrograms
plot_heatmapsOnSel(icaSet = icaSetCarbayo, selCutoff = 2, level = "features", keepVar = keepVar,
                   keepComp=1, doSamplesDendro = TRUE, doGenesDendro = TRUE,
                   heatmapCol = maPalette(low = "blue",high = "red", mid = "yellow", k=44),
                   file = "heatmapWithDendro_zval3.pdf")



}
}
\author{
  Anne Biton
}
\seealso{
  \code{heatmap.plus}, \code{\link{image}},
  \code{\link{annot2Color}},
  \code{\link{build_sortHeatmap}}
}

