\name{plotPosAnnotInComp}
\alias{plotPosAnnotInComp}
\title{Histograms of sample contributions for each annotation level}
\usage{
  plotPosAnnotInComp(icaSet, params,
    keepVar = varLabels(icaSet),
    keepComp = indComp(icaSet),
    keepSamples = sampleNames(icaSet), pathPlot = NULL,
    breaks = 20, colAll = "grey74", colSel, resClus,
    funClus = c("Mclust", "kmeans"), nbClus = 2,
    by = c("annot", "component"),
    typeImage = c("pdf", "png", "none"), ...)
}
\arguments{
  \item{icaSet}{An object of class \code{IcaSet}}

  \item{params}{A \code{MineICAParams} object}

  \item{keepVar}{The variable labels to be considered, i.e
  a subset of the column labels of the pheno data of icaSet
  available in (\code{varLabels(icaSet)})}

  \item{keepComp}{A subset of components available in
  \code{indComp(icaSet)}; by default, all components are
  used}

  \item{keepSamples}{A subset of samples, must be available
  in \code{sampleNames(icaSet)}; by default, all samples
  are used}

  \item{pathPlot}{A character specifying the path where the
  plots will be saved}

  \item{breaks}{The number of breaks to be used in the
  histograms}

  \item{colSel}{The colour of the histogram of the group of
  interest, default is "red"}

  \item{colAll}{The colour of the global histogram, default
  is "grey74"}

  \item{resClus}{A list containing the outputs of function
  \code{clusterSamplesByComp}, which consists of sample
  clustering applied to matrix A of argument \code{icaSet}.
  If missing, the clustering is performed by the function.}

  \item{funClus}{The clustering method to be used, either
  \code{"Mclust"} or \code{"kmeans"}. If \code{resClus} is
  not missing, equals \code{resClus$funClus}.}

  \item{nbClus}{If \code{resClus} is missing, it provides
  the number of clusters to be computed by \code{funClus},
  default is 2}

  \item{by}{Either \code{"annot"} to plot the histograms of
  each variable across all components, or
  \code{"component"} to plot the histograms for each
  component across variables. When \code{by="annot"} one
  pdf file is created by variable name, while when
  \code{annot="component"}, one pdf file is created by
  component.}

  \item{typeImage}{The type of image to be created, either
  "pdf" (default) or "png". "png" is not recommended,
  unless there are at the most 4 histograms to be plotted,
  because it does not allow to deal with multiple pages of
  plots.}

  \item{...}{Additional parameters for function
  \code{\link{hist}}}
}
\value{
  NULL
}
\description{
  This function plots the positions of groups of samples
  formed by the variables (i.e the sample annotations)
  across all the components of an object of class
  \code{\link[MineICA:IcaSet-class]{icaSet}}. For each
  variable level (e.g for each tumor stage) this function
  plots the positions of the corresponding samples (e.g the
  subset of samples having this tumor stage) within the
  histogram of the global sample contributions. The plots
  are saved in pdf file, one file is created per variable.
  The pdf files are names 'variable.pdf' and save either in
  \code{pathPlot} if specified or the current directory.
}
\details{
  The plotted values are the sample contributions across
  the components, i.e across the columns of
  \code{A(icaSet)}.

  If argument \code{resClus} is missing, the function
  computes the clustering of the samples on each component
  (i.e on each column of \code{A(icaSet)}) using
  \code{funClus} and \code{nbClus}.

  The association between the clusters and the considered
  sample group is tested using a chi-square test. The
  p-values of these tests are available in the title of
  each plot.

  When \code{by="annot"} this function plots the histograms
  of each variable across all components, to plot the
  histograms for each component across variables, please
  use \code{by="component"}.
}
\examples{
\dontrun{
## load an example of IcaSet
data(icaSetCarbayo)

## Use icaSetCarbayo, look at the available annotations
varLabels(icaSetCarbayo)

## Plot positions of samples in components according to annotations 'SEX' and 'STAGE'
# plots are saved in files SEX.pdf and STAGE.pdf created in the current directory
plotPosAnnotInComp(icaSet=icaSetCarbayo, keepVar=c("SEX","STAGE"), keepComp=1:2,  funClus="Mclust")
# specifiy arg 'pathPlot' to save the pdf in another directory, but make sure it exists before
# specifiy arg 'by="comp"' to create one pdf file per component
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{plotPosSamplesInComp}}, \code{chisq.test}
}

