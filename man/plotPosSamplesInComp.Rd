\name{plotPosSamplesInComp}
\alias{plotPosSamplesInComp}
\title{Histograms of sample subsets}
\usage{
  plotPosSamplesInComp(samplesByGroup, labGroups = NULL,
    icaSet, keepComp = indComp(icaSet), file = NULL,
    breaks = 20, colAll = "grey74", colSel = "red",
    titlesup = "", resClus,
    funClus = c("Mclust", "kmeans"), ...)
}
\arguments{
  \item{samplesByGroup}{A list whose elements are vector of
  sample names, these sample names must be available in
  \code{sampleNames(icaSet)}. The list should be indexed by
  the name of the corresponding groups.}

  \item{labGroups}{A vector of group names, will be used to
  add names to \code{sampleByGroup} if
  \code{names(samplesByGroup)} is NULL.}

  \item{icaSet}{An object of class
  \code{\link[MineICA:class-IcaSet]{IcaSet}}}

  \item{keepComp}{A subset of components available in
  \code{indComp(icaSet)}, if NULL (default) all components
  are used}

  \item{file}{A pdf file}

  \item{breaks}{The number of breaks to be used in the
  histograms}

  \item{colSel}{The colour of the histogram of the group of
  interest, default is "red"}

  \item{colAll}{The colour of the global histogram, default
  is "grey74"}

  \item{resClus}{A list containing the outputs of function
  \code{\link{clusterSamplesByComp}}, which consists of
  results of clustering applied to matrix A of argument
  \code{icaSet}.}

  \item{funClus}{Specifies the clustering method used,
  either \code{"Mclust"} or \code{"kmeans"}. If
  \code{resClus} is not missing, equals
  \code{resClus$funClus}.}

  \item{titlesup}{Additional title for the histograms}

  \item{...}{Additional parameters for function
  \code{\link{hist}}}
}
\value{
  NULL
}
\description{
  This function plots the positions of several groups of
  samples across all the components of an
  \code{\link[MineICA:class-IcaSet]{icaSet}} object.
}
\details{
  For each subgroup of samples this function plots their
  positions within the histogram of the global sample
  contributions.

  The values of interest are the sample contributions
  across the components, i.e across the columns
  \code{A(icaSet)}.

  If argument \code{resClus} is not missing, the
  association between the clusters and the sub-groups of
  samples is tested using a chi-square test. The p-values
  of these tests are available in the title of each plot.
}
\examples{
\dontrun{
## load an example of IcaSet
data(icaSetCarbayo)

## selection of sample groups according to annotations STAGE
samplesByGroup <- lapply(split(pData(icaSetCarbayo),pData(icaSetCarbayo)[c("STAGE")]), rownames)
# select groups including at least 2 samples
samplesByGroup <- samplesByGroup[which(unlist(lapply(samplesByGroup,length))>1)]

## clustering of samples according to A using Mclust imposing two Gaussian
resClus <- clusterSamplesByComp(icaSet=icaSetCarbayo,funClus="Mclust", nbClus=2, clusterOn="A")

## Plot positions of the groups in 5th component
pdf(file="stageOnIC5.pdf", height = 8.267717, width = 29.7/2.54, paper = 'a4r', title="stageOnIC5")
plotPosSamplesInComp(samplesByGroup=samplesByGroup, icaSet=icaSetCarbayo, funClus="Mclust",
                     resClus = resClus, keepComp=5)
dev.off()
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{hist}},
  \code{\link[MineICA:class-IcaSet]{IcaSet}}
}
\keyword{internal}

