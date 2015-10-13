\name{clusterSamplesByComp}
\alias{clusterSamplesByComp}
\title{Cluster samples from an IcaSet}
\usage{
  clusterSamplesByComp(icaSet, params,
    funClus = c("Mclust", "kmeans", "pam", "pamk", "hclust", "agnes"),
    filename, clusterOn = c("A", "S"),
    level = c("genes", "features"), nbClus,
    metric = "euclidean", method = "ward", ...)
}
\arguments{
  \item{icaSet}{An \code{IcaSet} object}

  \item{params}{A \code{MineICAParams} object}

  \item{funClus}{The function to be used for clustering,
  must be one of
  \code{c("Mclust","kmeans","pam","pamk","hclust","agnes")}}

  \item{filename}{A file name to write the results of the
  clustering in}

  \item{clusterOn}{Specifies the matrix used to apply
  clustering: \describe{ \item{\code{"A"}:}{the clustering
  is performed in one dimension, on the vector of sample
  contributions,} \item{"S":}{the clustering is performed
  on the original data restricted to the contributing
  individuals.}}}

  \item{level}{The level of projections to be used when
  \code{clusterOn="S"}, either \code{"features"} or
  \code{"genes"}.}

  \item{nbClus}{The number of clusters to be computed,
  either a single number or a numeric vector whose length
  equals the number of components. If missing (only allowed
  if \code{funClus} is one of \code{c("Mclust","pamk")})}

  \item{metric}{Metric used in \code{pam} and
  \code{hclust}, default is \code{"euclidean"}}

  \item{method}{Method of hierarchical clustering, used in
  \code{hclust} and \code{agnes}}

  \item{...}{Additional parameters required by the
  clustering function \code{funClus}.res <-
  clusterSamplesByComp(icaSet=icaSetCarbayo, params=params,
  funClus="kmeans",}
}
\value{
  A list consisting of three elements
  \describe{\item{clus:}{a list specifying the sample
  clustering for each component,}\item{resClus:}{the
  complete output of the clustering
  function,}\item{funClus:}{the function used to perform
  the clustering.}}. When \code{clusterOn="S"}, if some
  components were not used because no contributing elements
  is selected using the cutoff, the icaSet with the
  corresponding component deleted is also returned.
}
\description{
  This function allows to cluster samples according to the
  results of an ICA decomposition. One clustering is run
  independently for each component.
}
\examples{
data(icaSetCarbayo)
params <- buildMineICAParams(resPath="carbayo/", selCutoff=4)

## cluster samples according to their contributions
# using Mclust without a number of clusters
res <- clusterSamplesByComp(icaSet=icaSetCarbayo, params=params, funClus="Mclust",
                            clusterOn="A", filename="clusA")

# using kmeans
res <- clusterSamplesByComp(icaSet=icaSetCarbayo, params=params, funClus="kmeans",
                            clusterOn="A", nbClus=2, filename="clusA")
}
\author{
  Anne Biton
}
\seealso{
  \code{Mclust}, \code{kmeans}, \code{pam}, \code{pamk},
  \code{hclust}, \code{agnes}, \code{cutree}
}

