\name{clusterSamplesByComp_multiple}
\alias{clusterSamplesByComp_multiple}
\title{Cluster samples from an IcaSet}
\usage{
  clusterSamplesByComp_multiple(icaSet, params,
    funClus = c("Mclust", "kmeans", "pam", "pamk", "hclust", "agnes"),
    filename, clusterOn = c("A", "S"),
    level = c("genes", "features"), nbClus,
    metric = "euclidean", method = "ward", ...)
}
\arguments{
  \item{icaSet}{An \code{IcaSet} object}

  \item{params}{A \code{MineICAParams} object}

  \item{funClus}{The function to be used for clustering,
  must be several of
  \code{c("Mclust","kmeans","pam","pamk","hclust","agnes")}}

  \item{filename}{A file name to write the results of the
  clustering in}

  \item{clusterOn}{Specifies the matrix used to apply
  clustering, can be several of: \describe{
  \item{\code{"A"}:}{the clustering is performed in one
  dimension, on the vector of sample contributions,}
  \item{\code{"S"}:}{the clustering is performed on the
  original data restricted to the contributing
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
  clustering function \code{funClus}.}
}
\value{
  A list consisting of three elements
  \describe{\item{clus:}{a data.frame specifying the sample
  clustering for each component using the different ways of
  clustering,}\item{resClus:}{the complete output of the
  clustering function(s),}\item{comparClus:}{the adjusted
  Rand indices, used to compare the clusterings obtained
  for a same component.}}
}
\description{
  This function allows to cluster samples according to the
  results of an ICA decomposition. Several clustering
  functions and several levels of data for clustering can
  be performed by the function.
}
\details{
  One clustering is run independently for each component.
}
\examples{
data(icaSetCarbayo)
params <- buildMineICAParams(resPath="carbayo/", selCutoff=3)

## compare kmeans clustering applied to A and data restricted to the contributing genes
## on components 1 to 3
res <- clusterSamplesByComp_multiple(icaSet=icaSetCarbayo[,,1:3], params=params, funClus="kmeans",
                                     nbClus=2, clusterOn=c("A","S"), level="features")
head(res$clus)
}
\author{
  Anne
}
\seealso{
  \code{Mclust}, \code{adjustedRandIndex}, \code{kmeans},
  \code{pam}, \code{pamk}, \code{hclust}, \code{agnes},
  \code{cutree}
}

