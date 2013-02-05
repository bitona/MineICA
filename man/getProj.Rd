\name{getProj}
\alias{getProj}
\title{Extract projection values}
\usage{
  getProj(icaSet, ids, keepComp,
    level = c("features", "genes"))
}
\arguments{
  \item{icaSet}{An object of class \code{\link{IcaSet}}}

  \item{ids}{feature or gene IDs}

  \item{keepComp}{Index of the components to be conserved,
  must be in \code{indComp(icaSet)}}

  \item{level}{The level of projections to be extracted,
  either \code{"features"} or \code{"genes"}}
}
\value{
  A vector or a list of projection values
}
\description{
  Extract projection values of a given set of IDs on a
  subset of components.
}
\examples{
## load an example of IcaSet
data(icaSetCarbayo)

##get the projection of your favorite proliferation genes
#on all components
getProj(icaSetCarbayo, ids=c("TOP2A","CDK1","CDC20"), level="genes")

#on some components
getProj(icaSetCarbayo, ids=c("TOP2A","CDK1","CDC20"),
keepComp=c(1,6,9,12),level="genes")

##get the gene projection values on the sixth component
getProj(icaSetCarbayo, keepComp=6,level="genes")
}
\author{
  Anne Biton
}

