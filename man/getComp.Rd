\name{getComp}
\alias{getComp}
\alias{getComp,IcaSet-method}
\alias{getComp,IcaSet,character,numeric}
\alias{getComp,IcaSet,character,numeric-method}


\title{Retrieve feature and sample values on a component stored in an \code{IcaSet} object.}
\description{This generic function retrieves, from an \code{\link{IcaSet}} object,
  the feature projections (contained in attribute \code{S}) and
  sample contributions (contained in attribute \code{A})
  corresponding to a specific component.
}
\usage{
getComp(object, level, ind)
}
\arguments{
  \item{object}{Object of class \code{\link{IcaSet}}.}
  \item{level}{Either "features" to retrieve projections contained in
    \code{S}, or "genes" to retrieve projections contained in
    \code{SByGene}.}
  \item{ind}{The index of the component to be retrieved.}
  
}
\value{
  \code{getComp}  returns a list containing two elements:
  \describe{\item{proj:}{the feature or gene
      projections on the given component,}\item{contrib:}{the sample
      contributions on the given component.}}
}
\author{Anne Biton}

\seealso{\code{\link{IcaSet-class}}}
