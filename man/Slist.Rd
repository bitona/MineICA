\name{Slist}
\alias{Slist}
\alias{SlistByGene}


\title{Retrieve feature/gene projections stored in an \code{\link{IcaSet}} object
  as a list.}
\description{These generic functions retrieve, from an IcaSet object,
  the feature and gene projections contained in
  the attribute \code{S} and \code{SByGene}  as
  a list where feature and gene IDs are preserved.
}
\usage{
Slist(object)
SlistByGene(object)
}
\arguments{
  \item{object}{Object of class \code{IcaSet}.}
}
\value{
  \code{Slist} and \code{SlistByGene} return a list whose length equals the number of
  components contained in the \code{IcaSet} object. Each element of this
  list contains a vector of feature or gene projections 
  indexed by the feature or gene IDs.
}
\author{Anne Biton}

\seealso{\code{\link{IcaSet-class}}}
