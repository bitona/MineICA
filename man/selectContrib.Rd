\name{selectContrib}
\alias{selectContrib}
\alias{selectContrib,IcaSet-method}
\alias{selectContrib,IcaSet,numeric,character-method}
\alias{selectContrib,list,numeric,ANY}
\alias{selectContrib,list,numeric,ANY-method}


\title{Select contributing features/genes}
\usage{
  selectContrib(object, cutoff, level, ...)
}
\arguments{
  \item{object}{Either an \code{IcaSet} object, or a list
  of projection vectors, e.g the list of feature or gene
  projections on each component.}

  \item{cutoff}{The threshold according to which the
    elements will be selected. Must be either of length 1 and the same treshold is
  applied to all components, or of length equal to the number of
  components in order to use a specific threshold for each component.}
  
  \item{level}{The level of the selection: either \code{"genes"} to select
  contributing genes using \code{SByGene(icaSet)}, or \code{"features"}
  to select contributing features using \code{S(icaSet)}. }
  \item{...}{...}
}
\value{
  A list of projections restricted to the elements that are
  higher than \code{cutoff}.
}
\description{
  This function selects elements whose absolute scaled
  values exceed a given threshold.
}
\details{
  Each vector is first scaled and then only elements with
  an absolute scaled value higher than \code{cutoff} are
  kept.
}
\examples{
\dontrun{
## load an example of icaSet
data(icaSetCarbayo)

##### =========
####  When arg 'object' is an IcaSet object
##### =========

## select contributing genes
selectContrib(object=icaSetCarbayo, cutoff=3, level="genes")

## select contributing features
selectContrib(object=icaSetCarbayo, cutoff=3, level="features")

##### =========
#### When arg 'object' is a list
##### =========
 c1 <- rnorm(100); names(c1) <- 100:199
 c2 <- rnorm(100); names(c2) <- 1:99
 selectContrib(object=list(c1,c2), cutoff= 0.5)

## select contributing features
contribFlist <- selectContrib(Slist(icaSetCarbayo), 3)

## select contributing genes
contribGlist <- selectContrib(SlistByGene(icaSetCarbayo), 3)


}
}
\author{
  Anne Biton
}

