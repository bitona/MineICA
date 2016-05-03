\name{indComp}
\alias{indComp}
\alias{compNames}
\alias{witGenes}
\alias{indComp<-}
\alias{indComp<-,IcaSet-method}
\alias{compNames<-}
\alias{compNames<-,IcaSet-method}
%\alias{sampleNames<-}
%\alias{sampleNames<-,IcaSet-method}
\alias{witGenes<-}
\alias{witGenes<-,IcaSet-method}
\alias{setIndComp}
\alias{setLabelsComp}
\alias{setWitGenes}
\alias{getIndComp}
\alias{getLabelsComp}
\alias{getWitGenes}


\title{Retrieve and set component labels, indices, and witness genes from IcaSet}
\description{
  These generic functions access and set the attributes \code{compNames}, \code{indComp} and \code{witGenes}
   stored in an object of class \code{IcaSet}.
}
\usage{
indComp(object)
indComp(object) <- value
compNames(object)
compNames(object) <- value
witGenes(object)
witGenes(object) <- value
}
\arguments{
  \item{object}{object of class \code{IcaSet}}
  \item{value}{Numeric vector for \code{indComp}, character vector for \code{compNames} and
    \code{witGenes}, with length equal to \code{ncol(A(object))}
    and containing: component indices (for
    \code{indComp}), labels (for \code{compNames}), or gene witness IDs (for \code{witGenes}).}
}

\value{
  \code{indComp} returns a numeric vector containing  component indices;
  \code{compNames} returns a character vector containing component labels;
  \code{witGenes} returns a character vector containing  witness genes IDs.
  
}
\author{Anne Biton}

