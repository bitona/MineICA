\name{A}

\alias{S}
\alias{getS}
\alias{SByGene}
\alias{getSByGene}

\alias{A}
\alias{getA}
\alias{A<-}
\alias{S<-}
\alias{SByGene<-}

\alias{setA<-}
\alias{setS<-}
\alias{setSByGene<-}

\alias{A,IcaSet-method}
\alias{A<-,IcaSet,data.frame-method} 
\alias{A<-,IcaSet-method} 
\alias{S,IcaSet-method}
\alias{S<-,IcaSet,data.frame-method} 
\alias{S<-,IcaSet-method} 
\alias{SByGene,IcaSet-method}
\alias{SByGene<-,IcaSet,data.frame-method} 
\alias{SByGene<-,IcaSet-method}

\alias{nbComp}
\alias{nbComp,IcaSet-method}

\alias{getA,IcaSet-method}
\alias{setA,IcaSet-method}
\alias{getS,IcaSet-method}
\alias{setS,IcaSet-method}
\alias{getSByGene,IcaSet-method}
\alias{setSByGene,IcaSet-method}



% \alias{[}
% \alias{[,ANY,ANY,IcaSet-method}
% \alias{[,IcaSet,ANY-method} 
% \alias{[,IcaSet,ANY,ANY-method} 
% \alias{[,IcaSet,ANY,ANY,ANY-method} 
% \alias{[<-}
% \alias{[<-,IcaSet,ANY,ANY,ANY,ANY-method} 
% \alias{[<-,IcaSet,ANY,ANY,ANY-method} 
% \alias{[<-,IcaSet,ANY,ANY-method} 


\title{Retrieve and set Source S and Mixing matrix A from IcaSet}
\description{
  These generic functions access and set the attributes \code{S, SByGene} and \code{A}
   stored in an object of class \code{IcaSet}.
}
\usage{
S(object)
S(object) <- value
SByGene(object)
SByGene(object) <- value
A(object)
A(object) <- value
nbComp(object)
}
\arguments{
  \item{object}{object of class \code{IcaSet}}
  \item{value}{Data.frame with rows representing: features (for
    \code{S}), genes (for \code{SByGene}), or samples (for \code{A}) and
    columns representing components.}
}

\value{
  \code{S} returns a data.frame containing feature projection values;
  \code{SByGene} returns a data.frame containing gene projection values;
  \code{A} returns a data.frame containing  sample contribution values.
  \code{nbComp} returns the number of components, i.e the number of
  columns of \code{A}.

}
\author{Anne Biton}


