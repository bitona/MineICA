\name{Alist}
\alias{Alist}


\title{Retrieve sample contributions stored in an \code{\link{IcaSet}} object
  as a list.}
\description{This generic function retrieves, from an IcaSet object,
  the sample contributions contained in
  the attribute \code{A} as
  a list where sample IDs are preserved.
}
\usage{
Alist(object)
}
\arguments{
  \item{object}{Object of class \code{IcaSet}.}
}
\value{
  \code{Alist} returns a list whose length equals the number of
  components contained in the \code{IcaSet} object. Each element of this
  list contains a vector of sample contributions 
  indexed by the sample IDs.
}
\author{Anne Biton}

\seealso{\code{\link{IcaSet-class}}}
