\name{readS}
\alias{readS}
\title{read S}
\usage{
  readS(Sfile, datfile, dat, annot = TRUE)
}
\arguments{
  \item{Sfile}{The file which contains the matrix of
  feature projections. It must be a txt file where the
  separator is \code{white space}, that is one or more
  spaces, tabs, newlines or carriage returns.}

  \item{datfile}{The file which contains the matrix (of
  dimension features x samples) based on which the matrix S
  was calculated. It must be a txt file where the separator
  is \code{white space}, that is one or more spaces, tabs,
  newlines or carriage returns.}

  \item{dat}{The data based on which the matrix A was
  calculated (features x samples)}

  \item{annot}{TRUE (default) if the Afile contains
  rownames of matrix A, FALSE if the rownames has to be
  extracted from dat}
}
\value{
  This function returns a matrix of dimension features x
  components with rownames filled with feature IDs.
}
\description{
  This function reads and annotates matrix S.
}
\details{
  The matrix dat must be the one on which the matrix S was
  calculated. It is assumed that the number of components
  is lower than the number of features, the matrix will be
  transposed to have dimension 'features x components'
  according to this assumption. If \code{annot} is FALSE,
  rownames of dat are used to annotate rownames of S.
}
\author{
  Anne Biton
}
\keyword{internal}

