\name{annot2Color}
\alias{annot2Color}
\title{Association of a colour with each annotation level}
\usage{
  annot2Color(annot)
}
\arguments{
  \item{annot}{a data.frame containing the sample
  annotations (of dimension 'samples x annotations').}
}
\value{
  A vector of colours indexed by the annotation levels.
}
\description{
  Given a data.frame consisting of sample annotations, this
  function returns a vector which gives a colour per
  annotation level.
}
\details{
  Arbitrary colours are attributed to some specific
  annotations met by the author, and for the remaining
  annotation levels, the colours are attributed using
  packages \code{RColorBrewer} and \code{rcolorspace}.
}
\author{
  Anne Biton
}

