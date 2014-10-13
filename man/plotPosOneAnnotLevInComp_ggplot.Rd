\name{plotPosOneAnnotLevInComp_ggplot}
\alias{plotPosOneAnnotLevInComp_ggplot}
\title{Plots the position of a subset of samples in the histogram of all samples using \code{\link{ggplot2}}.}
\usage{
  plotPosOneAnnotLevInComp_ggplot(annot, colAnnot, selLev,
    comp, title = NULL, colSel = "red", colAll = "grey74",
    binwidth = 0.1, geneExpr = NULL, geneRef = NULL, ...)
}
\arguments{
  \item{annot}{a data.frame of dimensions 'samples x
  annotations' with one column corresponding to the
  component to trait ("comp" column) and one column
  corresponding to the groups of interest ("interest"
  column)}

  \item{colAnnot}{the name of a column of the argument
  \code{annot} with the groups of interest}

  \item{selLev}{the name of the group of interest}

  \item{comp}{a vector of sample contributions}

  \item{colSel}{the colour of the histogram of the group of
  interest, default is "red"}

  \item{colAll}{the colour of the global histogram}

  \item{geneExpr}{a vector of values representative of the
  component, e.g the expression of the witness gene of the
  component}

  \item{geneRef}{the ID of the feature/gene \code{geneExpr}
  corresponds to, e.g the name of the witness gene}

  \item{title}{A title for the plot}

  \item{binwidth}{set the width of the bins, see
  \code{\link[ggplot2]{geom_histogram}}}

  \item{...}{other parameters given to
  \code{\link[ggplot2]{geom_histogram}}}
}
\value{
  An object of class ggplot2 containing the histogram
}
\description{
  Given a sample annotation (e.g a tumor specific stage),
  this function plots the positions of the corresponding
  samples (e.g the subset of samples having this tumor
  stage) within the histogram of the global sample
  contributions. This function is called by
  \code{\link{plotPosOneAnnotInComp_ggplot}} and is only
  dedicated to the plot of the histogram using the package
  \code{\link{ggplot2}}.
}
\author{
  Anne Biton
}
\seealso{
  \code{\link[ggplot2]{geom_histogram}}
}
\keyword{internal}

