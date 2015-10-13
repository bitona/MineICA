\name{plotDens2classInComp_plotOnly}
\alias{plotDens2classInComp_plotOnly}
\title{Plots the densities or boxplots of the component contributions using \code{\link{ggplot2}}.}
\usage{
  plotDens2classInComp_plotOnly(annot, colAnnot, global,
    keepLev, comp.label = NULL, colours,
    legend.title = NULL, pval, test, title.add = NULL,
    data_ref = NULL, geneExpr = NULL, geneRef = NULL,
    ylab = NULL, trace_globalExpression = FALSE,
    trace_groupExpression = TRUE,
    typePlot = c("density", "boxplot"), addPoints = FALSE)
}
\arguments{
  \item{annot}{a data.frame of dimensions 'samples x
  annotations' with one column corresponding to the
  component to trait ("comp" column) and one column
  corresponding to the groups of interest ("interest"
  column)}

  \item{colAnnot}{the name of a column of the argument
  \code{annot} with the groups of interest}

  \item{global}{a vector with the global distribution, e.g
  the contribution values of all samples on the component}

  \item{keepLev}{the groups of interest, i.e the levels of
  the annotation \code{colAnnot} to be considered}

  \item{comp.label}{the label of the component}

  \item{colours}{a vector of colours indexed by the names
  of the groups of interest}

  \item{legend.title}{the title of the legend, if NULL
  (default) colAnnot is used}

  \item{pval}{the p-value of the test, will be written in
  the title}

  \item{test}{name of test that gave the p-value}

  \item{title.add}{a title to add to the automatically
  generated title}

  \item{data_ref}{a data.frame similar to the argument
  \code{annot} but restricted to a set of reference
  samples}

  \item{geneExpr}{a vector of values representative of the
  component, e.g the expression of the witness gene of the
  component}

  \item{geneRef}{the ID of the feature/gene \code{geneExpr}
  corresponds to, e.g the name of the witness gene}

  \item{ylab}{A label for the y-axis (character)}

  \item{trace_globalExpression}{if TRUE, geneExpr is
  plotted below the graph as a set of points whose colour
  is representative of the amount of expression, default is
  FALSE}

  \item{trace_groupExpression}{if TRUE (default), geneExpr
  is plotted below the graph, by group, as a set of points
  whose colour is representative of the amount of
  expression}

  \item{typePlot}{The type of plot, either "density" or
  "boxplot"}

  \item{addPoints}{If TRUE, points are superimposed on the
  boxplots}
}
\value{
  A
}
\description{
  This internal function is called by
  \code{\link{plotDensOneAnnotInAllComp}} and
  \code{\link{qualVarAnalysis}} and is dedicated to the
  plot of the densities or boxplots using the package
  \code{\link{ggplot2}}.
}
\author{
  Anne Biton
}
\seealso{
  \code{\link[ggplot2]{geom_density}},
  \code{\link[ggplot2]{geom_boxplot}},
  \code{\link[ggplot2]{geom_point}}
}
\keyword{internal}

