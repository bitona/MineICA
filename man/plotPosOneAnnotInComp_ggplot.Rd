\name{plotPosOneAnnotInComp_ggplot}
\alias{plotPosOneAnnotInComp_ggplot}
\title{Tests if groups of samples are differently distributed on the components and do the corresponding plots.}
\usage{
  plotPosOneAnnotInComp_ggplot(icaSet, params, colAnnot,
    keepLev = NULL, keepComp, samples, colAll = "grey74",
    binwidth = 0.1, addExpr = TRUE, file = NULL, ...)
}
\arguments{
  \item{icaSet}{An object of class \code{IcaSet}}

  \item{params}{An object of the class \code{MineICAParams}
  containing the parameters of the analysis}

  \item{colAnnot}{a variable label, i.e one of the
  variables available in (\code{varLabels(icaSet)})
  containing the groups of interest}

  \item{samples}{a subset of sample names available in
  \code{samplenames(icaSet)}, if NULL (default) all samples
  are used}

  \item{keepComp}{a subset of components available in
  \code{indComp(icaSet)}, if NULL (default) all components
  are used}

  \item{keepLev}{the groups of interest, i.e the levels of
  the variable \code{colAnnot} to be considered}

  \item{colAll}{The colour of the global histogram, default
  is "grey74"}

  \item{file}{the file where the histograms will be
  plotted}

  \item{addExpr}{if TRUE (default) the expression profiles
  of the witness genes of each component are added below
  the plot}

  \item{binwidth}{binwidth of the histogram (default is
  0.1)}

  \item{...}{other parameters for geom_histogram function
  from ggplot2 package}
}
\value{
  NULL
}
\description{
  Given a variable of the phenoData, this function tests if
  the groups of samples formed by this variable are
  differently distributed, in terms of contribution value
  (i.e of values in matrix \code{A(icaSet)}), on the
  components. The distribution of the groups on the
  components are represented using density plots. It is
  possible to restrict the tests and the plots to a subset
  of samples and/or components.
}
\details{
  Wilcoxon or Kruskal-Wallis tests are applied depending on
  the number of groups of interest from the considered
  annotation (argument \code{keepLev}). One png image is
  created by plot and located in \code{path}. Each image is
  named by component-of-component_colAnnot.png.
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{plotPosOneAnnotLevInComp_ggplot}},
  \code{\link{geom_histogram}}
}
\keyword{internal}

