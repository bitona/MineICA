\name{plotDensOneAnnotInAllComp}
\alias{plotDensOneAnnotInAllComp}
\title{Tests if groups of samples are differently distributed on the components and do the corresponding plots.}
\usage{
  plotDensOneAnnotInAllComp(icaSet, keepVar, path = NULL,
    samples, keepComp, keepLev = NULL, colours = NULL,
    legend.title = NULL, doPlot = TRUE, cutoff = 0.05,
    onlySign = TRUE, resTests)
}
\arguments{
  \item{icaSet}{an object of class
  \code{\link[MineICA:IcaSet-class]{IcaSet}}}

  \item{keepVar}{a variable label, i.e the label of a
  column of the pheno data of icaSet available in
  (\code{varLabels(icaSet)}) wich contains the groups of
  interest}

  \item{path}{the directory where the plots will be
  located}

  \item{samples}{a subset of sample names available in
  \code{samplenames(icaSet)}, if NULL (default) all samples
  are used}

  \item{keepComp}{a subset of components available in
  \code{indComp(icaSet)}, if NULL (default) all components
  are used}

  \item{keepLev}{the groups of interest, i.e the levels of
  the annotation \code{keepVar} to be considered}

  \item{colours}{A vector of colours indexed by the
  elements of \code{keepLev}, if NULL the colours are
  generated automatically using annot2Color}

  \item{legend.title}{title of the legend}

  \item{cutoff}{The threshold p-value for statistical
  significance}

  \item{doPlot}{if TRUE (default), the plots are drawn,
  else if FALSE only test results are returned}

  \item{onlySign}{if TRUE (default), only the significant
  results are plotted}

  \item{resTests}{a vector of p-values per component, if
  NULL (default) the p-values are calculated using Wilcoxon
  or Kruskal-Wallis test}
}
\value{
  Returns a data.frame of dimensions 'components x 1'
  containing the results of the non-parametric tests
  (Wilcoxon or Kruskal-Wallis tests) that test if the
  groups of interest are differently distributed on the
  components
}
\description{
  Given a variable of the phenotype data (i.e vector of
  sample annotations), this function tests if the groups of
  samples formed by this variable are differently
  distributed on the components, in terms of contribution
  values. The distribution of the groups on the components
  are represented using density plots. It is possible to
  restrict the tests and the plots to a subset of samples
  and/or components.
}
\details{
  Wilcoxon or Kruskal-Wallis tests are applied depending on
  the number of groups of interest from the considered
  annotation (argument \code{keepLev}). The plots are saved
  in individual files (one file per component) in arg
  'path' if specified or in the current directory if not
  specificied. Ech individual file is nameb
  'index-of-component_colAnnot.png.' Recall that the
  sample-contribution values are contained in
  \code{A(icaSet)}, and the sample annotations in
  \code{pData(icaSet)}.

  One png image is created by plot and located in
  \code{path}. Each image is named by
  'index-of-component_keepVar.png'.
}
\examples{
\dontrun{
## load an example of IcaSet
data(icaSetCarbayo)

## have a look at the sample annotations which are available
varLabels(icaSetCarbayo)

## with doPlot=TRUE trace the contributions of the samples according
## to their grade on the components
restests <- plotDensOneAnnotInAllComp(icaSet=icaSetCarbayo, keepVar="GRADE",
                                                doPlot=FALSE)
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{wilcoxOrKruskalOnA}},
  code{\link{writeHtmlResTestsByAnnot}},
  code{\link{wilcox.test}}, code{\link{kruskal.test}}
}
\keyword{internal}

