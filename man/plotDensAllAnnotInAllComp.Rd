\name{plotDensAllAnnotInAllComp}
\alias{plotDensAllAnnotInAllComp}
\title{Tests if groups of samples are differently distributed on the components according and do the corresponding plots.}
\usage{
  plotDensAllAnnotInAllComp(icaSet, params, path,
    keepVar = NULL, keepComp, samples,
    legend.title_list = NULL,
    colours = params["annot2col"], doPlot = TRUE,
    pval.cutoff = params["pvalCutoff"], typeImage = "png",
    filename = NULL, onlySign = TRUE)
}
\arguments{
  \item{icaSet}{an object of class
  \code{\link[MineICA:IcaSet-class]{IcaSet}}}

  \item{params}{An object of the class
  \code{\link[MineICA:MineICAParams-class]{MineICAParams}}
  containing the parameters of the analysis}

  \item{path}{the directory where the plots will be
  located}

  \item{keepVar}{The variable labels to be considered, i.e
  a subset of (\code{varLabels(icaSet)})}

  \item{samples}{a subset of sample names available in
  \code{samplenames(icaSet)}, if NULL (default) all samples
  are used}

  \item{keepComp}{a subset of components available in
  \code{indComp(icaSet)}, if NULL (default) all components
  are used}

  \item{legend.title_list}{A list of titles for each
  component, indexed by elements of argument
  \code{keepVar}, default is NULL}

  \item{colours}{A vector of colours indexed by the
  variable levels, if missing the colours are automatically
  generated using \code{annot2Color}}

  \item{doPlot}{if TRUE (default), the plots are drawn,
  else if FALSE only the tests are performed}

  \item{pval.cutoff}{The threshold p-value for statistical
  significance}

  \item{typeImage}{The type of image file where each plot
  is saved}

  \item{filename}{A file where the results will be
  displayed in format HTML, if NULL no file is created}

  \item{onlySign}{if TRUE (default), only the significant
  results are plotted}
}
\value{
  Returns a data.frame of dimensions 'components x
  variables' containing the p-values of the non-parametric
  tests (Wilcoxon or Kruskal-Wallis tests) wich test if the
  samples groups defined by each variable are differently
  distributed on the components
}
\description{
  This function tests if the groups of samples formed by
  the variables (i.e sample annotations) are differently
  distributed on the components, in terms of contribution
  value (i.e of values in matrix \code{A(icaSet)}). The
  distribution of the groups on the components are
  represented using density plots. It is possible to
  restrict the tests and the plots to a subset of samples
  and/or components.
}
\details{
  This function writes an HTML file containing the results
  of the tests and links to the corresponding density
  plots.  One png image is created by plot and located in
  the sub-directory plots of \code{path}. Each image is
  named by index-of-component_var.png. Wilcoxon or
  Kruskal-Wallis tests are applied depending on the number
  of groups of interest from the considered annotation
  (argument \code{keepLev}).
}
\examples{
\dontrun{
## load an example of IcaSet
data(icaSetCarbayo)
## have a look at the sample annotations which are available
varLabels(icaSetCarbayo)
## create parameters, specifying the result path
params <- buildMineICAParams(resPath="carbayo/")

## trace the contributions of the samples according to their cancer stages and gender on the components
## make sure the arg 'path' exists in the directory contained in resPath(params)!
restests <- plotDensAllAnnotInAllComp(icaSet=icaSetCarbayo, keepVar=c("stage","SEX"),
                                      params=params, path="testPlotDens")
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{wilcoxOrKruskalOnA}},
  \code{\link{writeHtmlResTestsByAnnot}},
  \code{\link{plotDensOneAnnotInAllComp}}
}
\keyword{internal}

