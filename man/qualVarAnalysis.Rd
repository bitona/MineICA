\name{qualVarAnalysis}
\alias{qualVarAnalysis}
\title{Tests association between qualitative variables and components.}
\usage{
  qualVarAnalysis(params, icaSet, keepVar,
    keepComp = indComp(icaSet),
    keepSamples = sampleNames(icaSet),
    adjustBy = c("none", "component", "variable"),
    method = "BH", doPlot = TRUE, typePlot = "density",
    addPoints = FALSE, onlySign = TRUE,
    cutoff = params["pvalCutoff"],
    colours = annot2col(params), path = "qualVarAnalysis/",
    filename = "qualVar", typeImage = "png")
}
\arguments{
  \item{params}{An object of class
  \code{\link[MineICA:MineICAParams-class]{MineICAParams}}
  providing the parameters of the analysis.}

  \item{icaSet}{An object of class
  \code{\link[MineICA:IcaSet-class]{IcaSet}}.}

  \item{keepVar}{The variable labels to be considered, must
  be a subset of \code{varLabels(icaSet)}.}

  \item{keepComp}{A subset of components, must be included
  in \code{indComp(icaSet)}. By default, all components are
  used.}

  \item{keepSamples}{A subset of samples, must be included
  in \code{sampleNames(icaSet)}. By default, all samples
  are used.}

  \item{adjustBy}{The way the p-values of the Wilcoxon and
  Kruskal-Wallis tests should be corrected for multiple
  testing: \code{"none"} if no p-value correction has to be
  done, \code{"component"} if the p-values have to be
  corrected by component, \code{"variable"} if the p-values
  have to be corrected by variable}

  \item{method}{The correction method, see
  \code{\link{p.adjust}} for details, default is
  \code{"BH"} for Benjamini & Hochberg.}

  \item{doPlot}{If TRUE (default), the plots are done, else
  only tests are performed.}

  \item{addPoints}{If TRUE, points are superimposed on the
  boxplot.}

  \item{typePlot}{The type of plot, either \code{"density"}
  or \code{"boxplot"}.}

  \item{onlySign}{If TRUE (default), only the significant
  results are plotted.}

  \item{cutoff}{A threshold p-value for statistical
  significance.}

  \item{colours}{A vector of colours indexed by the
  variable levels, if missing the colours are automatically
  generated using \code{\link{annot2Color}}.}

  \item{path}{A directory _within resPath(params)_ where
  the files containing the plots and the p-value results
  will be located. Default is \code{"qualVarAnalysis/"}.}

  \item{typeImage}{The type of image file to be used.}

  \item{filename}{The name of the HTML file containing the
  p-values of the tests, if NULL no file is created.}
}
\value{
  Returns A data.frame of dimensions 'components x
  variables' containing the p-values of the non-parametric
  tests (Wilcoxon or Kruskal-Wallis tests) wich test if the
  samples groups defined by each variable are differently
  distributed on the components.
}
\description{
  This function tests if the groups of samples formed by
  the variables are differently distributed on the
  components, in terms of contribution value (i.e of values
  in matrix \code{A(icaSet)}). The distribution of the
  samples on the components are represented using either
  density plots of boxplots. It is possible to restrict the
  tests and the plots to a subset of samples and/or
  components.
}
\details{
  This function writes an HTML file containing the results
  of the tests as a an array of dimensions 'variables *
  components' containing the p-values of the tests. When a
  p-value is considered as significant according to the
  threshold \code{cutoff}, it is written in bold and filled
  with a link pointing to the corresponding plot. One image
  is created by plot and located into the sub-directory
  "plots/" of \code{path}. Each image is named by
  index-of-component_var.png. Wilcoxon or Kruskal-Wallis
  tests are performed depending on the number of groups of
  interest in the considered variable (argument
  \code{keepLev}).
}
\examples{
## load an example of IcaSet
data(icaSetCarbayo)

## build MineICAParams object
params <- buildMineICAParams(resPath="carbayo/")

## Define the directory containing the results
dir <- paste(resPath(params), "comp2annot/", sep="")

## Run tests, make no adjustment of the p-values,
# for variable grade and components 1 and 2,
# and plot boxplots when 'doPlot=TRUE'.
qualVarAnalysis(params=params, icaSet=icaSetCarbayo, adjustBy="none", typePlot="boxplot",
                keepVar="GRADE", keepComp=1:2, path=dir, doPlot=FALSE)
}
\author{
  Anne Biton
}
\seealso{
  , \code{\link{qualVarAnalysis}}, \code{\link{p.adjust}},
  \code{link{writeHtmlResTestsByAnnot}},
  \code{wilcox.test}, \code{kruskal.test}
}

