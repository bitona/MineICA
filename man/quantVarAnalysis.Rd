\name{quantVarAnalysis}
\alias{quantVarAnalysis}
\title{Correlation between variables and components.}
\usage{
  quantVarAnalysis(params, icaSet, keepVar,
    keepComp = indComp(icaSet),
    keepSamples = sampleNames(icaSet),
    adjustBy = c("none", "component", "variable"),
    method = "BH", typeCor = "pearson", doPlot = TRUE,
    onlySign = TRUE, cutoff = 0.4,
    cutoffOn = c("cor", "pval"), colours,
    path = "quantVarAnalysis/", filename = "quantVar",
    typeImage = "png")
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

  \item{onlySign}{If TRUE (default), only the significant
  results are plotted.}

  \item{cutoff}{A threshold p-value for statistical
  significance.}

  \item{cutoffOn}{The value the cutoff is applied to,
  either "cor" for correlation or "pval" for p-value}

  \item{typeCor}{the type of correlation to be used, one of
  \code{c("pearson","spearman","kendall")}.}

  \item{colours}{A vector of colours indexed by the
  variable levels, if missing the colours are automatically
  generated using \code{\link{annot2Color}}.}

  \item{path}{A directory _within resPath(params)_ where
  the files containing the plots and the p-value results
  will be located. Default is \code{"quantVarAnalysis/"}.}

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
  This function tests if numeric variables are correlated
  with components.
}
\details{
  This function writes an HTML file containing the
  correlation values and test p-values as a an array of
  dimensions 'variables * components' containing the
  p-values of the tests. When a p-value is considered as
  significant according to the threshold \code{cutoff}, it
  is written in bold and filled with a link pointing to the
  corresponding plot. One image is created by plot and
  located into the sub-directory "plots/" of \code{path}.
  Each image is named by index-of-component_var.png.
}
\examples{
## load an example of IcaSet
data(icaSetCarbayo)

# build MineICAParams object
params <- buildMineICAParams(resPath="carbayo/")

# Define the directory containing the results
dir <- paste(resPath(params), "comp2annottest/", sep="")

# Check which variables are numeric looking at the pheno data, here only one  -> AGE
# pData(icaSetCarbayo)

## Perform pearson correlation tests and plots association corresponding
# to correlation values larger than 0.2
quantVarAnalysis(params=params, icaSet=icaSetCarbayo, keepVar="AGE", keepComp=1:2,
                 adjustBy="none", path=dir, cutoff=0.2, cutoffOn="cor")

\dontrun{
## Perform Spearman correlation tests and do scatter plots for all pairs
quantVarAnalysis(params=params, icaSet=icaSetCarbayo, keepVar="AGE", adjustBy="none", path=dir,
                 cutoff=0.1, cutoffOn="cor", typeCor="spearman", onlySign=FALSE)

## Perform pearson correlation tests and plots association corresponding
# to p-values lower than 0.05 when 'doPlot=TRUE'
quantVarAnalysis(params=params, icaSet=icaSetCarbayo, keepVar="AGE", adjustBy="none", path=dir,
                 cutoff=0.05, cutoffOn="pval", doPlot=FALSE)
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{qualVarAnalysis}}, \code{\link{p.adjust}},
  \code{link{writeHtmlResTestsByAnnot}}, \code{code}
}

