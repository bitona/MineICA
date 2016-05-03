\name{clusVarAnalysis}
\alias{clusVarAnalysis}
\title{Tests association between clusters of samples and variables}
\usage{
  clusVarAnalysis(icaSet, params, resClus, keepVar,
    keepComp, funClus = "",
    adjustBy = c("none", "component", "variable"),
    method = "BH", doPlot = FALSE,
    cutoff = params["pvalCutoff"],
    path = paste(resPath(params), "clus2var/", sep = ""),
    onlySign = TRUE, typeImage = "png",
    testBy = c("variable", "level"), filename)
}
\arguments{
  \item{icaSet}{An object of class
  \code{\link[MineICA:IcaSet-class]{IcaSet}}}

  \item{params}{An object of class
  \code{\link[MineICA:MineICAParams-class]{MineICAParams}}
  providing the parameters of the analysis}

  \item{resClus}{A list of numeric vectors indexed by
  sample IDs, which specifies the sample clusters. There
  must be one clustering by component of \code{icaSet}. The
  names of the list must correspond to the component
  indices.}

  \item{keepVar}{The variable labels to be considered, i.e
  a subset of the variables of icaSet available in
  \code{varLabels(icaSet)}.}

  \item{keepComp}{A subset of components available in
  \code{indComp(icaSet)} to be considered, if missing all
  components are used.}

  \item{funClus}{The name of the function used to perform
  the clustering (just for text in written files).}

  \item{adjustBy}{The way the p-values of the Wilcoxon and
  Kruskal-Wallis tests should be corrected for multiple
  testing: \code{"none"} if no p-value correction has to be
  done, \code{"component"} if the p-values have to be
  corrected by component, \code{"variable"} if the p-values
  have to be corrected by variable.}

  \item{testBy}{Chi-square tests of association can be
  performed either by \code{"variable"} (one test by
  variable, default) or by variable \code{"level"} (as many
  tests as there are annotation levels).}

  \item{method}{The correction method, see
  \code{\link{p.adjust}} for details, default if
  \code{"BH"} for Benjamini & Hochberg.}

  \item{doPlot}{If TRUE, the barplots showing the
  distribution of the annotation levels among the clusters
  are plotted and the results are provided in an HTML file
  'cluster2annot.htm', else no plot is created.}

  \item{cutoff}{The threshold for statistical
  significance.}

  \item{filename}{File name for test results, if
  \code{doPlot=TRUE} will be an HTML file else will be a
  'txt' file. If missing when \code{doPlot=TRUE}, will be
  "clusVar".}

  \item{path}{A directory _within resPath(params)_ where
  the outputs are saved if \code{doPlot=TRUE}, default is
  \code{'cluster2annot/'}.}

  \item{onlySign}{If TRUE (default), only the significant
  results are plotted.}

  \item{typeImage}{The type of image file where each plot
  is saved.}
}
\value{
  This function returns a list whose each element gives,
  for each component, the results of the association
  chi-squared tests between the clusters and the annotation
  levels.
}
\description{
  From a clustering of samples performed according to their
  contribution to each component, this function computes
  the chi-squared test of association between each variable
  level and the cluster, and summarizes the results in an
  HTML file.
}
\details{
  When \code{doPlot=TRUE}, this function writes an HTML
  file containing the results of the tests as a table of
  dimension 'variable levels x components' which contains
  the p-values of the tests. When a p-value is considered
  as significant according to the threshold \code{cutoff},
  it is written in bold and filled with a link pointing to
  the corresponding barplot displaying the distribution of
  the clusters across the levels of the variables.

  One image is created by plot and located into the
  sub-directory "plots/" of \code{path}. Each image is
  named by index-of-component_var.png
}
\examples{
## load an example of IcaSet
data(icaSetCarbayo)
## build object of class MineICAParams
params <- buildMineICAParams(resPath="carbayo/")

## cluster samples according to the columns of the mixing matrix A with kmeans in 2 groups
resClus <- clusterSamplesByComp(icaSet=icaSetCarbayo, params=params, funClus="kmeans",
                               clusterOn="A", nbClus=2)$clus

## specify directory for the function outputs (here same directory as the default one)
## this directory will be created by the function in resPath(params)
dir <- "clus2var/"

## compute chi-square tests of association, p-value are not adjusted (adjustBy="none"),
# test results are written in txt format (doPlot=FALSE and filename not missing)
resChi <- clusVarAnalysis(icaSet=icaSetCarbayo, params=params, resClus=resClus, funClus="kmeans",
                          adjustBy="none", doPlot=FALSE, path=dir, filename="clusVarTests")

\dontrun{
## compute chi-square tests of association, p-value are not adjusted (adjustBy="none"),
# write results and plots in HTML files (doPlot=TRUE)
resChi <- clusVarAnalysis(icaSet=icaSetCarbayo, params=params, resClus=resClus, funClus="kmeans",
                          path=dir, adjustBy="none", doPlot=TRUE, filename="clusVarTests")

## compute chi-square tests of association by only considering a subset of components and variables,
# adjust p-values by component  (adjustBy="component"),
# do not write results (doPlot=FALSE and filename is missing).
resChi <- clusVarAnalysis(icaSet=icaSetCarbayo, params=params, resClus=resClus, keepComp = 1:10,
                          keepVar=c("GENDER","STAGE"), funClus="kmeans", adjustBy="none",
                          doPlot=FALSE)
}
}
\author{
  Anne Biton
}
\seealso{
  \code{clusterSamplesByComp}
}

