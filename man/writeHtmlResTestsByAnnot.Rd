\name{writeHtmlResTestsByAnnot}
\alias{writeHtmlResTestsByAnnot}
\title{Tests if groups of samples are differently distributed on the components according and do the corresponding plots.}
\usage{
  writeHtmlResTestsByAnnot(params, icaSet, res, res2,
    nameres = "p", nameres2 = "cor", onlySign = TRUE,
    cutoff = params["pvalCutoff"],
    cutoffDir = c("<=", ">="), path, pathplot = "plots/",
    filename = NULL, typeImage = "png", caption = "",
    keepVar)
}
\arguments{
  \item{params}{An object of class
  \code{\link[MineICA:MineICAParams-class]{MineICAParams}}
  containing the parameters of the analysis}

  \item{icaSet}{An object of class
  \code{\link[MineICA:IcaSet-class]{IcaSet}}}

  \item{res}{A matrix or data.frame of dimension
  'components x variables' containing numeric values that
  quantify the association of the components with sample
  variables (e.g p-values, FDR, correlation values). This
  is the matrix used to select the significant results
  according to \code{cutoff} and \code{cutoffDir}.}

  \item{res2}{A matrix or data.frame of dimension
  'components x variables' containing numeric values that
  quantify the association of the components with sample
  annotations (e.g p-values, FDR, correlation values). It
  is only used as an additional result displayed in the
  output.}

  \item{nameres}{Name of the values contained in
  \code{res}, default is "p"}

  \item{nameres2}{Name of the values contained in
  \code{res2}, default is "cor"}

  \item{onlySign}{If TRUE (default), only the significant
  results are plotted}

  \item{cutoff}{The threshold p-value for statistical
  significance}

  \item{path}{A directory for the HTML file containing the
  p-value results}

  \item{pathplot}{A directory for the plots}

  \item{filename}{The name of the file where the results
  will be displayed in format HTML, if NULL no file is
  created}

  \item{typeImage}{The type of image file where each plot
  is saved}

  \item{caption}{The title of the HTML table}

  \item{cutoffDir}{The direction to be used with the
  cutoff: \code{"inf"} for "<=" and \code{"sup"} for ">="}

  \item{keepVar}{The variable labels to be considered, i.e
  a subset of the variables of icaSet available in
  \code{varLabels(icaSet)}.}
}
\value{
  Returns a data.frame of dimensions 'components x
  variables' containing the p-values of the non-parametric
  tests (Wilcoxon or Kruskal-Wallis tests) wich test if the
  samples groups defined by each variable are differently
  distributed on the components
}
\description{
  This internal function creates an HTML file containing a
  table of dimensions 'variables x components' with
  p-values. When a p-value is considered as significant
  according to the threshold \code{cutoff}, it is written
  in bold and filled with a link pointing to the
  corresponding plot. These plots are contained in images
  located in the path \code{pathplot}. To be identified by
  the function, the file syntax of each image file must be
  "index-of-component_colAnnot.typeImage".
}
\details{
  If argument \code{onlySign} is TRUE, then only links to
  plots that are significant according to the given
  threshold are provided.

  When \code{res2} is not missing, the values contained in
  \code{res2} are pasted to the values contained in res in
  the output array. \code{nameres} and \code{nameres2} are
  used such as every element in the ouput array contains
  two indexed values: \code{nameres}=x, \code{nameres2}=y.
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{p.adjust}}, \code{\link{qualVarAnalysis}},
  \code{\link{quantVarAnalysis}}
}
\keyword{internal}

