\name{hypergeoAn}
\alias{hypergeoAn}
\title{Runs an enrichment analysis per component using package \code{\link{GOstats}}.}
\usage{
  hypergeoAn(icaSet, params,
    path = paste(resPath(params), "GOstatsEnrichAnalysis/", sep = "/"),
    SlistSel, hgCutoff = 0.01, db = "go", onto = "BP",
    cond = TRUE, universe, entrez2symbol)
}
\arguments{
  \item{icaSet}{An object of class \code{IcaSet}}

  \item{params}{An object of class
  \code{\link[MineICA:MineICAParams-class]{MineICAParams}}
  containing the parameters of the analysis}

  \item{path}{The path where results will be saved}

  \item{SlistSel}{A list of contributing gene projection
  values per component.  Each element of the list
  corresponds to a component and is restricted to the
  features or genes exceeding a given threshold. If
  missing, is computed by the function.}

  \item{hgCutoff}{The p-value threshold}

  \item{db}{The database to be used (\code{"GO"} or
  \code{"KEGG"})}

  \item{onto}{A character specifying the GO ontology to
  use. Must be one of \code{"BP"}, \code{"CC"}, or
  \code{"MF"}, see
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.
  Only used when argument \code{db} is \code{"GO"}.}

  \item{cond}{A logical indicating whether the calculation
  should conditioned on the GO structure, see
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.}

  \item{universe}{The universe for the hypergeometric
  tests, see
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.}

  \item{entrez2symbol}{A vector of all gene Symbols
  involved in the analysis indexed by their Entrez Gene
  IDs. It is only used when \code{annotation(params)} is
  empty, and allows to associate gene sets to Symbols.}
}
\description{
  Runs an enrichment analysis of the contributing genes
  associated with each component, using the function
  \code{hyperGTest} of package
  \code{\link[GOstats]{GOstats}}. The easiest way to run
  enrichment analysis is to use function
  \code{\link{runEnrich}}.
}
\details{
  An annotation package must be available in
  \code{annotation(icaSet)} to provide the contents of the
  gene sets. If none corresponds to the technology you deal
  with, please choose the org.*.eg.db package according to
  the organism (for example org.Hs.eg.db for Homo sapiens).
  Save results of the enrichment tests in a '.rda' file
  located in
  \code{path}/\code{db}/\code{onto}/\code{zvalCutoff(params)}.
}
\examples{
\dontrun{
## load an example of IcaSet
data(icaSetCarbayo)

## define params
# Use threshold 3 to select contributing genes.
# Results of enrichment analysis will be written in path 'resPath(params)/GOstatsEnrichAnalysis'
params <- buildMineICAParams(resPath="~/resMineICACarbayo/", selCutoff=3)

## Annotation package for IcaSetCarbayo is hgu133a.db.
# check annotation package
annotation(icaSetCarbayo)

## Define universe, i.e the set of EntrezGene IDs mapping to the feature IDs of the IcaSet object.
universe <- as.character(na.omit(unique(unlist(AnnotationDbi::mget(featureNames(icaSetCarbayo),
                         hgu133aENTREZID, ifnotfound = NA)))))

## Apply enrichement analysis (of the contributing genes) to the first components using gene sets from KEGG.
# Since an annotation package is available, we don't need to fill arg 'entrez2symbol'.
# run the actual enrichment analysis
hypergeoAn(icaSet=icaSetCarbayo[,,1], params=params, db="GO",onto="BP", universe=universe)
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{runEnrich}}, \code{\link[xtable]{xtable}},
  \code{\link[biomaRt]{useMart}},
  \code{\link[Category]{hyperGTest}},
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}},
  \code{\link{mergeGostatsResults}}
}

