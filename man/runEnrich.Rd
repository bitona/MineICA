\name{runEnrich}
\alias{runEnrich}
\title{Enrichment analysis through \link{GOstats}}
\usage{
  runEnrich(icaSet, params, dbs = c("KEGG", "GO"),
    ontos = c("BP", "CC", "MF"), cond = TRUE,
    hgCutoff = params["pvalCutoff"])
}
\arguments{
  \item{icaSet}{An object of class
  \code{\link[MineICA:class-IcaSet]{IcaSet}}}

  \item{params}{An object of class
  \code{\link[MineICA:class-MineICAParams]{MineICAParams}}
  providing the parameters of the analysis}

  \item{dbs}{The database to use, default is
  \code{c("GO","KEGG")}}

  \item{ontos}{A string specifying the GO ontology to use.
  Must be one of \code{"BP"}, \code{"CC"}, or \code{"MF"},
  see \code{\link[Category]{GOHyperGParams-class}}. Only
  used when argument \code{dbs} includes \code{"GO"}.}

  \item{cond}{A logical indicating whether the calculation
  should condition on the GO structure, see
  \code{\link[Category]{GOHyperGParams-class}}. Only used
  when argument \code{dbs} includes \code{"GO"}.}

  \item{hgCutoff}{The threshold p-value for statistical
  significance, default is \code{pvalCutoff(params)}}
}
\value{
  NULL
}
\description{
  This function tests the enrichment of the components of
  an \code{\link{IcaSet}} object using package
  \code{GOstats} through function \code{hyperGTest}.
}
\details{


  An annotation package should be available in
  \code{annotation(icaSet)} to provide the contents of the
  gene sets. If none corresponds to the technology you deal
  with, please choose the org.*.eg.db package according to
  the organism (for example org.Hs.eg.db for Homo sapiens).
  By default, if \code{annotation(icaSet)} is empty and
  organism is one of
  \code{c("Human","HomoSapiens","Mouse","Mus Musculus")},
  then either \code{org.Hs.eg.db} or \code{org.Mm.eg.db} is
  used.

  Use of \code{GOstats} requires the input IDs to be Entrez
  Gene, this function will therefore annotate either the
  feature names or the gene names into Entrez Gene ID using
  either the annotation package (\code{annotation(icaSet)})
  or \code{biomaRt}.

  Three types of enrichment tests are computed for each
  component: the threshold is first used to select gene
  based on their absolute projections, then positive and
  negative projections are treated individually.

  For each database \code{db} (each ontology if \code{db}
  is "GO"), this function writes an HTML file containing
  the outputs of the enrichment tests computed through the
  function \code{\link[Category]{hyperGTest}}. The
  corresponding files are located in
  \code{resPath(icaSet)}/GOstatsEnrichAnalysis/byDb/. The
  results obtained for each database/ontology are then
  merged into an array for each component, this array is
  written as an HTML file in the directory
  \code{resPath(icaSet)}/GOstatsEnrichmentAnalysis/ (this
  directory is first deleted if it already exists). This
  file is the one the user should look at.

  The outputs of \code{\link[Category]{hyperGTest}} that
  are given in each table are: \describe{ \item{DB, ID,
  Term:}{the database, the gene set ID, and the gene Set
  name} \item{P-value:}{probability of observing the number
  of genes annotated for the gene set among the selected
  gene list, knowing the total number of annotated genes
  among the universe}, \item{Expected counts:}{expected
  number of genes in the selected gene list to be found at
  each tested category term/gene set,} \item{Odds
  ratio:}{odds ratio for each category term tested which is
  an indicator of the level of enrichment of genes within
  the list as against the universe,} \item{Counts:}{number
  of genes in the selected gene list that are annotated for
  the gene set,} \item{Size:}{number of genes from the
  universe annotated for the gene set.}}
}
\examples{
\dontrun{
# Load examples of IcaSet object
data(icaSetCarbayo)

## Define parameters
# Use threshold 3 to select contributing genes on which enrichment analysis will be applied
# Results of enrichment analysis will be written in path 'resPath(params)/GOstatsEnrichAnalysis'
params <- buildMineICAParams(resPath="carbayo/", selCutoff=3)

## Run enrichment analysis on the first two components contained in the icaSet object 'icaSetCarbayo'
runEnrich(params=params,icaSet=icaSetCarbayo[,,1:2],dbs="GO", ontos="BP")

}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link{buildIcaSet}},
  \code{\link[biomaRt]{useMart}},
  \code{\link[Category]{hyperGTest}},
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}},
  \code{\link{hypergeoAn}},
  \code{\link{mergeGostatsResults}}
}

