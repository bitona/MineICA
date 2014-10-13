\name{mergeGostatsResults}
\alias{mergeGostatsResults}
\title{Merge enrichment results obtained for different databases into one file per component.}
\usage{
  mergeGostatsResults(resPath, GOstatsPath,
    rdata = "hgres", cutoffs = NULL, hgCutoff = 0.01,
    cond = TRUE, pathGenes)
}
\arguments{
  \item{resPath}{The global path where results of ICA
  analysis are written}

  \item{GOstatsPath}{The path within argument
  \code{resPath} where files will be written}

  \item{rdata}{The name of the rdata file containing the
  enrichment analysis of all components}

  \item{cutoffs}{The threshold(s) used to select genes used
  in enrichment analysis}

  \item{hgCutoff}{The p-value threshold}

  \item{cond}{A logical indicating whether the calculation
  has been conditioned on the GO structure, see
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.}

  \item{pathGenes}{The path where HTML files containing
  gene projections for each component are located}
}
\value{
  NULL
}
\description{
  This function is internal and called by function
  \code{runEnrich}. It merges enrichment results obtained
  with either KEGG, GO, or both databases into one file.
}
\details{
  This function writes an HTML file per component,
  containing the outputs of the enrichment tests computed
  through the function \code{\link[Category]{hyperGTest}}.
  The results of the enrichment tests are loaded from .rda
  files located in
  \code{resPath(icaSet)}/GOstatsEnrichAnalysis/byDb/'db-name'/('ontology-name'/).
  The results obtained for the different
  databases/ontologies are then merged into an array for
  each component, this array is written as an HTML file in
  the directory
  \code{resPath(icaSet)}/\code{zvalCutoff(params)}. The
  arguments \code{hgCutoff} and \code{cond} have to be
  provided because they will be used in the file names of
  the resulting files.

  This function makes several important assumptions: only
  databases GO and KEGG have been tested, p-values are not
  available for gene sets that have not been selected as
  significant.

  The outputs of \code{\link[Category]{hyperGTest}} that
  are given in each table are: \describe{ \item{DB, ID,
  Term}{The database, the gene set ID, and the gene set
  name,} \item{P-value}{probability of observing the number
  of genes annotated for the gene set among the selected
  gene list, knowing the total number of annotated genes
  among the universe}, \item{Expected counts}{expected
  number of genes in the selected gene list to be found at
  each tested category term/gene set,} \item{Odds
  ratio}{odds ratio for each category term tested which is
  an indicator of the level of enrichment of genes within
  the list as against the universe,} \item{Counts}{number
  of genes in the selected gene list which are annotated
  for the gene set,} \item{Size}{number of genes from the
  universe annotated for the gene set.}}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link[xtable]{xtable}},
  \code{\link[biomaRt]{useMart}},
  \code{\link[Category]{hyperGTest}},
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}},
  \code{\link{hypergeoAn}},
  \code{\link{mergeGostatsResults}}
}
\keyword{internal}

