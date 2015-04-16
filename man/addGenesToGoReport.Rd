\name{addGenesToGoReport}
\alias{addGenesToGoReport}
\title{Add Symbol IDs to hyperGTest results}
\usage{
  addGenesToGoReport(hgOver, universe,
    db = c("GO", "KEGG"), onto = c("CC", "MF", "BP"),
    annotation = NULL, entrez2symbol = NULL)
}
\arguments{
  \item{hgOver}{Output of function
  \code{\link[Category]{hyperGTest}}}

  \item{universe}{A vector including all IDs on which
  enrichment analysis was applied}

  \item{db}{The database to use, default is c("GO","KEGG")}

  \item{onto}{A string specifying the GO ontology to use.
  Must be one of "BP", "CC", or "MF", see
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.
  Only used when argument \code{db} is "GO".}

  \item{annotation}{An annotation package}

  \item{entrez2symbol}{A vector indexed by Entrez Gene ID
  and filled with the corresponding Gene Symbols}
}
\value{
  A data.frame containing the summary of the output of
  function hyperGTest (\code{summary(hgOver)}) with an
  additional column providing the gene Symbols included in
  the significant gene sets.
}
\description{
  Add gene Symbols contained in gene sets selected as
  significant by \code{\link[Category]{hyperGTest}}
  function
}
\details{
  This function takes as inputs the outputs of
  \code{\link[Category]{hyperGTest}} which takes Entrez
  Gene IDs as inputs to perform the enrichment analysis.
  The goal of this function is to select the Entrez Gene
  IDs responsible for the significant enrichment of a given
  gene set and annotate them in to gene Symbol IDs. When
  the annotation package \code{annotation} was used to map
  feature IDs to Entrez Gene ID, it can also be used here
  to map Entrez and Symbol IDs. If the annotation package
  was not used, but the Entrez Gene IDs were directly
  provided to the hyperGtest function, \code{annotation} is
  expected to be NULL and \code{entrez2symbol} must be
  specified.

  This function returns the outputs of function
  \code{\link[Category]{hyperGTest}} which contain:
  \describe{ \item{DB, ID, Term}{The database, the gene set
  ID, and the gene Set name,} \item{P-value}{probability of
  observing the number of genes annotated for the gene set
  among the selected gene list, knowing the total number of
  annotated genes among the universe}, \item{Expected
  counts}{expected number of genes in the selected gene
  list to be found at each tested category term/gene set,}
  \item{Odds ratio}{odds ratio for each category term
  tested which is an indicator of the level of enrichment
  of genes within the list as against the universe,}
  \item{Counts}{number of genes in the selected gene list
  which are annotated for the gene set,} \item{Size}{number
  of genes from the universe annotated for the gene set.}}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link[Category]{hyperGTest}},
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}
}
\keyword{internal}

