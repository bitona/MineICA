\name{doEnrichment}
\alias{doEnrichment}
\title{Runs enrichment analysis of contributing genes}
\usage{
  doEnrichment(compSel, chip, onto, hgCutoff, cond,
    universe, path, db, pack.annot.EID, Slist, it, cutoff,
    entrez2symbol)
}
\arguments{
  \item{compSel}{A list containing three elements
  \describe{ \item{compSel}{the projection values of
  contributing genes that were selected based on their
  absolute projection} \item{compSel_neg}{the projection
  values of contributing genes that have negative
  projections} \item{compSel_pos}{the projection values of
  contributing genes that have positive projections}}}

  \item{chip}{The annotation package}

  \item{onto}{A string specifying the GO ontology to use.
  Must be one of 'BP', 'CC', or 'MF', see
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.
  Only used when argument \code{db} is 'GO'.}

  \item{hgCutoff}{The p-value threshold}

  \item{cond}{A logical indicating whether the calculation
  should conditioned on the GO structure, see
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.}

  \item{universe}{The universe for the hypergeometric
  tests, see
  \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.}

  \item{path}{The path where results will be saved}

  \item{db}{The used database to use ('GO' or 'KEGG')}

  \item{pack.annot.EID}{The name of the environment of the
  annotation package containing the annotation for Entrez
  Gene.}

  \item{Slist}{The list of gene projections across all
  components}

  \item{it}{The index of the component}

  \item{cutoff}{The threshold applied on the gene
  projections, used to select the contributing genes}

  \item{entrez2symbol}{A vector of all gene Symbols
  involved in the analysis indexed by their Entrez Gene
  IDs. It is only used when \code{annotation(params)} is
  empty, and allows to associate gene sets to Symbols.}
}
\value{
  Object of class \code{GOHyperGResult-class}
}
\description{
  doEnrichment This internal function is called by
  \code{hypergeoAn} and runs hypergeometric tests through
  function \code{hyperGTest} to associate the contributing
  genes of a component to gene sets.
}
\author{
  Anne Biton
}
\keyword{internal}

