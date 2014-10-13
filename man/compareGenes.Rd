\name{compareGenes}
\alias{compareGenes}
\title{Union and intersection of contributing genes}
\usage{
  compareGenes(keepCompByIcaSet, icaSets, lab, cutoff = 0,
    type = c("union", "intersection"), annotate = TRUE,
    file,
    mart = useMart("ensembl", "hsapiens_gene_ensembl"))
}
\arguments{
  \item{icaSets}{List of \code{IcaSet} objects, the
  \code{geneNames} of the \code{IcaSet} objects must be
  from the same type (e.g, gene Symbols).}

  \item{keepCompByIcaSet}{Indices of the components to be
  considered in each \code{IcaSet}.}

  \item{lab}{The names of the icaSets (e.g the names of the
  datasets they originate from).}

  \item{cutoff}{The cutoff (on the absolute centered and
  scaled projections) above which the genes have to be
  considered.}

  \item{type}{\code{"intersection"} to restrict the list of
  genes to the ones that are common between all datasets,
  or \code{"union"} to consider all the union of genes
  available across the datasets.}

  \item{annotate}{If TRUE (default) the genes are annotated
  using function \code{writeGenes}.}

  \item{file}{The HTML file name where the genes and their
  annotations are written, default is
  \code{type}Genes_lab1-i_lab2-j_... where i and j are the
  component indices contained in \code{keepCompByIcaSet}.}

  \item{mart}{The mart object (database and dataset) used
  for annotation, see function \code{useMart} of package
  \code{biomaRt}.}
}
\value{
  A data.frame containing \describe{
  \item{\code{typeID(icaSets[[1]])['geneID_biomart']}:}{the
  gene IDs,} \item{median_rank}{the median of the ranks of
  each gene across the \code{IcaSet} objects,}
  \item{analyses}{the labels of the \code{IcaSet} objects
  in which each gene is above the given \code{cutoff}}
  \item{min_rank}{the minimum of the ranks of each gene
  across the \code{IcaSet} objects,} \item{ranks}{the ranks
  of each gene in each \code{IcaSet} where it is
  available,} \item{scaled_proj}{the centered and reduced
  projection of each gene in each \code{IcaSet} where it is
  available.}}
}
\description{
  Compute and annotate the intersection or union between
  contributiong genes of components originating from
  different IcaSet objects.
}
\examples{
\dontrun{
data(icaSetCarbayo)
mart <- useMart("ensembl", "hsapiens_gene_ensembl")

## comparison of two components
## here the components come from the same IcaSet for convenience
## but they must come from different IcaSet in practice.
compareGenes(keepCompByIcaSet = c(9,4), icaSets = list(icaSetCarbayo, icaSetCarbayo),
             lab=c("Carbayo", "Carbayo2"), cutoff=3, type="union",  mart=mart)

}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link[MineICA]{writeGenes}}
}

