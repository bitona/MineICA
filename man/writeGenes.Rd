\name{writeGenes}
\alias{writeGenes}
\title{Description of features using package \code{biomaRt}.}
\usage{
  writeGenes(data, filename = NULL,
    mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"),
    typeId = "hgnc_symbol", typeRetrieved = NULL,
    sortBy = NULL, sortAbs = TRUE, colAnnot = NULL,
    decreasing = TRUE, highlight = NULL, caption = "")
}
\arguments{
  \item{data}{Either a data.frame whose rownames or one of
  its columns contain the IDs to be annotated, or a vector
  of IDs.}

  \item{filename}{The name of the HTML file where gene
  annotations are written.}

  \item{mart}{Output of function \code{useMart} from
  package \code{biomaRt}.}

  \item{typeId}{The type of IDs available in \code{data},
  in the biomaRt way (type \code{listFilters(mart)} to
  choose one).}

  \item{typeRetrieved}{The descriptors uses to annotate the
  features of \code{data} (type \code{listAttributes(mart)}
  to choose one or several).}

  \item{sortBy}{Name of a column of \code{data} used to
  order the output.}

  \item{sortAbs}{If TRUE absolute value of column
  \code{sortBy} is used to order the output.}

  \item{colAnnot}{The column containing the IDs to be
  annotated, if NULL or missing and argument \code{data} is
  a data.frame, then rownames of \code{data} must contain
  the IDs.}

  \item{decreasing}{If TRUE, the output is sorted by
  decreasing values of the \code{sortBy} column}

  \item{highlight}{IDs to be displayed in colour red in the
  returned table}

  \item{caption}{A title for the HTML table}
}
\value{
  This function returns a data.frame which contains
  annotations of the input data.
}
\description{
  This function annotates IDs (typically gene IDs) provided
  by the user and returns an html file with their
  description.
}
\details{
  \code{"hgnc_symbol", "ensembl_gene_id", "description",
  "chromosome_name", "start_position", "end_position",
  "band"}, and \code{"strand"}, are automatically added to
  the list of fields available in argument
  \code{typeRetrieved} queried on biomaRt. The web-links to
  www.genecards.org and www.proteinatlas.org are
  automatically added in the columns of the output
  respectively corresponding to \code{hgnc_symbol} and
  \code{ensembl_gene_id}.
}
\examples{
if (interactive()) {
## define the database to be used
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

### Describe:
## a set of hgnc symbols with default descriptions (typeRetrieved=NULL)
genes <- c("TOP2A","E2F3","E2F1","CDK1","CDC20","MKI67")
writeGenes(data=genes, filename="foo", mart=mart, typeId = "hgnc_symbol")

## a data.frame indexed by hngc symbols, sort output according to column "values", add a title to the HTML output
datagenes <- data.frame(values=rnorm(6),row.names = genes)
writeGenes(data=datagenes, filename="foo", sortBy = "values", caption = "Description of some proliferation genes.")

## a set of Entrez Gene IDs with default descriptions
genes <- c("7153","1871","1869","983","991","4288")
writeGenes(data=genes, filename="foo", mart=mart, typeId = "entrezgene")
}
\dontrun{
## add the GO category the genes belong to
## search in listAttributes(mart)[,1] which filter correspond to the Gene Ontology -> "go_id"
writeGenes(data=genes, filename="foo", mart=mart, typeId = "entrezgene", typeRetrieved = "go_id")
}
}
\author{
  Anne Biton
}
\seealso{
  \code{\link[biomaRt]{getBM}},
  \code{\link[biomaRt]{listFilters}},
  \code{\link[biomaRt]{listAttributes}},
  \code{\link[biomaRt]{useMart}}
}

