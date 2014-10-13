\name{annotFeaturesWithBiomaRt}
\alias{annotFeaturesWithBiomaRt}
\title{Annotation of features using \code{biomaRt}}
\usage{
  annotFeaturesWithBiomaRt(features, featureId, geneId,
    mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))
}
\arguments{
  \item{features}{Feature IDs to be annotated}

  \item{featureId}{The type of the feature IDs, in the
  \code{biomaRt} way (type \code{listFilters(mart)} to
  choose one)}

  \item{geneId}{The type of the gene IDs, in the
  \code{biomaRt} way (type \code{listAttributes(mart)} to
  choose one)}

  \item{mart}{The mart object (database and dataset) used
  for annotation, see function \code{useMart} of package
  \code{biomaRt}}
}
\value{
  A vector of gene IDs indexed by the feature IDs.
}
\description{
  This function annotates a set of features using
  \code{biomaRt}
}
\examples{
if (interactive()) {
# define the database to be queried by biomaRt
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# annotate a set of HG-U133a probe sets IDs into Gene Symbols
annotFeaturesWithBiomaRt(features = c("1007_s_at", "1053_at", "117_at", "121_at", "1255_g_at"),
featureId="affy_hg_u133a", geneId="hgnc_symbol", mart=mart)

# annotate a set of Ensembl Gene IDs into Gene Symbols
annotFeaturesWithBiomaRt(features = c("ENSG00000101412", "ENSG00000112242",
                                      "ENSG00000148773", "ENSG00000131747", "ENSG00000170312",
                                      "ENSG00000117399"), featureId="ensembl_gene_id", geneId="hgnc_symbol", mart=mart)
}
}
\author{
  Anne Biton
}

