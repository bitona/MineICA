\name{writeGostatsHtmltable}
\alias{writeGostatsHtmltable}
\title{Writes enrichment results in a HTML file}
\usage{
  writeGostatsHtmltable(d, label, side = "both", db, file,
    cutoff = 3)
}
\arguments{
  \item{d}{A data.frame describing enrichment results,
  output of function \code{\link[Category]{hyperGTest}}}

  \item{label}{The label of the data the results originate
  from}

  \item{side}{The side of the component used for enrichment
  analysis}

  \item{db}{The database used ("GO" or "KEGG")}

  \item{file}{File name for output}

  \item{cutoff}{The threshold used to select the genes used
  to run the enrichment analysis}
}
\value{
  NULL
}
\description{
  This function takes as input in argument \code{d} the
  output of function \code{\link{addGenesToGoReport}} whose
  goal is to add genes included in gene sets detected as
  significantly enriched by
  \code{\link[Category]{hyperGTest}} function. It writes
  the enrichment results in an HTML file which redirects
  each gene set ID to its web-description and each gene to
  its Gene Card web-page.
}
\examples{
hgOver <- structure(list(GOBPID = c("GO:0003012", "GO:0030049"),
                    Pvalue = c(1.70848789161935e-10, 6.62508415367712e-05),
                    OddsRatio = c(22.1043956043956, 26.4190476190476),
                    ExpCount = c(1.19549929676512, 0.246132208157525),
                    Count = c(12L, 4L), Size = c(68L, 14L),
                    Term = c("muscle system process", "muscle filament sliding"),
                    In_geneSymbols = c("ACTA2,ACTC1,ACTG2,CASQ2,CNN1,DES,MYH3,MYLK,PTGS1,TPM2,MYL9,LMOD1","ACTC1,DES,MYH3,TPM2")),
                    .Names = c("GOBPID", "Pvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term", "In_geneSymbols"),
                    class = "data.frame", row.names=1:2)

MineICA:::writeGostatsHtmltable(d=hgOver, label="Example of enrichment analysis", db="KEGG",
                                file="outputHyper_example.htm")
}
\author{
  Anne Biton
}
\seealso{
  \code{\link[xtable]{xtable}},
  \code{\link{addGenesToGoReport}},
  \code{\link[Category]{hyperGTest}}
}
\keyword{internal}

