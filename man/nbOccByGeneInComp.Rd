\name{nbOccByGeneInComp}
\alias{nbOccByGeneInComp}
\title{nbOccByGeneInComp}
\description{For each feature/gene, this function returns the indices of the components they contribute to. }
\usage{nbOccByGeneInComp(Slist, cutoff, sel)}
\arguments{
  \item{Slist}{A list whose each element contains projection values of
    features/genes on a component.}
  \item{cutoff}{A threshold to be used to define a gene as contributor}
  \item{sel}{A list whose each element contains projection values of
    contributing features/genes on a component (the difference with arg
   \code{Slist} is that \code{sel} is already restricted to the contributing genes).}
}

\value{This function returns a list which gives for each feature/gene the indices of the components it contributes to.}

\examples{
 c1 <- rnorm(100); names(c1) <- paste("g",100:199,sep="")
 c2 <- rnorm(100); names(c2) <- paste("g",1:99,sep="")
 MineICA:::nbOccByGeneInComp(Slist=list(c1,c2), cutoff= 0.5)

 }
\author{Anne Biton}





