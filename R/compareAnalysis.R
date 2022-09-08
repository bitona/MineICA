##' Compare \code{\link{IcaSet}} objects by computing the correlation between either
##' projection values of common features or genes, or contributions of common
##' samples.
##' 
##' The user must carefully choose the object on which the correlation will be
##' computed.
##' If \code{level='samples'}, the correlations are based on the mixing
##' matrices of the ICA decompositions (of dimension samples x components). \code{'A'}
##' will be typically chosen when the ICA decompositions were computed on the
##' same dataset, or on datasets that include the same samples.
##' If \code{level='features'} is
##' chosen, the correlation is calculated between the source matrices (of dimension
##' features x components) of the ICA decompositions. \code{'S'} will be typically
##' used when the ICA decompositions share common features (e.g same microarrays).
##' If \code{level='genes'}, the correlations are calculated
##' on the attributes \code{'SByGene'} which store the 
##' projections of the annotated features. \code{'SByGene'} will be typically chosen
##' when ICA were computed on datasets from different technologies, for which
##' comparison is possible only after annotation into a common ID, like genes.
##' 
##' \code{cutoff_zval} is only used when \code{level} is one of \code{c('genes','features')}, in
##' order to restrict the correlation to the contributing features or genes.
##' 
##' When \code{cutoff_zval} is specified, for each pair of components, genes or features that are
##' included in the circle of center 0 and radius \code{cutoff_zval} are excluded from
##' the computation of the correlation.
##'
##' It must be taken into account by the user that if
##' \code{cutoff_zval} is different from \code{NULL} or \code{0}, the computation will be much
##' slowler since each pair of component is treated individually.
##'
##' @title Comparison of IcaSet objects using correlation 
##' @param icaSets list of IcaSet objects, e.g results of ICA decompositions
##' obtained on several datasets.
##' @param labAn vector of names for each icaSet, e.g the the names of the
##' datasets on which were calculated the decompositions.
##' @param type.corr Type of correlation to compute, either
##' \code{'pearson'} or \code{'spearman'}.
##' @param cutoff_zval either NULL or 0 (default) if all genes are used to
##' compute the correlation between the components, or a threshold to compute
##' the correlation on the genes that have at least a scaled projection higher
##' than cutoff_zval. Will be used only when correlations are calculated on S or
##' SByGene.
##' @param level Data level of the \code{IcaSet} objects on which is applied the correlation.
##' It must correspond to a feature shared by the IcaSet objects:
##' \code{'samples'} if they were applied to common samples (correlations are computed between matrix \code{A}), \code{'features'} if they were applied to
##' common features (correlations are computed between matrix \code{S}), \code{'genes'} if they share gene IDs after
##' annotation into genes (correlations are computed between matrix \code{SByGene}).
##' @return A list whose length equals the number of pairs of \code{IcaSet} and whose elements
##' are outputs of function \code{\link{cor2An}}. 
##' @export
##' @examples
##'
##' dat1 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat1) <- paste("g", 1:1000, sep="")
##' colnames(dat1) <- paste("s", 1:10, sep="")
##' dat2 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat2) <- paste("g", 1:1000, sep="")
##' colnames(dat2) <- paste("s", 1:10, sep="")
##' 
##' ## run ICA
##' resJade1 <- runICA(X=dat1, nbComp=3, method = "JADE")
##' resJade2 <- runICA(X=dat2, nbComp=3, method = "JADE")
##' 
##' ## build params
##' params <- buildMineICAParams(resPath="toy/")
##'
##' ## build IcaSet object
##' icaSettoy1 <- buildIcaSet(params=params, A=data.frame(resJade1$A), S=data.frame(resJade1$S),
##'                           dat=dat1, alreadyAnnot=TRUE)$icaSet
##' icaSettoy2 <- buildIcaSet(params=params, A=data.frame(resJade2$A), S=data.frame(resJade2$S),
##'                           dat=dat2, alreadyAnnot=TRUE)$icaSet
##' 
##' listPairCor <- compareAn(icaSets=list(icaSettoy1,icaSettoy2), labAn=c("toy1","toy2"), 
##'                          type.corr="pearson", level="genes", cutoff_zval=0)
##'
##' 
##' \dontrun{
##' #### Comparison of 2 ICA decompositions obtained on 2 different gene expression datasets.
##' ## load the two datasets
##' library(breastCancerMAINZ)
##' library(breastCancerVDX)
##' data(mainz)
##' data(vdx)
##' 
##' ## Define a function used to build two examples of IcaSet objects
##' treat <- function(es, annot="hgu133a.db") {
##'    es <- selectFeatures_IQR(es,10000)
##'    exprs(es) <- t(apply(exprs(es),1,scale,scale=FALSE))
##'    colnames(exprs(es)) <- sampleNames(es)
##'    resJade <- runICA(X=exprs(es), nbComp=10, method = "JADE", maxit=10000)
##'    resBuild <- buildIcaSet(params=buildMineICAParams(), A=data.frame(resJade$A), S=data.frame(resJade$S),
##'                         dat=exprs(es), pData=pData(es), refSamples=character(0),
##'                         annotation=annot, typeID= typeIDmainz,
##'                         chipManu = "affymetrix", mart=mart)
##'    icaSet <- resBuild$icaSet
##' }
##' ## Build the two IcaSet objects
##' icaSetMainz <- treat(mainz)
##' icaSetVdx <- treat(vdx)
##' 
##' ## The pearson correlation is used as a measure of association between the gene projections
##' # on the different components (type.corr="pearson").
##' listPairCor <- compareAn(icaSets=list(icaSetMainz,icaSetVdx),
##' labAn=c("Mainz","Vdx"), type.corr="pearson", level="genes", cutoff_zval=0)
##' 
##' ## Same thing but adding a selection of genes on which the correlation between two components is computed:
##' # when considering pairs of components, only projections whose scaled values are not located within
##' # the circle of radius 1 are used to compute the correlation (cutoff_zval=1).
##' listPairCor <-  compareAn(icaSets=list(icaSetMainz,icaSetVdx),
##' labAn=c("Mainz","Vdx"), type.corr="pearson", cutoff_zval=1, level="genes")
##' }
##' @export
##' @author Anne Biton
##' @seealso \code{\link{cor2An}}

compareAn <- function (icaSets,
                       labAn,
                       type.corr = c("pearson","spearman"),
                       cutoff_zval=0,
                       level = c("samples","features","genes")
                       ) {
    
    if (missing(labAn))
        if (is.null(names(icaSets)))
            labAn <- c(1:nbAn)
        else
            labAn <- names(icaSets)
    
    nbAn <- length(icaSets)

    if (length(labAn) != nbAn)
        stop("labAn must have the same length than the list icaSets")
          
    level <- match.arg(tolower(level), choices = c("features","genes","samples"))
    switch(level,
           features={ object="S"},
           genes={object="SByGene"},
           samples={object="A"}
           )

    
 type.corr <- match.arg(type.corr)
  nbCompByAn <- lapply(icaSets, function(x) ncol(A(x)))
  
      
  matlist <-
      lapply(icaSets,
             function(icaSet, object) {
                 data.frame(icaSet[object],check.names = FALSE, stringsAsFactors = FALSE)
             }
             , object = object
             )

  
  
  allpairs <- combn(x=1:nbAn,m = 2, simplify = FALSE)
  pairnames <- unlist(lapply(allpairs, function(x,labAn) paste(labAn[x[1]],labAn[x[2]],sep="-"), labAn = labAn))

  listCorPair <- 
          lapply(allpairs,
             function(pa, matlist, labAn, cutoff_zval, type.corr) {
                 mat1 <- matlist[[pa[1]]]
                 mat2 <- matlist[[pa[2]]]
                 resCor <- cor2An(mat1, mat2, lab = c(labAn[pa[1]],labAn[pa[2]]), type.corr = type.corr, cutoff_zval = cutoff_zval)
                 
                 
                 
             }
             , matlist = matlist
             , labAn = labAn
             , type.corr = type.corr
             , cutoff_zval = if (missing(cutoff_zval)) NA else cutoff_zval
            )
  names(listCorPair) <- pairnames
 
  return(listCorPair)
}



##' This function measures the correlation between two matrices containing the results of
##' two decompositions.
##' 
##' Before computing the correlations, the components are scaled and restricted
##' to common row names.
##' 
##' It must be taken into account by the user that if \code{cutoff_zval} is different
##' from NULL or zero, the computation will be slowler since each pair of
##' component is treated individually.
##'
##' When \code{cutoff_zval} is specified, for each
##' pair of components, genes that are included in the circle of center 0 and
##' radius \code{cutoff_zval} are excluded from the computation of the correlation
##' between the gene projection of the two components.
##' @title Correlation between two matrices 
##' @param mat1 matrix of dimension features/genes x number of components, e.g
##' the results of an ICA decomposition
##' @param mat2 matrix of dimension features/genes x number of components, e.g
##' the results of an ICA decomposition
##' @param lab The vector of labels for mat1 and mat2, e.g the the names of the
##' two datasets on which were calculated the two decompositions
##' @param type.corr Type of correlation, either \code{'pearson'} or
##' \code{'spearman'}
##' @param cutoff_zval cutoff_zval: 0 (default) if all genes are
##' used to compute the correlation between the components, or a threshold to
##' compute the correlation on the genes that have at least a scaled projection
##' higher than cutoff_zval.
##' @return This function returns a list consisting of:
##' \item{cor}{a matrix of dimensions '(nbcomp1+nbcomp2) x (nbcomp1*nbcomp2)', 
##' containing the correlation values between each pair of components,}
##' \item{pval}{ a matrix of dimension '(nbcomp1+nbcomp2) x (nbcomp1*nbcomp2)', 
##' containing the p-value of the correlation tests for each
##' pair of components,}
##' \item{inter}{ the intersection between the features/genes of \code{mat1}
##' and \code{mat2},}
##' \item{labAn}{ the labels of the compared matrices.}
##' @author Anne Biton
##' @export
##' @examples 
##' cor2An(mat1=matrix(rnorm(10000),nrow=1000,ncol=10), mat2=matrix(rnorm(10000),nrow=1000,ncol=10),
##'        lab=c("An1","An2"), type.corr="pearson")
##' 
##' @seealso \code{rcorr}, \code{cor.test}, \code{\link{compareAn}}

cor2An <- function (mat1,
                    mat2,
                    lab,
                    type.corr = c("pearson","spearman"),
                    cutoff_zval = 0
                    ) {
    
    if (is.null(cutoff_zval))
        cutoff_zval <- 0

    if (cutoff_zval < 0)
        stop("'cutoff_zval' must be positive")

    comp1 <- comp2 <- NULL

    mat1 <- data.frame(scale(mat1), check.names = FALSE, stringsAsFactors = FALSE)
    mat2 <- data.frame(scale(mat2), check.names = FALSE, stringsAsFactors = FALSE)
    nbcomp1 <- ncol(mat1)
    nbcomp2 <- ncol(mat2)
    
    gg <- intersect(rownames(mat1),rownames(mat2))

    if (length(gg)<10)
        warning(paste("Less than 10 elements are common between the icaSet objects ", lab[1], " and ", lab[2], ", the comparison of these two decompositions should be reconsidered.", sep = ""))
    
    mat1 <- mat1[gg,]
    mat2 <- mat2[gg,]

    if (missing(lab))
        lab <- paste("an",c(1:2),sep = "")

    if (cutoff_zval == 0) {
        r <- signif(rcorr(as.matrix(mat1), as.matrix(mat2), type = type.corr)$r, digits = 5)
        dimnames(r) <- list(c(paste("comp",c(1:nbcomp1), "_", lab[1], sep = ""),  paste("comp",c(1:nbcomp2), "_", lab[2], sep = "")),c(paste("comp",c(1:nbcomp1), "_", lab[1], sep = ""),  paste("comp",c(1:nbcomp2), "_", lab[2], sep = "")))
    
        rpval <- NULL	

        rpval <- rcorr(as.matrix(mat1),as.matrix(mat2),type=type.corr)$P
        dimnames(rpval) <- list(c(paste("comp",c(1:nbcomp1), "_", lab[1], sep = ""),  paste("comp",c(1:nbcomp2), "_", lab[2], sep = "")),c(paste("comp",c(1:nbcomp1), "_", lab[1], sep = ""),  paste("comp",c(1:nbcomp2), "_", lab[2], sep = "")))

    }
    else {

            ## transform matrix into lists
            l1 <- lapply(as.list(mat1),
                        function(x,n) {
                              names(x) <- n
                              return(x)
                        }
                        , n = rownames(mat1) 
                  )
            l2 <- lapply(as.list(mat2),
                        function(x,n) {
                            names(x) <- n
                            return(x)
                        }
                        , n = rownames(mat2) 
                  )
            
            rpval <- matrix(ncol = nbcomp2+nbcomp1 , nrow =nbcomp2+nbcomp1)
            r <- matrix(ncol = nbcomp2+nbcomp1 , nrow =nbcomp2+nbcomp1)

            for (i in 1:nbcomp1) {
                   reslist <- 
                sapply(1:nbcomp2,
                   function(j, i, l1, l2, cutoff_zval, gg) {
                       dc <- data.frame(comp1=l1[[i]],comp2=l2[[j]])
                       rownames(dc) <- gg
                       g2del <- rownames(subset(dc, comp1^2 + comp2^2<= cutoff_zval   ))
                       inter <- gg[!(gg %in% g2del)]
                       res <- cor.test(as.matrix(mat1[inter,i]),as.matrix(mat2[inter,j]), method = type.corr)
                       return(list(cor=res$estimate,pval=res$p.value))
                    }
                   , i = i
                   , l1 = l1
                   , l2 = l2
                   , cutoff_zval = cutoff_zval
                   , gg = gg)
               
               reslist <- as.data.frame(reslist)
               r[i,(nbcomp1+1):(nbcomp1+nbcomp2)] <- r[(nbcomp1+1):(nbcomp1+nbcomp2),i] <- unlist(reslist["cor",])
               rpval[i,(nbcomp1+1):(nbcomp1+nbcomp2)] <- rpval[(nbcomp1+1):(nbcomp1+nbcomp2),i] <- unlist(reslist["pval",])

                for (j in 1:nbcomp1) {
                    dc <- data.frame(comp1=l1[[i]],comp2=l1[[j]])
                    rownames(dc) <- gg
                    g2del <- rownames(subset(dc, comp1^2 + comp2^2<= cutoff_zval   ))
                    inter <- names(l1[[1]])[!(names(l1[[1]]) %in% g2del)]
                    res <- cor.test(as.matrix(mat1[inter,i]),as.matrix(mat1[inter,j]), method = type.corr)
                    rpval[i,j] <- rpval[j,i] <- res$p.value
                    r[i,j] <- r[j,i] <- res$estimate
                }
            }

            for (i in 1:nbcomp2) {
                for (j in 1:nbcomp2) {
                    dc <- data.frame(comp1=l2[[i]],comp2=l2[[j]])
                    rownames(dc) <- gg
                    g2del <- rownames(subset(dc, comp1^2 + comp2^2<= cutoff_zval   ))
                    inter <- names(l2[[1]])[!(names(l2[[1]]) %in% g2del)]
                    res <- cor.test(as.matrix(mat2[inter,i]),as.matrix(mat2[inter,j]), method = type.corr)
                    rpval[nbcomp1+i,nbcomp1+j] <- rpval[nbcomp1+j,nbcomp1+i] <- res$p.value
                    r[nbcomp1+i,nbcomp1+j] <- r[nbcomp1+j,nbcomp1+i] <- res$estimate
                    
                }
            }
            
            dimnames(rpval) <- list(c(paste("comp",c(1:nbcomp1), "_", lab[1], sep = ""),  paste("comp",c(1:nbcomp2), "_", lab[2], sep = "")),c(paste("comp",c(1:nbcomp1), "_", lab[1], sep = ""),  paste("comp",c(1:nbcomp2), "_", lab[2], sep = "")))
            dimnames(r) <- list(c(paste("comp",c(1:nbcomp1), "_", lab[1], sep = ""),  paste("comp",c(1:nbcomp2), "_", lab[2], sep = "")),c(paste("comp",c(1:nbcomp1), "_", lab[1], sep = ""),  paste("comp",c(1:nbcomp2), "_", lab[2], sep = "")))

        }
      return(list(cor=r,pval=rpval,inter = gg, nbComp = c(nbcomp1,nbcomp2), labAn = lab))
}





      
correl2Comp <- function (
### This function computes the correlation between two components.
                         comp1
                         ### The first component, a vector of projections or contributions indexed by labels
                         ,comp2
                         ### The second component, a vector of projections or contributions indexed by labels
                         ,type.corr = "pearson"
                         ###  name of the type of correlation to compute, either 'pearson' or 'spearman'
                         ,plot = FALSE
                         ###  if TRUE the plot of comp1 vs comp2 is drawn
                         ,cutoff_zval = 0
                         ###  either NULL or 0 (default) if all genes are used to compute the correlation between the components, or a threshold to compute the correlation on the genes that have at least a scaled projection higher than cutoff_zval. 
                         ,test = FALSE
                         ### if TRUE the p-value of correlation is returned instead of the correlation value
                         ,alreadyTreat = FALSE
                         ### if TRUE comp1 and comp2 are considered as being already treated (i.e scaled and restricted to common elements) 
                         ) {
##details<< Before computing the correlation, the components are scaled and restricted to common labels.  
##  When \code{cutoff_zval} is different from 0, the elements that are included in the circle of center 0 and radius \code{cutoff_zval} are not taken into account during the computation of the correlation.

    if(!alreadyTreat) {
        comp1 <- (comp1-mean(comp1))/sd(comp1)
        comp2 <- (comp2-mean(comp2))/sd(comp2)
        gg <- intersect(names(comp1),names(comp2))
        comp1 <- comp1[gg]
        comp2 <- comp2[gg]
    }

            dc <- data.frame(comp1=comp1,comp2=comp2)
            rownames(dc) <- gg
            g2del <- rownames(subset(dc, comp1^2 + comp2^2<= cutoff_zval   ))
            genes.intersect <- gg[!(gg %in% g2del)]
        
        if (length(genes.intersect) < 4)
            if (test) return(1) else return(0)
	comp1.ord = comp1[genes.intersect]
	comp2.ord = comp2[genes.intersect]
        
	
	if (test) {
                
                r = cor.test(comp1.ord,comp2.ord,method = type.corr, alternative = "two.sided")$p.value
	}
	else {
                r = cor(comp1.ord,comp2.ord,method = type.corr)
                r = signif(r, digits = 3)
		
	}
	if (plot) plot(comp1.ord,comp2.ord,xlab = paste("comp1"),ylab=paste("comp2"),pch = 16,sub = paste("corr_",type.corr," = ",r,sep=""))
	return(r)
### This function returns either the correlation value or the p-value of the correlation test.
}


##' This function builds a data.frame describing for each node of the graph
##' its ID and which analysis/data it comes from.
##'
##' The created file is used in Cytoscape.
##' @title Generate node attributes
##' @param nbAn Number of analyses being considered, i.e number of IcaSet
##' objects
##' @param nbComp Number of components by analysis, if of length 1 then
##' it is assumed that each analysis has the same number of components.
##' @param labAn Labels of the analysis, if missing it will be generated
##' as an1, an2, ...
##' @param labComp List containing the component labels indexed by analysis, if missing
##' will be generated as comp1, comp2, ...
##' @param file File where the description of the node attributes will be
##' written
##' @return A data.frame describing each node/component
##' @export
##' @examples
##' ## 4 datasets, 20 components calculated in each dataset, labAn    
##' nodeAttrs(nbAn=4, nbComp=20, labAn=c("tutu","titi","toto","tata"))
##' 
##' @author Anne Biton

nodeAttrs <- function(nbAn,
                      nbComp,
                      labAn,
                      labComp,
                      file
                      ) {

    nb <- an <- lab <- comp <- NULL
    
    if (missing(labAn))
        labAn <- paste(rep("an",nbAn),1:nbAn,sep = "")
    else if (length(labAn) != nbAn)
        stop("The length of 'labAn' must equal 'nbAn'.")

    
    if (missing(labComp)) {
        if (length(nbComp) == 1) {
            labid <- 
                paste(
                      rep(paste(rep("comp",nbComp),c(1:nbComp),sep = ""),nbAn),
                      paste(rep("an",nbComp*nbAn),sapply(c(1:nbAn),rep,nbComp),sep=""),
                      sep = "_"
                      )
            labComp <-rep(paste(rep("comp",nbComp),c(1:nbComp),sep = ""),nbAn)
                      
            nbComp <- rep(nbComp,nbAn)
        }
        else {
            if (length(nbComp) != nbAn)
                stop("The provided number of components must equal the number of analyses")
           
            labid <- foreach(nb = nbComp, an = labAn) %do% {paste(paste(rep("comp",nb),c(1:nb),sep = ""),rep(an,nb),sep="_")}
            labid <- unlist(labid)
            labComp <- unlist(foreach(nb = nbComp, an = 1:nbAn) %do% {paste(rep("comp",sum(nb)),c(1:nb),sep = "")})
        }
    }
    else {
        if (!is.list(labComp))
            stop("The component labels must be provided using a list.")
        if (length(labComp) != nbAn)
            stop("The length of the list containing the component labels 'labComp' must have the same length than the list containing the analysis labels 'labAn'.")
        names(labComp) <- names(labAn)

        foreach(lab = labComp, nb = nbComp) %do% {
            if (length(lab) != nb)
                stop("The length of the list 'labComp' containing component labels does not map to 'nbComp'.")
        }
        labid <- unlist(foreach(nb = nbComp, an = labAn, comp = labComp) %do% {paste(comp,rep(an,nb),sep="_")})
        labComp <- unlist(labComp)

        
    }

   compnb <- unlist(lapply(nbComp,function(x) c(1:x)))
   names(nbComp) <- labAn
   an <- unlist(lapply(names(nbComp), function(x,nb) rep(x,nb[x]), nb = nbComp))

   dataAttr <- data.frame(id = labid, labComp = labComp, indComp = compnb, labAn = an, stringsAsFactors = FALSE, check.names = FALSE)

    if (!missing(file))
        if (!is.null(file))
            write.table(dataAttr,file=file,sep = " ", row.names = FALSE, quote = FALSE)

    return(dataAttr)
}

##' compareAn2graphfile
##' 
##' This function builds a correlation graph from the outputs of function \code{\link{compareAn}}.
##'
##' When correlations are considered (\code{useVal}="cor"), absolute values
##' are used since the components have no direction.
##' 
##' If \code{useMax} is \code{TRUE} each component is linked to the most correlated component of
##' each different \code{IcaSet}.
##' 
##' If \code{cutoff} is specified, only
##' correlations exceeding this value are taken into account during the graph construction.
##' For example, if \code{cutoff} is 1, only relationships between components
##' that correspond to a correlation value larger than 1 will be included.
##'
##' When \code{useVal="pval"} and \code{useMax=TRUE}, the minimum value is taken
##' instead of the maximum. 
##' 
##' @param listPairCor The output of the function \code{\link{compareAn}}, containing the
##' correlation between several pairs of objects of class \code{\link{IcaSet}}.
##' @param useMax If TRUE, the graph is restricted to edges that correspond to
##' maximum score, see details
##' @param cutoff Cutoff used to select pairs that will be included in the
##' graph.
##' @param useVal The value on which is based the graph, either \code{"cor"} for
##' correlation or \code{"pval"} for p-values of correlation tests.
##' @param file File name.
##' @return A data.frame with the graph description, has
##' two columns \code{n1} and \code{n2} filled with node IDs, each row denotes that there is an edge from \code{n1} to \code{n2}. Additional columns quantify the strength of association: correlation (\code{cor}), p-value (\code{pval}), (\code{1-abs(cor)}) (\code{distcor}), log10-pvalue (\code{logpval}).
##' @author Anne Biton
##' @seealso \code{\link{compareAn}}, \code{\link{cor2An}}
##' @export
##' @examples
##'
##' dat1 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat1) <- paste("g", 1:1000, sep="")
##' colnames(dat1) <- paste("s", 1:10, sep="")
##' dat2 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat2) <- paste("g", 1:1000, sep="")
##' colnames(dat2) <- paste("s", 1:10, sep="")
##' 
##' ## run ICA
##' resJade1 <- runICA(X=dat1, nbComp=3, method = "JADE")
##' resJade2 <- runICA(X=dat2, nbComp=3, method = "JADE")
##' 
##' ## build params
##' params <- buildMineICAParams(resPath="toy/")
##'
##' ## build IcaSet object
##' icaSettoy1 <- buildIcaSet(params=params, A=data.frame(resJade1$A), S=data.frame(resJade1$S),
##'                           dat=dat1, alreadyAnnot=TRUE)$icaSet
##' icaSettoy2 <- buildIcaSet(params=params, A=data.frame(resJade2$A), S=data.frame(resJade2$S),
##'                           dat=dat2, alreadyAnnot=TRUE)$icaSet
##' 
##' resCompareAn <- compareAn(icaSets=list(icaSettoy1,icaSettoy2), labAn=c("toy1","toy2"), 
##'                          type.corr="pearson", level="genes", cutoff_zval=0)
##'
##' ## Build a graph where edges correspond to maximal correlation value (useVal="cor"),
##' compareAn2graphfile(listPairCor=resCompareAn, useMax=TRUE, useVal="cor", file="myGraph.txt")
##' 
##' 
##' \dontrun{
##' #### Comparison of 2 ICA decompositions obtained on 2 different gene expression datasets.
##' ## load the two datasets
##' library(breastCancerMAINZ)
##' library(breastCancerVDX)
##' data(mainz)
##' data(vdx)
##' 
##' ## Define a function used to build two examples of IcaSet objects
##' treat <- function(es, annot="hgu133a.db") {
##'    es <- selectFeatures_IQR(es,10000)
##'    exprs(es) <- t(apply(exprs(es),1,scale,scale=FALSE))
##'    colnames(exprs(es)) <- sampleNames(es)
##'    resJade <- runICA(X=exprs(es), nbComp=10, method = "JADE", maxit=10000)
##'    resBuild <- buildIcaSet(params=buildMineICAParams(), A=data.frame(resJade$A), S=data.frame(resJade$S),
##'                         dat=exprs(es), pData=pData(es), refSamples=character(0),
##'                         annotation=annot, typeID= typeIDmainz,
##'                         chipManu = "affymetrix", mart=mart)
##'    icaSet <- resBuild$icaSet
##' }
##' ## Build the two IcaSet objects
##' icaSetMainz <- treat(mainz)
##' icaSetVdx <- treat(vdx)
##' 
##' ## Compute correlation between every pair of IcaSet objects.
##' resCompareAn <- compareAn(icaSets=list(icaSetMainz,icaSetVdx),
##' labAn=c("Mainz","Vdx"), type.corr="pearson", level="genes", cutoff_zval=0)
##' 
##' ## Same thing but adding a selection of genes on which the correlation between two components is computed:
##' # when considering pairs of components, only projections whose scaled values are not located within
##' # the circle of radius 1 are used to compute the correlation (cutoff_zval=1).
##' resCompareAn <-  compareAn(icaSets=list(icaSetMainz,icaSetVdx),
##' labAn=c("Mainz","Vdx"), type.corr="pearson", cutoff_zval=1, level="genes")
##'
##' ## Build a graph where edges correspond to maximal correlation value (useVal="cor"),
##' ## i.e, component A of analysis i is linked to component B of analysis j,
##' ## only if component B is the most correlated component to A amongst all component of analysis j.
##' compareAn2graphfile(listPairCor=resCompareAn, useMax=TRUE, useVal="cor", file="myGraph.txt")
##' 
##' ## Restrict the graph to correlation values exceeding 0.4
##' compareAn2graphfile(listPairCor=resCompareAn, useMax=FALSE, cutoff=0.4,  useVal="cor", file="myGraph.txt")
##'
##' }
compareAn2graphfile <- function (listPairCor,
                                 useMax = TRUE,
                                 cutoff = NULL,
                                 useVal = c("cor","pval"),
                                 file = NULL
                                 ) {
    useVal <- match.arg(useVal)
    dataGraphs <- 
    lapply(listPairCor,
           function(res, useMax, cutoff, useVal) {
               dd <- abs(res[[useVal]])
               
               nbcomp1 <- res$nbComp[1]
               nbcomp2 <- res$nbComp[2]
               attrGraph <- c("cor",if ("pval" %in% names(res)) "pval")
               if (useMax) {
                   if (useVal == "pval") 
                       op <- "which.min"
                   else 
                       op <- "which.max"
                   
                  
                   ## find correlation max both direction an1 <-> an2 since the max is not necessarily reciprocal
                   subdd <- dd[(nbcomp1+1):(nbcomp1+nbcomp2), 1:nbcomp1]
                   ind <- apply(subdd,2,op)
                   coord <- as.list(data.frame(t(cbind(ind,1:nbcomp1))))
                   
                   dataGraph <- data.frame(n1=rownames(dd)[1:nbcomp1],
                                           n2=rownames(dd)[(nbcomp1+1):(nbcomp1+nbcomp2)][ind])
                   dataGraph[[useVal]] <- unlist(lapply(coord, function(x,m) m[x[1],x[2]], m = subdd))
                   addval <- setdiff(attrGraph,useVal)

                   if (length(addval)==1)
                       dataGraph[[addval]] <- unlist(lapply(coord, function(x,m) m[x[1],x[2]], m = res[[addval]][(nbcomp1+1):(nbcomp1+nbcomp2), 1:nbcomp1]))

                   ## an2 -> an1
                   subdd <- dd[1:nbcomp1,(nbcomp1+1):(nbcomp1+nbcomp2)]
                   ind <- apply(subdd,2,op)
                   coord <- as.list(data.frame(t(cbind(ind,1:nbcomp2))))
                   
                   dataGraph2 <- data.frame(n1=rownames(dd)[(nbcomp1+1):(nbcomp1+nbcomp2)],
                                           n2=rownames(dd)[1:nbcomp1][ind])

                   dataGraph2[[useVal]] <- unlist(lapply(coord, function(x,m) m[x[1],x[2]], m = subdd))

                   addval <- setdiff(attrGraph,useVal)

                   if (length(addval)==1)
                       dataGraph2[[addval]] <- unlist(lapply(coord, function(x,m) m[x[1],x[2]], m = res[[addval]][1:nbcomp1,(nbcomp1+1):(nbcomp1+nbcomp2)]))

                   dataGraph <- rbind(dataGraph,dataGraph2)

                   pval <- NULL
                   if (!is.null(cutoff)) {
                       if (useVal == "pval")
                           dataGraph <- subset(dataGraph, pval <= cutoff)
                       else
                           dataGraph <- subset(dataGraph, cor >= cutoff)
                           
                   }
               }
               else if (!is.null(cutoff)) {
                   subdd <- dd[1:nbcomp1, (nbcomp1+1):(nbcomp1+nbcomp2)]
                   dataGraph <- 
                       lapply(as.list(colnames(subdd)), 
                         function(nn, cutoff, subdd,  useVal, attrGraph, res,op, subInd)   {
                             x <- subdd[,nn]

                             ind <- which(eval(parse(text=paste("x",op,"cutoff")))) 
                             if (length(ind) > 0) {
                                 dataGraph <- data.frame(n1 = rep(nn,length(ind)), n2 = rownames(subdd)[ind])
                                 dataGraph[[useVal]] <- subdd[ind,nn]
                                 addval <- setdiff(attrGraph,useVal)
                                 if (length(addval)==1) {
                                    coord <- as.list(data.frame(t(cbind(dataGraph$n1,dataGraph$n2))))
                                    dataGraph[[addval]] <-  unlist(lapply(coord, function(x,m) m[x[1],x[2]], m = res[[addval]][subInd[[1]],subInd[[2]]]))
                                 }
                                 
                                 return(dataGraph)
                                 
                             }
                          }
                          , cutoff = cutoff
                          , subdd = subdd
                          , useVal = useVal
                          , attrGraph = attrGraph
                          , res = res
                          , op = if (useVal == "pval") "<=" else ">="
                          , subInd = list(1:nbcomp1, (nbcomp1+1):(nbcomp1+nbcomp2))
                          )
                   dataGraph <- do.call(rbind,dataGraph)
                   

               }
               return(dataGraph)                        

           }
           , useMax = useMax
           , cutoff = if (missing(cutoff)) NULL else cutoff 
           , useVal = useVal 
           )
    dataGraphs <- do.call(rbind,dataGraphs)

    
   dataGraphs$invcor = 1/(10*as.numeric(dataGraphs$cor))
   dataGraphs$distcor = 1-abs(as.numeric(dataGraphs$cor))
   dataGraphs[,c("cor","pval","invcor","distcor")] <- apply(dataGraphs[,c("cor","pval","invcor","distcor")],2,as.numeric)
   dataGraphs <- annotReciprocal(dataGraph = dataGraphs, keepOnlyReciprocal = FALSE)# file  
   dataGraphs[,c("cor","pval","invcor","distcor")] <- apply(dataGraphs[,c("cor","pval","invcor","distcor")],2,as.numeric)

    if (!missing(file))
        if (!is.null(file))
            write.table(dataGraphs, file = file, row.names = FALSE, quote = FALSE, sep = "\t")
    
   return(dataGraphs)

}

##' annotReciprocal
##' 
##' This function notes edges of a graph as reciprocal or not.
##' 
##' 
##' @param dataGraph data.frame which contains the graph description, must have
##' two columns \code{n1} and \code{n2} filled with node IDs, each row denoting there is an edge from \code{n1} to \code{n2}.
##' @param file file where the graph description is written
##' @param keepOnlyReciprocal if TRUE \code{dataGraph} is restricted to
##' reciprocal edges, else all edges are kept (default).
##' @return This function returns the argument \code{dataGraph} with an
##' additional column named 'reciprocal' which contains TRUE if the edge
##' described by the row is reciprocal, and FALSE if it is not reciprocal.
##' @examples 
##' dg <- data.frame(n1=c("A","B","B","C","C","D","E","F"),n2=c("B","G","A","B","D","C","F","E"))
##' annotReciprocal(dataGraph=dg)
##' 
##' @export
##' @author Anne Biton
annotReciprocal <- function (dataGraph,
                             file,
                             keepOnlyReciprocal = FALSE
                             ) {

	names(dataGraph)[1:2] <- c("n1","n2") 
        dataGraph$keep <- rep(NA, length=nrow(dataGraph))

	dataGraph_keep <- apply( dataGraph, MARGIN = 1 , 
					   function (r, dataGraph) {
						dataN2 = subset(dataGraph, dataGraph$n1 == r["n2"])
						if (r["n1"] %in% dataN2$n2) r["keep"] = TRUE
						else r["keep"] = FALSE
						return(r)		
					    }
					    , dataGraph = dataGraph
	)

        dataGraph_keep <- as.data.frame(t(dataGraph_keep), check.names = FALSE, stringsAsFactors = FALSE)

        names(dataGraph_keep)[ncol(dataGraph_keep)] <- "reciprocal"
        
        if (keepOnlyReciprocal) {
          dataGraph_keep<- subset(dataGraph_keep, dataGraph_keep$reciprocal == TRUE )
          dataGraph_keep <- dataGraph_keep[,-ncol(dataGraph_keep)]
        }
        
        if (!missing(file))
            if(!is.null(file)) 
                write.table(dataGraph_keep, file  = file, sep = "\t", quote = FALSE, row.names = FALSE)
        
	return(dataGraph_keep)
}
 

##' This function plots the
##' correlation graph in an interactive device using function \code{tkplot}.
##'
##' You have to slighly move the nodes to see cliques because strongly related nodes are often superimposed. 
##' The \code{edgeWeight} column is used to weight the edges within the
##' fruchterman.reingold layout available in the package \code{igraph}.
##'
##' The argument
##' \code{nodeCol} typically denotes the column containing the names of the datasets.
##' Colors are automatically
##' attributed to the nodes using palette Set3 of package \code{RColorBrewer}. The
##' corresponding colors can be directly specified in the 'col' argument. In
##' that case, 'col' must be a vector of colors indexed by the unique elements
##' contained in \code{nodeCol} column (e.g dataset ids).
##'
##' As for colors, one can define
##' the column of \code{nodeAttrs} that is used to define the node shapes.  The
##' corresponding shapes can be directly specified in the \code{shape} argument. In
##' that case, \code{shape} must be one of \code{c("circle","square", " vcsquare", "rectangle", "crectangle", "vrectangle")} and must be 
##' indexed by the unique elements of \code{nodeShape} column.
##'
##' Unfortunately, shapes
##' can't be taken into account when tkplot is TRUE (interactive plot).
##'
##' If \code{reciproCol} is not missing, it is used to color the edges, either in grey if the
##' edge is not reciprocal or  in black if the edge is reciprocal.
##'
##' 
##' @title Plots graph using 
##' @param dataGraph A data.frame containing the graph description. It must
##' have two columns \code{n1} and \code{n2}, each row denoting that there is an edge from n1
##' to n2.  Node labels in columns \code{n1} and \code{n2} of \code{dataGraph} must correspond to
##' node IDs in column \code{id} of \code{nodeAttrs}.
##' @param edgeWeight The column of dataGraph used to weight edges.
##' @param nodeAttrs A data.frame with node description, see function \code{nodeAttrs}.
##' @param nodeShape Denotes the column of \code{nodeAttrs} used to attribute the node shapes. 
##' @param nodeCol Denotes the column of \code{nodeAttrs} used to
##' color the nodes in the graph.
##' @param nodeName Denotes the column of \code{nodeAttrs} used
##' as labels for the nodes in the graph.
##' @param col A vector of colors, for the nodes, indexed by the unique elements of \code{nodeCol}
##' column from \code{nodeAttrs}. If missing, colors will be automatically attributed.
##' @param shape A vector of shapes indexed by the unique elements of
##' column \code{nodeShape} from \code{nodeAttrs}. If missing, shapes will be automatically
##' attributed.
##' @param title Title for the plot
##' @param reciproCol Denotes the column of \code{dataGraph} containing \code{TRUE} if the
##' row defines a reciprocal node, else \code{FALSE}. See \code{\link{annotReciprocal}}.
##' @param tkplot If TRUE, performs interactive plot with function \code{tkplot}, else uses \code{plot.igraph}.
##' @param \dots Additional parameters as required by \code{tkplot}.
##' @return A list consisting of \describe{
##' \item{dataGraph}{a data.frame defining the correlation
##' graph}
##' \item{nodeAttrs}{a data.frame describing the node
##' of the graph}
##' \item{graph}{the graph as an object of class \code{igraph}}
##' \item{graphid}{the id of the graph plotted using \code{tkplot}} }
##' @export
##' @author Anne Biton
##' @seealso \code{\link{compareAn}}, \code{\link{nodeAttrs}}, \code{\link{compareAn2graphfile}}, \code{\link{runCompareIcaSets}}
##' @examples
##'
##' dat1 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat1) <- paste("g", 1:1000, sep="")
##' colnames(dat1) <- paste("s", 1:10, sep="")
##' dat2 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat2) <- paste("g", 1:1000, sep="")
##' colnames(dat2) <- paste("s", 1:10, sep="")
##' 
##' ## run ICA
##' resJade1 <- runICA(X=dat1, nbComp=3, method = "JADE")
##' resJade2 <- runICA(X=dat2, nbComp=3, method = "JADE")
##' 
##' ## build params
##' params <- buildMineICAParams(resPath="toy/")
##'
##' ## build IcaSet object
##' icaSettoy1 <- buildIcaSet(params=params, A=data.frame(resJade1$A), S=data.frame(resJade1$S),
##'                           dat=dat1, alreadyAnnot=TRUE)$icaSet
##' icaSettoy2 <- buildIcaSet(params=params, A=data.frame(resJade2$A), S=data.frame(resJade2$S),
##'                           dat=dat2, alreadyAnnot=TRUE)$icaSet
##' icaSets <- list(icaSettoy1, icaSettoy2)
##' 
##' resCompareAn <- compareAn(icaSets=list(icaSettoy1,icaSettoy2), labAn=c("toy1","toy2"), 
##'                          type.corr="pearson", level="genes", cutoff_zval=0)
##'
##' ## Build a graph where edges correspond to maximal correlation value (useVal="cor"),
##' dataGraph <- compareAn2graphfile(listPairCor=resCompareAn, useMax=TRUE, useVal="cor", file="myGraph.txt")
##' 
##' ## construction of the data.frame with the node description
##' nbComp <- rep(3,2) #each IcaSet contains 3 components
##' nbAn <- 2 # we are comparing 2 IcaSets
##' # labels of components created as comp*i* 
##' labComp <- foreach(icaSet=icaSets, nb=nbComp, an=1:nbAn) %do% {
##'                   paste(rep("comp",sum(nb)),1:nbComp(icaSet),sep = "")}
##' 
##' # creation of the data.frame with the node description
##' nodeDescr <- nodeAttrs(nbAn = nbAn, nbComp = nbComp, labComp = labComp,
##'                        labAn = c("toy1","toy2"), file = "nodeInfo.txt")
##'
##' ## Plot correlation graph, slightly move the attached nodes to make the cliques visible
##' ## use tkplot=TRUE to have an interactive graph
##' res <- plotCorGraph(title = "Compare toy 1 and 2", dataGraph = dataGraph, nodeName = "indComp", tkplot = FALSE, 
##'                  nodeAttrs = nodeDescr, edgeWeight = "cor", nodeShape = "labAn", reciproCol = "reciprocal")
##'
##' 
##' \dontrun{
##' ## load two microarray datasets
##' library(breastCancerMAINZ)
##' library(breastCancerVDX)
##' data(mainz)
##' data(vdx)
##' 
##' ## Define a function used to build two examples of IcaSet objects
##' treat <- function(es, annot="hgu133a.db") {
##'    es <- selectFeatures_IQR(es,10000)
##'    exprs(es) <- t(apply(exprs(es),1,scale,scale=FALSE))
##'    colnames(exprs(es)) <- sampleNames(es)
##'    resJade <- runICA(X=exprs(es), nbComp=10, method = "JADE", maxit=10000)
##'    resBuild <- buildIcaSet(params=buildMineICAParams(), A=data.frame(resJade$A), S=data.frame(resJade$S),
##'                         dat=exprs(es), pData=pData(es), refSamples=character(0),
##'                         annotation=annot, typeID= typeIDmainz,
##'                         chipManu = "affymetrix", mart=mart)
##'    icaSet <- resBuild$icaSet
##' }
##' ## Build the two IcaSet objects
##' icaSetMainz <- treat(mainz)
##' icaSetVdx <- treat(vdx)
##' 
##' icaSets <- list(icaSetMainz, icaSetVdx)
##' labAn <- c("Mainz", "Vdx")
##'
##' ## correlations between gene projections of each pair of IcaSet
##' resCompareAn <- compareAn(icaSets = icaSets, level = "genes", type.corr= "pearson",
##'                           labAn = labAn, cutoff_zval=0)
##'
##' ## construction of the correlation graph using previous output
##' dataGraph <- compareAn2graphfile(listPairCor=resCompareAn, useMax=TRUE, file="corGraph.txt")
##'
##' ## construction of the data.frame with the node description
##' nbComp <- rep(10,2) #each IcaSet contains 10 components
##' nbAn <- 2 # we are comparing 2 IcaSets
##' # labels of components created as comp*i* 
##' labComp <- foreach(icaSet=icaSets, nb=nbComp, an=1:nbAn) %do% {
##'                   paste(rep("comp",sum(nb)),1:nbComp(icaSet),sep = "")}
##' 
##' # creation of the data.frame with the node description
##' nodeDescr <- nodeAttrs(nbAn = nbAn, nbComp = nbComp, labComp = labComp,
##'     labAn = labAn, file = "nodeInfo.txt")
##'
##' ## Plot correlation graph, slightly move the attached nodes to make the cliques visible
##' res <- plotCorGraph(title = "Compare two ICA decomsitions obtained on \n two
##'                  microarray-based data of breast tumors", dataGraph = dataGraph, nodeName = "indComp", 
##'                  nodeAttrs = nodeDescr, edgeWeight = "cor", nodeShape = "labAn", reciproCol = "reciprocal")
##'
##' }
plotCorGraph <- function(
                      dataGraph,
                      edgeWeight = "cor",
                      nodeAttrs,
                      nodeShape,
                      nodeCol = "labAn",
                      nodeName = "indComp",
                      col,
                      shape,
                      title = "",
                      reciproCol = "reciprocal",
                      tkplot = FALSE,
                      ...
                      ) {
    if (!(edgeWeight %in% colnames(dataGraph))) {
        stop (paste(edgeWeight,"is not available within the columns of dataGraph."))
    }
    nodeName <- match.arg(nodeName, choices = colnames(nodeAttrs))
    

    if (missing(nodeShape)) 
        nodeAttrs$shape <- rep("circle",nrow(nodeAttrs))
    else if (!(nodeShape %in% colnames(nodeAttrs))) {
        nodeAttrs$shape <-  rep("circle",nrow(nodeAttrs))
        warning(paste("The column ", nodeShape, " is not available in 'dataGraph'.",sep = ""))
    }
    else {
        potShapes <- c("circle","square", "csquare", "rectangle", "crectangle", "vrectangle")

        autoShape <- TRUE
        if (!missing(shape)) {
            autoShape <- FALSE
            if (length(intersect(names(shape),unique(nodeAttrs[[nodeShape]]))) != length(unique(nodeAttrs[[nodeShape]]))) { 
                warning("'shape' argument is incorrect, some elements of nodeAttrs$'nodeShape' are not indexed. Shapes will be attributed automatically.")
                autoShape <- TRUE
            }
            
            if (length(setdiff(shape,potShapes))>0) {
                warning("'shape' argument is incorrect, the available shapes are: 'circle','square', 'csquare', 'rectangle', 'crectangle', and 'vrectangle'. Shapes will be attributed automatically.")
                autoShape <- TRUE
                
            }
        }
        if (autoShape) {
            levs <- unique(nodeAttrs[[nodeShape]])
            if (length(levs) > length(potShapes)) {
                warning("R graphs can only handle 6 node shapes at the most.")
                shape <- rep(potShapes,6)[1:length(levs)]
            }
            else
                shape <- potShapes[1:length(levs)]
            names(shape) <- levs
            nodeAttrs$shape <- shape[nodeAttrs[[nodeShape]]]
        }
        else 
            nodeAttrs$shape <- shape[nodeAttrs[[nodeShape]]]
        
    }
        
    
    if (!(nodeCol %in% colnames(nodeAttrs))) { 
        stop(paste("The column ",nodeCol, " is not available in nodeAttrs.",sep=""))
    }
    else {
        nbAn <- length(unique(nodeAttrs[[nodeCol]]))
        
        if (!missing(col) && !is.null(col)) {
            if (length(intersect(names(col),unique(nodeAttrs[[nodeCol]]))) != length(unique(nodeAttrs[[nodeCol]]))) {
                warning("'col' argument is incorrect, some elements of nodeAttrs$'nodeCol' are not indexed. Colors will be attributed automatically.")
                colAn <- brewer.pal(nbAn,"Set3")
                names(colAn) <- unique(nodeAttrs[[nodeCol]])
                nodeAttrs$col = colAn[nodeAttrs[[nodeCol]]]
            }
            else
                nodeAttrs$col <- col[nodeAttrs[[nodeCol]]]
        }
        else {
            colAn <- brewer.pal(nbAn,"Set3")
            names(colAn) <- unique(nodeAttrs[[nodeCol]])
            nodeAttrs$col = colAn[nodeAttrs[[nodeCol]]]
        }
    }
    
    
    n1 <- NULL
    nodes <- as.character(unique(as.character(dataGraph$n1)))
    edges <- llply(nodes,
                   function(n,data,edgeWeight) {
                       ll <- list(edges = as.character(subset(data, n1 == n)$n2), weights = subset(data, n1 == n)[[edgeWeight]])
                   },
                   data = dataGraph, edgeWeight = edgeWeight)
    names(edges) <- nodes
    
    
    
    
    g <- new("graphNEL", nodes = nodes, edgeL = edges, edgemode = "directed")
                                        # build igraph object from graphNEL object
    ig <- igraph.from.graphNEL(g, name = TRUE, weight = TRUE)
    
    V(ig)$color <- nodeAttrs$col[match(V(ig)$name,nodeAttrs$id)]
    V(ig)$shape <- nodeAttrs$shape[match(V(ig)$name,nodeAttrs$id)]
    E(ig)$width <- 10^E(ig)$weight
    indHighCor <- which(abs(E(ig)$weight)>0.4)
    indLowCor <- which(abs(E(ig)$weight)<=0.4)
    E(ig)$weight <- 1/abs(log2(abs(E(ig)$weight)))
    E(ig)$weight[indHighCor] <- E(ig)$weight[indHighCor]+6
    #E(ig)$weight[indLowCor] <- E(ig)$weight[indLowCor]-0.25

   if (!missing(reciproCol)) {
       if (reciproCol %in% colnames(dataGraph)) {

           nodes <- as.character(unique(as.character(dataGraph$n1)))
           colEdges <- llply(nodes,
                             function(n,data,reciproCol) {
                              c('TRUE' = "black", 'FALSE' = "gray50")[as.character(subset(data, n1 == n)[[reciproCol]])]
                             },
                             data = dataGraph, reciproCol = reciproCol)

           colEdges <- unlist(colEdges)
           E(ig)$color <- colEdges
       }
       else {
            warning(paste("No column of 'dataGraph' has the name",reciproCol))
            E(ig)$color <- "black"
        }

    }
    else
        E(ig)$color <- "black"
    
    lay <- layout.fruchterman.reingold(ig,niter=500,area=vcount(ig)^2.3,repulserad=vcount(ig)^100, weights = E(ig)$weight)

    graph <- ig
            
    if (!capabilities()[["X11"]])
        tkplot <- FALSE
    
    if (tkplot) {
        
        if(capabilities()[["X11"]] & tolower(Sys.info()["sysname"])!="windows") {
            dimScreen <- system("xdpyinfo | grep dimensions", intern = TRUE)
            dimScreen <- as.numeric(strsplit(gsub(gsub(gsub(dimScreen, pattern = " ", replacement = ""), pattern = "[A-Za-z0-9]*:", replacement = ""), pattern = "pixels\\([A-Za-z0-9]*\\)", replacement = ""), split = "x")[[1]])
        }
        else
            dimScreen <- c(450,450)
        
        graphid <- tkplot(ig,
                          layout = lay,
                          vertex.label =  nodeAttrs[[nodeName]][match(V(ig)$name,nodeAttrs$id)],
                          vertex.shape = nodeAttrs$shape[match(V(ig)$name,nodeAttrs$id)],
                          canvas.width=dimScreen[1], canvas.height=dimScreen[2],
                          ...)
        
        
        tkplot.fit.to.screen(graphid)
        
    }
    else {
        
        graph <- plot(ig,
                      layout = lay,
                      vertex.label = nodeAttrs[[nodeName]][match(V(ig)$name,nodeAttrs$id)],
                      vertex.shape = nodeAttrs$shape[match(V(ig)$name,nodeAttrs$id)],
                      main = title,
                      ...
                      )
        graph <- ig
        graphid = ""

        
    }
    return(list(dataGraph=dataGraph, nodeAttrs = nodeAttrs, graph = graph, graphid=graphid))

}



##' runCompareIcaSets
##' 
##' This function encompasses the comparison of several IcaSet objects using correlations
##' and the plot of the corresponding correlation graph.
##' The IcaSet objects are compared by calculating the correlation between either
##' projection values of common features or genes, or contributions of common
##' samples.
##' 
##' This function calls four functions: \code{\link{compareAn}} which computes the
##' correlations, \code{\link{compareAn2graphfile}} which builds the graph,
##' \code{\link{nodeAttrs}} which builds the node description data, and \code{\link{plotCorGraph}}
##' which uses tkplot to plot the graph in an interactive device.
##' 
##' If the user wants to see the correlation graph in Cytoscape, he must fill
##' the arguments \code{fileDataGraph} and \code{fileNodeDescr}, in order to
##' import the graph and its node descriptions as a .txt file in Cytoscape.
##' 
##' When \code{labAn} is missing, each element i of \code{icaSets} is labeled as
##' 'Ani'.
##' 
##' The user must carefully choose the data level used in the comparison:
##' If \code{level='samples'}, the correlations are based on the mixing
##' matrices of the ICA decompositions (of dimension samples x components). \code{'A'}
##' will be typically chosen when the ICA decompositions were computed on the
##' same dataset, or on datasets that include the same samples.
##' If \code{level='features'} is
##' chosen, the correlation is calculated between the source matrices (of dimension
##' features x components) of the ICA decompositions. \code{'S'} will be typically
##' used when the ICA decompositions share common features (e.g same microarrays).
##' If \code{level='genes'}, the correlations are calculated
##' on the attributes \code{'SByGene'} which store the 
##' projections of the annotated features. \code{'SByGene'} will be typically chosen
##' when ICA were computed on datasets from different technologies, for which
##' comparison is possible only after annotation into a common ID, like genes.
##'
##' \code{cutoff_zval}
##' is only used when \code{level} is one of \code{c('features','genes')}, in
##' order to restrict the correlation to the contributing features or genes.
##' 
##' When \code{cutoff_zval} is specified, for each pair of components, genes or features that are
##' included in the circle of center 0 and radius \code{cutoff_zval} are excluded from
##' the computation of the correlation.
##'
##' It must be taken into account by the user that if cutoff_zval
##' is different from NULL or zero, the computation will be much slowler since
##' each pair of component is treated individually.
##' 
##' Edges of the graph are built based on the correlation values between the
##' components. Absolute values of correlations are used since
##' components have no direction.
##' 
##' If \code{useMax} is \code{TRUE} each component will be linked to only one component of
##' each other IcaSet that corresponds to the most correlated component among
##' all components of the same IcaSet. If \code{cutoff_graph} is specified, only
##' correlations exceeding this value are taken into account to build the graph.
##' For example, if \code{cutoff} is 1, only relationships between components
##' that correspond to a correlation value higher than 1 will be included.
##' Absolute correlation values are used since the components have no direction.
##' 
##' The contents of the returned list are \describe{
##' \item{dataGraph:}{\code{dataGraph} data.frame that describes the correlation
##' graph,} \item{nodeAttrs:}{\code{nodeAttrs} data.frame that describes the node
##' of the graph} \item{graph}{\code{graph} the graph as an igraph-object,}
##' \item{graphid:}{\code{graphid} the id of the graph plotted using tkplot.} }
##' 
##' @param icaSets List of \code{\link{IcaSet}} objects, e.g results of ICA decompositions
##' obtained on several datasets.
##' @param labAn Vector of names for each icaSet, e.g the the names of the
##' datasets on which were calculated the decompositions.
##' @param type.corr Type of correlation to compute, either
##' \code{'pearson'} or \code{'spearman'}.
##' @param cutoff_zval Either NULL or 0 (default) if all genes are used to
##' compute the correlation between the components, or a threshold to compute
##' the correlation using the genes that have at least a scaled projection higher
##' than cutoff_zval. Will be used only when \code{level} is one of \code{c("features","genes")}.
##' @param level Data level of the \code{IcaSet} objects on which is applied the correlation.
##' It must correspond to a data level shared by the IcaSet objects:
##' \code{'samples'} if they were applied to common samples (correlations are computed between matrix
##' \code{A}), \code{'features'} if they were applied to common features (correlations are computed between matrix \code{S}), \code{'genes'}
##' if they share gene IDs after
##' annotation into genes (correlations are computed between matrix \code{SByGene}).
##' @param fileNodeDescr File where node descriptions are saved (useful when the
##' user wants to visualize the graph using Cytoscape).
##' @param fileDataGraph File where graph description is saved (useful when the
##' user wants to visualize the graph using Cytoscape).
##' @param plot if \code{TRUE} (default) plot the correlation graph
##' @param title title of the graph
##' 
##' @param col vector of colors indexed by elements of labAn; if missing, colors
##' will be automatically attributed
##' @param cutoff_graph the cutoff used to select pairs that will be included in
##' the graph
##' @param useMax if \code{TRUE}, the graph is restricted to edges that correspond to
##' maximum correlation between components, see details
##' @param tkplot If TRUE, performs interactive plot with function \code{tkplot}, else uses \code{plot.igraph}
##' @return A list consisting of \describe{
##' \item{dataGraph:}{a data.frame defining the correlation
##' graph}
##' \item{nodeAttrs:}{a data.frame describing the node
##' of the graph,}
##' \item{graph:}{the graph as an object of class \code{igraph},}
##' \item{graphid}{the id of the graph plotted with \code{tkplot}}. }
##' @export
##' @author Anne Biton
##' @seealso \code{\link{compareAn2graphfile}}, \code{\link{compareAn}}, \code{\link{cor2An}}, \code{\link{plotCorGraph}}
##' @examples
##'
##' dat1 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat1) <- paste("g", 1:1000, sep="")
##' colnames(dat1) <- paste("s", 1:10, sep="")
##' dat2 <- data.frame(matrix(rnorm(10000),ncol=10,nrow=1000))
##' rownames(dat2) <- paste("g", 1:1000, sep="")
##' colnames(dat2) <- paste("s", 1:10, sep="")
##' 
##' ## run ICA
##' resJade1 <- runICA(X=dat1, nbComp=3, method = "JADE")
##' resJade2 <- runICA(X=dat2, nbComp=3, method = "JADE")
##' 
##' ## build params
##' params <- buildMineICAParams(resPath="toy/")
##'
##' ## build IcaSet objects
##' icaSettoy1 <- buildIcaSet(params=params, A=data.frame(resJade1$A), S=data.frame(resJade1$S),
##'                           dat=dat1, alreadyAnnot=TRUE)$icaSet
##' icaSettoy2 <- buildIcaSet(params=params, A=data.frame(resJade2$A), S=data.frame(resJade2$S),
##'                           dat=dat2, alreadyAnnot=TRUE)$icaSet
##'
##' ## compare IcaSet objects
##' ## use tkplot=TRUE to get an interactive graph
##' rescomp <- runCompareIcaSets(icaSets=list(icaSettoy1, icaSettoy2), labAn=c("toy1","toy2"),
##'                              type.corr="pearson", level="genes", tkplot=FALSE)
##'
##'
##' \dontrun{
##' ## load the microarray-based gene expression datasets
##' ## of breast tumors 
##' library(breastCancerMAINZ)
##' library(breastCancerVDX)
##' data(mainz)
##' data(vdx)
##' 
##' ## Define a function used to build two examples of IcaSet objects
##' ## and annotate the probe sets into gene Symbols
##' treat <- function(es, annot="hgu133a.db") {
##'    es <- selectFeatures_IQR(es,10000)
##'    exprs(es) <- t(apply(exprs(es),1,scale,scale=FALSE))
##'    colnames(exprs(es)) <- sampleNames(es)
##'    resJade <- runICA(X=exprs(es), nbComp=10, method = "JADE", maxit=10000)
##'    resBuild <- buildIcaSet(params=buildMineICAParams(), A=data.frame(resJade$A), S=data.frame(resJade$S),
##'                         dat=exprs(es), pData=pData(es), refSamples=character(0),
##'                         annotation=annot, typeID= typeIDmainz,
##'                         chipManu = "affymetrix", mart=mart)
##'    icaSet <- resBuild$icaSet
##' }
##' ## Build the two IcaSet objects
##' icaSetMainz <- treat(mainz)
##' icaSetVdx <- treat(vdx)
##' 
##' ## compare the IcaSets
##' runCompareIcaSets(icaSets=list(icaSetMainz, icaSetVdx), labAn=c("Mainz","Vdx"), type.corr="pearson", level="genes")
##' }
runCompareIcaSets <- function(icaSets,
                              labAn,
                              type.corr = c("pearson","spearman"),
                              cutoff_zval=0,
                              level = c("genes","features","samples"),
                              fileNodeDescr = NULL,
                              fileDataGraph = NULL,
                              plot = TRUE,
                              title = "",
                              col,
                              cutoff_graph = NULL,
                              useMax = TRUE,
                              tkplot = FALSE
                              ) {

    nb <- NULL

    if (missing(labAn))
        labAn <- paste("An",1:length(icaSets))
    else
        if (length(labAn) != length(icaSets))
            stop("Length of 'labAn' is different from length of 'icaSets'.")
    
    if (!missing(col))
        if (length(col) != length(icaSets))
            stop("Length of 'col' is different from length of 'icaSets'.")
        else
            if (is.null(names(col)) | length(intersect(names(col), labAn)) != length(labAn))
                names(col) <- labAn
    
    type.corr <- type.corr[1]
    
    resCompareAn <- compareAn(icaSets = icaSets, level = level,type.corr= type.corr, labAn = labAn, cutoff_zval=cutoff_zval)
    dataGraph <- compareAn2graphfile(listPairCor=resCompareAn, useMax = useMax, file = fileDataGraph, cutoff = cutoff_graph)
    nbComp <- unlist(lapply(icaSets,function(x) length(indComp(x))))
    nbAn <- length(icaSets)
    icaSet <- NULL
    labComp <- (foreach(icaSet = icaSets, nb = nbComp, an = 1:nbAn) %do% {paste(rep("comp",sum(nb)),indComp(icaSet),sep = "")})

    nodeDescr <- nodeAttrs(nbAn = nbAn, nbComp = nbComp, labComp = labComp, labAn = labAn, file = fileNodeDescr)

    if (plot) {
        res <- plotCorGraph(title = title, dataGraph = dataGraph, nodeName = "indComp", nodeAttrs = nodeDescr, edgeWeight = "cor", col = if(missing(col)) NULL else col, nodeShape = "labAn", reciproCol = "reciprocal", tkplot = tkplot)
    return(list(dataGraph=res$dataGraph, 
                nodeAttrs=res$nodeAttrs,
                graph=res$graph, 
                graphid=res$graphid 
                ))

    }
    else 
        return(list(dataGraph=dataGraph, nodeAttrs = nodeAttrs))
}


