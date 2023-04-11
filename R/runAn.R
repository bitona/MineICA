##' This function runs the fastICA algorithm several times with random initializations.
##' The obtained components are clustered and
##' the medoids of these clusters are used as the final estimates. The returned estimates are ordered by decreasing Iq values which measure the compactness of the clusters (see details).
##' 
##'
##' This function implements in R fastICA iterations followed by a clustering step, as defined in the matlab package 'icasso'.
##' Among the indices computed by icasso,
##' only the Iq index is currently computed. As defined in 'icasso', the Iq index measures the difference between the intra-cluster similarity and the extra-cluster similiarity.
##' No visualization of the clusters is yet available.
##'
##' If \code{bootstrap=TRUE} a bootstrap (applied to the observations) is used to perturb the data before each iteration, then function \code{fastICA} is applied with random initializations.
##' 
##' By default, in 'icasso', agglomerative hierarchical clustering with average linkage is performed. To use the same clustering, please use \code{funClus="hclust"} and \code{method="average"}. But this function also allows you to apply the clustering of your choice among \code{kmeans, pam, hclust, agnes} by specifying \code{funClus} and adding the adequat additional parameters.   
##'
##' See details of the functions \code{\link[fastICA]{fastICA}}.
##' @title Run of fastICA and JADE algorithms 
##' @param X A data matrix with n rows representing observations (e.g genes) and p columns representing variables (e.g samples).
##' @param nbComp The number of components to be extracted.
##' @param nbIt The number of iterations of FastICA
##' @param alg.type If \code{alg.type="parallel"} the components are extracted simultaneously (the default), if \code{alg.type="deflation"} the components are extracted one at a time, see \code{\link[fastICA]{fastICA}}.
##' @param fun The functional form of the G function used in the approximation to neg-entropy (see 'details' of the help of function \code{\link[fastICA]{fastICA}}).
##' @param row.norm a logical value indicating whether rows of the data matrix \code{X} should be standardized beforehand (see help of function \code{fastICA})
##' @param maxit The maximum number of iterations to perform.
##' @param tol A positive scalar giving the tolerance at which the un-mixing matrix is considered to have converged.
##' @param funClus The clustering function to be used to cluster the estimates
##' @param bootstrap if TRUE the data is bootstraped before each fastICA iteration, else (default) only random initializations are done
##' @param ... Additional parameters for code{funClus}
##' @return A list consisting of: \describe{\item{A}{the estimated mixing matrix} \item{S}{the estimated source matrix}, item{W}{the estimated unmixing matrix}, \item{Iq}{Iq indices.}}
##' @author Anne Biton
##' @export
##' @examples
##' ## generate a data 
##' set.seed(2004); 
##' M <- matrix(rnorm(5000*6,sd=0.3),ncol=10)
##' M[1:100,1:3] <- M[1:100,1:3] + 2
##' M[1:200,1:3] <- M[1:200,4:6] +1 
##'   
##' ## Random initializations are used for each iteration of FastICA
##' ## Estimates are clustered using hierarchical clustering with average linkage
##' res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
##'                           nbIt=3, funClus="hclust", method="average")
##'
##' ## Data are boostraped before each iteration and random initializations
##' ## are used for each iteration of FastICA
##' ## Estimates are clustered using hierarchical clustering with ward
##' res <- clusterFastICARuns(X=M, nbComp=2, alg.type="deflation",
##'                           nbIt=3, funClus="hclust", method="ward")
##' 
##' 
clusterFastICARuns <- function(X, nbComp, nbIt=100, alg.type = c("deflation", "parallel"), fun = c("logcosh","exp"), maxit = 500, tol = 10^-6, funClus= c("hclust","agnes","pam","kmeans"), row.norm = FALSE, bootstrap=FALSE, ...) {

    message(paste("FastICA iteration 1"))
    resICA <- fastICA(X = X, n.comp = nbComp, alg.typ = alg.type, fun = fun, maxit = maxit, tol = tol, row.norm=row.norm)

    ## whitening matrix
    whit <- resICA$K
    ## dewhitening matrix 
    dewhit <- solve(t(resICA$K)%*%resICA$K) %*% t(resICA$K)

    it <-  clus <- NULL
    ## run fastICA x times after data are bootstrapped (at the gene level)
    allW <- 
        foreach(it=2:nbIt, .combine=cbind) %dopar% {
            message(paste("FastICA iteration "),it)
            if (bootstrap)
                Xbis <- X[sample(1:ncol(X), replace=TRUE),]
            else
                Xbis <- X
            res <- fastICA(X = Xbis, n.comp = nbComp, alg.typ = alg.type, fun = fun, maxit = maxit, tol = tol, row.norm=row.norm)
            res$W
        }

    
    allW <- cbind(resICA$W,allW)
    ## project W in original space (dewhitening)
    allWdewith <- t(allW)%*%dewhit
    allWdewith <- (apply(allWdewith,2,scale))

    
    ## compute similarity between W = absolute correlation values
    sim <- abs(cor(t(allWdewith)))
    dsim <- 1-sim
    centrotypes <- c()

    switch(funClus,
           hclust={
               resClus <- hclust(d=as.dist(dsim), ...)
               partition <- cutree(resClus, k=nbComp)
           },
           agnes={
               resClus <- agnes(x=dsim, diss=TRUE, ...)
               partition <-  cutree(as.hclust(resClus), k=nbComp)
           },
           pam={
               resClus <- pam(x=dsim, diss=TRUE, k=nbComp, keep.diss=FALSE, ...)
               partition <- resClus$clustering
               centrotypes <- resClus$medoids
           },
           kmeans={
               resClus <- kmeans(x=allW, centers=nbComp, ...)
               partition <- resClus$cluster
           }
           )

    ## compute Iq indices and extract centrotypes
    # Iq=avg(intra-cluster similarity) - avg(extra-cluster similarity)
    
    getIqCentr <-  function(clus, partition, sim, funClus) {
	      	indC <- which(partition==clus)
            	if (length(indC)>1) {
                   if (funClus != "pam")
                      centrotypes <- indC[which.max(apply(sim[indC,indC],1,sum))]
                   internalSim <- mean(sim[indC,indC])
            	} else {
                   if (funClus != "pam")
                      centrotypes <- indC
                      internalSim <- sim[indC,indC]                                
            	}
            	externalSim <- mean(sim[indC,setdiff(1:ncol(sim),indC)])
            	iq <- internalSim-externalSim

            	return(c(centrotype=centrotypes,iq=iq))
    }
		

    if (requireNamespace("future", quietly = TRUE) & requireNamespace("future.apply", quietly = TRUE)) {
       future::plan(future::multisession) ## => parallelize on your local computer
       iqcentr <-
       	  future.apply::future_lapply(unique(partition), getIqCentr, partition=partition, sim=sim, funClus=funClus)
    } else {
        iqcentr <-
       	  lapply(unique(partition), getIqCentr, partition=partition, sim=sim, funClus=funClus)    
    }
    
    iqcentr <- do.call(rbind, iqcentr)
    Iq <- iqcentr[,"iq"]
    if (funClus != "pam") centrotypes <- iqcentr[,"centrotype"]
    
	

#    Iq <- 
#        foreach(clus=unique(partition), .combine=c) %dopar% {
#            indC <- which(partition==clus)
#            if (length(indC)>1) {
#                if (funClus != "pam")
#                    centrotypes <- c(centrotypes,indC[which.max(apply(sim[indC,indC],1,sum))])
#                internalSim <- mean(sim[indC,indC])
#            }
#            else {
#                if (funClus != "pam")
#                    centrotypes <- c(centrotypes,indC)
#                internalSim <- sim[indC,indC]                                
#            }
#            externalSim <- mean(sim[indC,setdiff(1:ncol(sim),indC)])
#            iq <- internalSim-externalSim
#
#            return(iq)
#        }

    
    ## Extract W including the centrotypes of each cluster
    W <- whit %*% allW[,centrotypes]
    A <- solve(t(W)%*%W) %*% t(W)
    S <- X%*%W
    rownames(S) <- rownames(X)
    colnames(A) <- colnames(X)

    orderIq <- order(Iq, decreasing=TRUE)
    return(list(A=t(A)[,orderIq],S=S[,orderIq],W=W[,orderIq],Iq=Iq[orderIq]))

}


##' This function performs ICA decomposition of a matrix using functions \code{\link[fastICA]{fastICA}} and \code{\link[JADE]{JADE}}.
##'
##' See details of the functions \code{\link[fastICA]{fastICA}} and \code{\link[JADE]{JADE}}.
##' @title Run of fastICA and JADE algorithms 
##' @param method The ICA method to use, either "JADE" (the default) or "fastICA".
##' @param X A data matrix with n rows representing observations (e.g genes) and p columns representing variables (e.g samples).
##' @param nbComp The number of components to be extracted.
##' @param alg.type If \code{alg.type="parallel"} the components are extracted simultaneously (the default), if \code{alg.type="deflation"} the components are extracted one at a time, see \code{\link[fastICA]{fastICA}}.
##' @param fun The functional form of the G function used in the approximation to neg-entropy (see 'details' of the help of function \code{\link[fastICA]{fastICA}}).
##' @param maxit The maximum number of iterations to perform.
##' @param tol A positive scalar giving the tolerance at which the un-mixing matrix is considered to have converged.
##' @param ... Additional parameters for \code{fastICA} and \code{JADE}
##' @return A list, see outputs of \code{\link[fastICA]{fastICA}} and \code{\link[JADE]{JADE}}. This list includes at least three elements: \describe{\item{A}{the estimated mixing matrix} \item{S}{the estimated source matrix}, item{W}{the estimated unmixing matrix}}
##' @author Anne Biton
##' @export
##' @examples 
##' set.seed(2004); 
##' M <- matrix(rnorm(5000*6,sd=0.3),ncol=10)
##' M[1:10,1:3] <- M[1:10,1:3] + 2
##' M[1:100,1:3] <- M[1:100,1:3] +1 
##' resJade <- runICA(X=M, nbComp=2, method = "JADE", maxit=10000) 
##' 
runICA <- function(method = c("fastICA","JADE"), X, nbComp, alg.type = c("deflation", "parallel"), fun = c("logcosh","exp"), maxit = 500, tol = 10^-6, ...) {

    method <- match.arg(method)
    Xc <- as.data.frame(t(apply(X,1,function(x) x-mean(x))), stringsAsFactors = FALSE, row.names = rownames(X))
    colnames(Xc) <- colnames(X)

    if (method == "JADE") {
        resICA <- JADE(X, n.comp = nbComp, maxiter = maxit, eps = tol, ...)
        rownames(resICA$A) <- colnames(X)
        rownames(resICA$S) <- rownames(X)
    }
    else if (method == "fastICA") {
        resICA <- fastICA(X = X, n.comp = nbComp, alg.typ = alg.type, fun = fun, maxit = maxit, tol = tol, verbose = TRUE, ...)
        colnames(resICA$A) <- colnames(X)
        rownames(resICA$S) <- rownames(X)
    }

    return(resICA)

}



##' This function runs the analysis of an ICA decomposition contained in an IcaSet object, according to the parameters entered by the user and contained in a MineICAParams.
##'
##' This function calls functions of the MineICA package depending on the arguments: \describe{
##' \item{\code{\link{writeProjByComp}} (if \code{writeGenesByComp=TRUE} or \code{writeFeaturesByComp})}{which writes in html files the description of the
##' features/genes contributing to each component, and their projection values on all the components.}
##'  \item{\code{\link{plot_heatmapsOnSel}} (if \code{plotHeatmap=TRUE})}{which plots heatmaps of the data restricted to the contributing features/genes of each component.}
##'  \item{\code{\link{plotPosAnnotInComp}} (if \code{plotHist=TRUE})}{which plots, within the histogram of the sample contribution values of every component, the position of groups of samples formed according to the sample annotations contained in \code{pData(icaSet)}.}
##'  \item{\code{\link{clusterSamplesByComp}} (if \code{runClustering=TRUE})}{which clusters the samples according to each component.}
##'  \item{\code{\link{clusVarAnalysis}} (if \code{runClustering=TRUE})}{which computes the chi-squared test of association between a given clustering of the samples and each annotation level contained in \code{pData(icaSet)}, and summarizes the results in an HTML file. }
##'  \item{\code{\link{runEnrich}} (if \code{runGOstats=TRUE})}{which perforns enrichment analysis of the contributing genes of the components using package \link{GOstats}.}
##'  \item{\code{\link{qualVarAnalysis}} and \code{\link{quantVarAnalysis}} (if \code{varAnalysis=TRUE})}{which tests if the groups of samples formed according to sample annotations contained in \code{pData(icaSet)} are differently distributed on the components, in terms of contribution value. }
##' }
##'
##' Several directories containing the results of each analysis are created by the function:
##' \describe{
##' \item{ProjByComp:}{contains the annotations of the features or genes, one file per component;}
##' \item{varAnalysisOnA:}{contains two directories: 'qual/' and 'quant/' which respectively contain the results of the association between components qualitative and quantitative variables;}
##' \item{Heatmaps:}{contains the heatmaps (one pdf file per component) of contributing genes by component;}
##' \item{varOnSampleHist:}{contains athe histograms of sample contributions superimposed with the histograms of the samples grouped by variable;}
##' \item{cluster2var:}{contains the association between a clustering of the samples performed on the mixing matrix \code{A} and the variables.} 
##' }
##' @title Run analysis of an IcaSet object
##' @param params An object of class \code{\link[MineICA:MineICAParams-class]{MineICAParams}} containing the parameters of the analysis.
##' @param icaSet An object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}.
##' @param keepVar  The variable labels to be considered, i.e a subset of the annotation variables available in  (\code{varLabels(icaSet)}).
##' @param keepSamples  The samples to be considered, i.e a subset of  (\code{sampleNames(icaSet)}).
##' @param heatmapCutoff The cutoff (applied to the scaled feature/gene projections contained in S/SByGene) used to select the contributing features/genes.
##' @param funClus The function to be used to cluster the samples, must be one of \code{c("Mclust","kmeans","pam","pamk","hclust","agnes")}. Default is \code{"Mclust"}.
##' @param nbClus The number of clusters to be computed when applying \code{funClus}. Can be missing (default) if \code{funClus="Mclust"} or \code{funClus="pamk"}.
##' @param keepComp The indices of the components to be analyzed, must be included in \code{indComp(icaSet)}. If missing, all components are treated.
##' @param adjustBy The way the p-values of the Wilcoxon and Kruskal-Wallis tests should be corrected for multiple testing: \code{"none"} if no p-value correction has to be done, \code{"component"} if the p-values have to be corrected by component, \code{"annotation"} if the p-values have to be corrected by variable
##' @param typePlot The type of plot used to show distribution of sample-groups contributions, either "density" or "boxplot"
##' @param mart A mart object used for annotation, see function \code{\link[biomaRt]{useMart}} 
##' @param dbGOstats The used database to use ('GO' and/or 'KEGG'), default is both.
##' @param ontoGOstats A string specifying the GO ontology to use. Must be one of 'BP', 'CC', or 'MF', see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}. Only used when argument \code{dbGOstats} is 'GO'.
##' @param condGOstats A logical indicating whether the calculation should conditioned on the GO structure, see \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}.
##' @param cutoffGOstats The p-value threshold used for selecting enriched gene sets, default is params["pvalCutoff"]
##' @param writeGenesByComp If TRUE (default) the gene projections (\code{SByGene(icaSet)}) are written in an html file and annotated using \code{biomaRt} for each component.
##' @param writeFeaturesByComp  If TRUE (default) the feature projections (\code{S(icaSet)}) are written in an html file and annotated using \code{biomaRt}  for each component.
##' @param runGOstats If TRUE the enrichment analysis of the contributing genes is run for each component using package \code{GOstats} (default is TRUE). 
##' @param plotHist If TRUE the position of the sample annotations within the histograms of the sample contributions are plotted. 
##' @param plotHeatmap If TRUE the heatmap of the contributing features/genes are plotted for each component. 
##' @param runClustering If TRUE the potential associations between a clustering of the samples (performed according to the components), and the sample annotations, are tested using chi-squared tests.
##' @param runVarAnalysis If TRUE the potential associations between sample contributions (contained in \code{A(icaSet)}) are tested using Wilcoxon or Kruskal-Wallis tests.
##' @param onlySign If TRUE (default), only the significant results are plotted in functions \code{qualVarAnalysis, quantVarAnalysis, clusVarAnalysis}, else all plots are done.
##' @param selCutoffWrite The cutoff  applied to the absolute feature/gene projection values to select the features/genes that will be annotated using package \code{biomaRt}, default is 2.5.
##' @param clusterOn Specifies the matrix used to apply clustering if \code{runClustering=TRUE}: \describe{
##' \item{\code{"A"}:}{the clustering is performed in one dimension, on the vector of sample contributions,}
##' \item{"S":}{the clustering is performed on the original data restricted to the contributing individuals,}
##' \item{"AS":}{the clustering is performed on the matrix formed by the product of the column of A and the row of S.}}
##' @return NULL
##' @seealso \code{\link{writeProjByComp}}, 
##' @author Anne Biton
##' @export
##' @examples \dontrun{
##'
##' ## load an example of IcaSet
##' data(icaSetCarbayo)
##' ## make sure the 'mart' attribute is correctly defined
##' mart(icaSetCarbayo) <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
##' 
##' ## creation of an object of class MineICAParams
##' ## here we use a low threshold because 'icaSetCarbayo' is already
##' # restricted to the contributing features/genes 
##' params <- buildMineICAParams(resPath="~/resMineICACarbayotestRunAn/", selCutoff=2, pvalCutoff=0.05)
##' require(hgu133a.db)
##' 
##' runAn(params=params, icaSet=icaSetCarbayo)
##' }
runAn <- function(params,
                  icaSet,
                  keepVar,
                  heatmapCutoff = params["selCutoff"],
                  funClus = c("Mclust","kmeans"),
                  nbClus,
                  clusterOn = "A",
                  keepComp,
                  keepSamples,
                  adjustBy = c("none","component","variable"),
                  typePlot = c("boxplot","density"), #or density
                  mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"),
                  dbGOstats = c("KEGG","GO"),
                  ontoGOstats = "BP",
                  condGOstats = TRUE,
                  cutoffGOstats = params["pvalCutoff"],
                  writeGenesByComp = TRUE,
                  writeFeaturesByComp = FALSE ,
                  selCutoffWrite = 2.5,
                  runVarAnalysis=TRUE,
                  onlySign=T,
                  runClustering=FALSE,
                  runGOstats = TRUE ,
                  plotHist = TRUE ,
                  plotHeatmap = TRUE
      
                  
                  ) {

            
    if (missing(keepVar))
        keepVar <- varLabels(icaSet)
    else {
        difa <- setdiff(keepVar,varLabels(icaSet))
        if (length(difa)>0) {
            if (length(difa)==length(keepVar))
                stop("The variables contained in arg 'keepVar' are not available in 'varLabels(icaSet)'")
            else
                warning(paste("The variables",difa,"are not available in'varLabels(icaSet)'"))
        }
        keepVar <- intersect(keepVar,varLabels(icaSet))
    }
    

    if (ncol(pData(icaSet))>0) {

        x <- NULL
        quantVar <- which(unlist(foreach(x=pData(icaSet)) %do% is.numeric(x)))
        
        ind <- which(apply(as.matrix(pData(icaSet)[,quantVar]),2,function(x) length(unique(x)))<2)
        
        if(length(ind)>0)
            quantVar <- quantVar[-ind]
        
        quantVar <- intersect(keepVar,varLabels(icaSet)[quantVar])
        qualVar <- setdiff(keepVar,quantVar)

        nbLevByQual <- sapply(qualVar, function(v,annot) {length(unique(as.character(annot[,v])))}, annot=pData(icaSet))
        indDel <- which(nbLevByQual %in% c(length(sampleNames(icaSet)),length(sampleNames(icaSet))-1,1))
        if (length(indDel)>0)
            qualVar <- qualVar[-indDel]
        
        if (length(annot2col(params))==0)
            colours <- annot2Color(pData(icaSet)[,qualVar,drop=FALSE])
        else
            colours <- annot2col(params)


    }
    else {
        quantVar <- NULL
        qualVar <- NULL
    }

    if (length(qualVar)==0)
        message("No qualitative variable is found in pheno data.")
    else
        message(paste("The detected qualitative variables are:",paste(qualVar,collapse=", ")))
    
    if (length(quantVar)==0)
        message("No quantitative variable is found in pheno data.")
    else
        message(paste("The detected quantitative variables are",paste(quantVar,collapse=", ")))
        


    if (missing(keepComp))
        keepComp <- indComp(icaSet)
    
    if (missing(keepSamples))
        keepSamples <- sampleNames(icaSet)

    if ((length(selCutoff(params))>1) & (length(selCutoff(params))<nbComp(icaSet))) {
        warning("The length of 'selCutoff' must be either 1 or equal to the number of components. Only the first value is used.")
        params["selCutoff"] <- selCutoff(params)[1]
    }

    
    icaSet <- icaSet[,keepSamples, keepComp]

    funClus <- match.arg(funClus)
    typePlot <- match.arg(typePlot)
    adjustBy <- match.arg(adjustBy)

    
    if (writeGenesByComp) {
        message("..Write gene contributions..")
        resAnnotGenes <-
            writeProjByComp(icaSet = icaSet,
                            params = params,
                            mart = mart, 
                            typeRetrieved = if (length(typeID(icaSet)>0)) {if (typeID(icaSet)["geneID_biomart"] == "hgnc_symbol") NULL} else "hgnc_symbol", 
                            addNbOcc = TRUE,
                            selectionByComp = NULL,
                            level = if (nrow(SByGene(icaSet))>0) "genes" else "features",
                            typeId = typeID(icaSet)["geneID_biomart"],
                            selCutoffWrite = selCutoffWrite
                            )
    }
    if (writeFeaturesByComp) {    
        message("..Write feature contributions..")
        resAnnotFeatures <-
            writeProjByComp(icaSet = icaSet,
                            params = params,
                            mart = mart, 
                            typeRetrieved = if (length(typeID(icaSet)>0)) {if (typeID(icaSet)["featureID_biomart"]  == "hgnc_symbol") NULL} else "hgnc_symbol", 
                            addNbOcc = TRUE,
                            selectionByComp = NULL,
                            level = if (nrow(SByGene(icaSet))>0) "genes" else "features",
                            typeId = typeID(icaSet)["featureID_biomart"]
                            )
    }
    
    if (plotHeatmap) {
        message("..Plot heatmaps..")
        
        heatmapsPath <- paste(params["resPath"], "Heatmaps/",sep="");
        heatmapsPathWithDendro <- paste(heatmapsPath, "withDendro/",sep="");    heatmapsPathWithoutDendro <- paste(heatmapsPath, "withoutDendro/",sep="");
        system(paste("mkdir", heatmapsPath, heatmapsPathWithoutDendro, heatmapsPathWithDendro), ignore.stderr=TRUE)
        
        plot_heatmapsOnSel(icaSet = icaSet,
                           selCutoff = heatmapCutoff,
                           level = "genes",
                           keepVar = if(!is.null(qualVar)) qualVar,
                           doSamplesDendro = TRUE,
                           doGenesDendro = TRUE,
                           heatmapCol = maPalette(low = "blue",high = "red", mid = "yellow", k=44),
                           file = paste(heatmapsPathWithDendro,"withDendro_",sep=""),
                           annot2col = annot2col(params))
        
        plot_heatmapsOnSel(icaSet = icaSet,
                           selCutoff = heatmapCutoff,
                           level = "genes",
                           keepVar = if(!is.null(qualVar)) qualVar, 
                           doSamplesDendro = FALSE,
                           doGenesDendro = FALSE,
                           heatmapCol = maPalette(low = "blue",high = "red", mid = "yellow", k=44),
                           file = paste(heatmapsPathWithoutDendro,"orderedByContrib_",sep=""),
                           annot2col=annot2col(params))
    }

    if (runVarAnalysis) {    
        dir_var = paste(resPath(params), "varAnalysisOnA",if (adjustBy=="none") "_noAdjust" else paste("_adjustBy",adjustBy,sep=""),"/",sep="")
        system(paste("mkdir", dir_var),ignore.stderr=TRUE)
        
        
        if (length(qualVar)>0) {
            dir_qualVar = paste(dir_var, "qual/", sep = "")
            system(paste("mkdir", dir_qualVar),ignore.stderr=TRUE)
            resQualVar <- qualVarAnalysis(params=params, icaSet = icaSet, keepVar = qualVar, adjustBy = adjustBy, path = dir_qualVar, doPlot = TRUE, colours = colours, cutoff = pvalCutoff(params),  typeImage = "png", filename =  "qualVar", typePlot = typePlot, onlySign=onlySign)
            
        }
        
        if (length(quantVar)>0) {
            dir_quantVar = paste(dir_var, "quant/", sep = "")
            system(paste("mkdir", dir_quantVar),ignore.stderr=TRUE)
            resQuantvar <- quantVarAnalysis(params=params, icaSet = icaSet, keepVar = quantVar, adjustBy = adjustBy, path = dir_quantVar, doPlot = TRUE, colours = colours, cutoff = 0.3, cutoffOn="cor",  typeImage = "png", filename =  "quantVar", onlySign=onlySign)
            
        }
        
        
    }
    
    
    
    if (plotHist) {
        if (length(qualVar)>0) {
            
            message("..Plot histograms of sample contributions..")
            dir_Hist = paste(resPath(params), "varOnSampleHist/", sep = "")
            dir_Histcomp2annot = paste(dir_Hist, "comp2annot/", sep = "")
            dir_Histannot2comp = paste(dir_Hist, "annot2comp/", sep = "")
            system(paste("mkdir", dir_Hist, dir_Histcomp2annot, dir_Histannot2comp),ignore.stderr=TRUE)
            
            plotPosAnnotInComp(icaSet = icaSet, breaks = 15, keepVar = qualVar, colAll = "grey74", colSel = colours,funClus = funClus, pathPlot = dir_Histcomp2annot, nbClus = nbClus, by="comp")
            
            plotPosAnnotInComp(icaSet = icaSet, breaks = 15, keepVar = qualVar, colAll = "grey74", colSel = colours,funClus = funClus, pathPlot = dir_Histannot2comp, nbClus = nbClus, by="annot")
            
        }
    }
    
    
    if (runClustering) {
        message("..Clustering on components and test association of clusters with variables..")
        
        dir_cluster2annot = paste(resPath(params), "cluster2var/", sep = "")
        system(paste("mkdir", dir_cluster2annot),ignore.stderr=TRUE)
        resClus <-clusterSamplesByComp(params=params, icaSet=icaSet, funClus=funClus, clusterOn=clusterOn, filename=paste(dir_cluster2annot,"clusters_",funClus,"On",clusterOn,sep=""), nbClus=nbClus)
        
        
        if (length(qualVar)>0) {
            
            resClusTests <- clusVarAnalysis(icaSet = if ("icaSet" %in% names(resClus)) resClus$icaSet else icaSet, params = params, keepVar = qualVar, resClus = resClus$clus, funClus = funClus, adjustBy = adjustBy, doPlot = TRUE, path = dir_cluster2annot, filename=paste("clusVar_",funClus,"On",clusterOn,sep=""), onlySign=onlySign)#, if (!missing(keepComp))  keepComp = keepComp else ...)
            
        }
        
        
    }
    
    
    if (runGOstats) {
        resEnrich <- runEnrich(icaSet = icaSet, params = params, dbs = dbGOstats, ontos = ontoGOstats, cond = condGOstats, hgCutoff = cutoffGOstats) 
    }    
    
    return(NULL)
    
}


##' This function tests the enrichment of the components of an \code{\link{IcaSet}} object using package \code{GOstats} through function \code{hyperGTest}.   
##'
##'
##' 
##' An annotation package should be available in \code{annotation(icaSet)} to provide the contents of the gene sets. If none corresponds to the technology you deal with, please choose the org.*.eg.db package according to the organism (for example org.Hs.eg.db for Homo sapiens). By default, if \code{annotation(icaSet)} is empty and organism is one of \code{c("Human","HomoSapiens","Mouse","Mus Musculus")}, then either \code{org.Hs.eg.db} or \code{org.Mm.eg.db} is used.
##' 
##' Use of \code{GOstats} requires the input IDs to be Entrez Gene, this function will therefore annotate either the feature names or the gene names into Entrez Gene ID using either the annotation package (\code{annotation(icaSet)}) or \code{biomaRt}.
##' 
##' Three types of enrichment tests are computed for each component: the threshold is first used to select gene based on their absolute projections, then positive and negative projections are treated individually.
##' 
##' For each database \code{db} (each ontology if \code{db} is "GO"), this function writes an HTML file containing the outputs of the enrichment tests computed through the function \code{\link[Category]{hyperGTest}}. The corresponding files are located in \code{resPath(icaSet)}/GOstatsEnrichAnalysis/byDb/.
##' The results obtained for each database/ontology are then merged into an array for each component, this array is written as an HTML file in the directory \code{resPath(icaSet)}/GOstatsEnrichmentAnalysis/ (this directory is first deleted if it already exists). This file is the one the user should look at.
##'
##' The outputs of \code{\link[Category]{hyperGTest}} that are given in each table are: 
##' \describe{
##' \item{DB, ID, Term:}{the database, the gene set ID, and the gene Set name}
##' \item{P-value:}{probability of observing the number of genes annotated for the gene set among the selected gene list, knowing the total number of annotated genes among the universe}, 
##' \item{Expected counts:}{expected number of genes in the selected gene list to be found at each tested category term/gene set,}
##' \item{Odds ratio:}{odds ratio for each category term tested which is an indicator of the level of enrichment of genes within the list as against the universe,} 
##' \item{Counts:}{number of genes in the selected gene list that are annotated for the gene set,}
##' \item{Size:}{number of genes from the universe annotated for the gene set.}}
##' 
##' @title Enrichment analysis through \link{GOstats}
##' @param icaSet An object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}
##' @param params An object of class \code{\link[MineICA:MineICAParams-class]{MineICAParams}} providing the parameters of the analysis
##' @param dbs The database to use, default is \code{c("GO","KEGG")}
##' @param ontos A string specifying the GO ontology to use. Must be one of \code{"BP"}, \code{"CC"}, or \code{"MF"}, see \code{\link[Category]{GOHyperGParams-class}}. Only used when argument \code{dbs} includes \code{"GO"}.
##' @param cond A logical indicating whether the calculation should condition on the GO structure, see \code{\link[Category]{GOHyperGParams-class}}. Only used when argument \code{dbs} includes \code{"GO"}.
##' @param hgCutoff The threshold p-value for statistical significance, default is \code{pvalCutoff(params)}
##' @return NULL
##' @author Anne Biton
##' @examples \dontrun{
##' # Load examples of IcaSet object
##' data(icaSetCarbayo)
##' 
##' ## Define parameters
##' # Use threshold 3 to select contributing genes on which enrichment analysis will be applied
##' # Results of enrichment analysis will be written in path 'resPath(params)/GOstatsEnrichAnalysis'
##' params <- buildMineICAParams(resPath="carbayo/", selCutoff=3)
##' 
##' ## Run enrichment analysis on the first two components contained in the icaSet object 'icaSetCarbayo'
##' runEnrich(params=params,icaSet=icaSetCarbayo[,,1:2],dbs="GO", ontos="BP")
##' 
##' }
##' @export
##' @seealso \code{\link{buildIcaSet}}, \code{\link[biomaRt]{useMart}}, \code{\link[Category]{hyperGTest}}, \code{\link[Category:GOHyperGParams-class]{GOHyperGParams}}, \code{\link{hypergeoAn}}, \code{\link{mergeGostatsResults}}
##' 
runEnrich <- function(icaSet, params, dbs = c("KEGG","GO"), ontos = c("BP","CC","MF"), cond = TRUE, hgCutoff = params["pvalCutoff"]) {


    entrezgene <- hgnc_symbol <- NULL

    
    mart <- mart(icaSet)

    if (length(params["selCutoff"])==1) {
        message(paste("..Enrichment analysis using GOstats, cutoff ", params["selCutoff"], "..",sep=""))
    }
    else {
        if (all(params["selCutoff"]==params["selCutoff"][1])) {
            params["selCutoff"] <- params["selCutoff"][1] 
            message(paste("..Enrichment analysis using GOstats, cutoff ", unique(params["selCutoff"]), "..",sep=""))
        }
        else {
            if (length(params["selCutoff"])!=nbComp(icaSet)) {
                warning("The length of slot 'selCutoff' of params must be either 1 or equal to the number of components. Only the first value is used.")
                params["selCutoff"] <-  params["selCutoff"][1]
                message(paste("..Enrichment analysis using GOstats, cutoff ", params["selCutoff"], "..",sep=""))

            }
            else
                message(paste("..Enrichment analysis using GOstats, cutoff ", paste(params["selCutoff"],collapse=","), "..",sep=""))
        }
    }
        

    ## Feature selection and annotation with Entrez id for GOstats
    
    ## if annotation package is available
    if (length(annotation(icaSet)) != 0 && annotation(icaSet) != "" &&  substr(annotation(icaSet), start = 1, stop = 3) != "org") {
        
        sel2entrez <- selectContrib(Slist(icaSet), cutoff = params["selCutoff"])
        
        pack.annot.EID <- eval(as.name(paste(gsub(".db", "", annotation(icaSet)), "ENTREZID", sep = "")))
    
        ## For each component : Annote and Select negative and positive projections in two different vectors.
        sel2entrez <- llply(sel2entrez, function(compSel) {
            compSel_neg =  na.omit(unique(unlist(AnnotationDbi::mget(names(compSel[which(compSel < 0)]), pack.annot.EID, ifnotfound = NA))))
            compSel_pos =  na.omit(unique(unlist(AnnotationDbi::mget(names(compSel[which(compSel > 0)]), pack.annot.EID, ifnotfound = NA))))
            compSel = na.omit(unique(unlist(AnnotationDbi::mget(names(compSel), pack.annot.EID, ifnotfound = NA))))
            return(list(compSel = compSel, compSel_neg = compSel_neg, compSel_pos = compSel_pos))
        })

        universe <- as.character(na.omit(unique(unlist(AnnotationDbi::mget(featureNames(icaSet), pack.annot.EID, ifnotfound = NA)))))#names(Slist(icaSet)[[1]])
        entrez2hugo <- NULL

    }
    else { ## else annotation with biomaRt
        sel2entrez <- selectContrib(SlistByGene(icaSet), cutoff = params["selCutoff"])

        if (mart@dataset == "")
            stop("Please fill the mart attribute of the 'icaSet' object using the function 'useMart'.")
        
        if (!(typeID(icaSet)["geneID_biomart"] %in%  listFilters(mart)[,1]))
            stop(paste(typeID(icaSet), "is not available in listFilters(mart)."))

        message(paste("BiomaRt is used to find the correspondence between ", typeID(icaSet)["geneID_biomart"], ", Symbols, and EntrezGene IDs.", sep = ""))
        universe1 <- getBM(values = geneNames(icaSet), filters = typeID(icaSet)["geneID_biomart"], attributes = c(typeID(icaSet)["geneID_biomart"],"entrezgene", "hgnc_symbol"), mart = mart)
        sumNA <- sum(is.na(universe1$entrezgene))
        if (sumNA>0)
            message(paste("According to biomaRt, ", sumNA, " gene names have no EntrezGene ID across the ", length(featureNames(icaSet)), " available.", sep = ""))
        sumNA <- sum(is.na(universe1$hgnc_symbol))
        if (sumNA>0)
            message(paste("According to biomaRt, ", sumNA, " gene names have no Symbol ID across the ", length(featureNames(icaSet)), " available.", sep = ""))
        
        universe2 <- subset(universe1, !is.na(entrezgene))
        universe <- as.character(universe2$entrezgene)
        names(universe) <- universe2[,1]
        universe2 <- subset(universe2, !is.na(hgnc_symbol))
        universe2 <- subset(universe2, hgnc_symbol != "")
        entrez2hugo <- as.character(universe2$hgnc_symbol)
        names(entrez2hugo) <- as.character(universe2$entrezgene)
      

        
        ## For each component : Annotate and Select negative and positive projections in two different vectors.
        sel2entrez <- llply(sel2entrez, function(compSel, icaSet, mart, univ) {
            compSel <- compSel[intersect(names(univ), names(compSel))]
            compSel_neg <- univ[names(compSel[compSel < 0])]
            compSel_pos <- univ[names(compSel[compSel > 0])]
            compSel <- univ[names(compSel)]
            return(list(compSel = compSel, compSel_neg = compSel_neg, compSel_pos = compSel_pos))
        }, icaSet = icaSet, mart = mart, univ = universe)        
    }

    
    if (length(annotation(icaSet)) == 0 || annotation(icaSet)=="") {
         pack <- c(HUMAN="org.Hs.eg.db","HOMO SAPIENS"="org.Hs.eg.db","MOUSE"="org.Mm.eg.db", "MUS MUSCULUS"="org.Mm.eg.db")[toupper(organism(icaSet))]
         if (is.na(pack))
             stop ("To run GOstats analysis, you need to provide an annotation package in package(icaSet). If none corresponds to the technology you deal with, choose the org.*.eg.db package according to the organism (for example org.Hs.eg.db for Homo sapiens).")
         else {
             Biobase:::annotation(icaSet) <- pack
             message(paste("Since no annotation package is available in icaSet, package",pack,"is used."))
         }
    }
 
                
    #message("..Computing enrichment analysis with  hypergeometric test using GOstats..")
    dir_fa <- paste(params["resPath"],"GOstatsEnrichAnalysis/",sep="")
    dir_fa_byDB <- paste(params["resPath"],"GOstatsEnrichAnalysis/byDB/",sep="")
    
    system(paste("rm -r", dir_fa),ignore.stderr=TRUE)
    system(paste("mkdir", dir_fa),ignore.stderr=TRUE)
    system(paste("mkdir", dir_fa_byDB),ignore.stderr=TRUE)

    resHg <- list()
    length(resHg) <- length(dbs)
    names(resHg) <- dbs
                ### GO
                for (db in dbs) {
                 message(paste("..Enrichment Analysis on", toupper(db), "database.."))
                 dir_fa_db <- paste(dir_fa_byDB, "/", toupper(db), sep = "")
		 system(paste("mkdir",dir_fa_db),ignore.stderr=TRUE)
                 
                 ## GO BP CC MF
                 if (db == "GO") {
                     resHg[[db]] <- list()
                     length(resHg[[db]]) <- length(ontos)
                     names(resHg[[db]]) <- ontos
                   for (onto in ontos) {
                       print(onto)
                       dir_fa_dbcat <- paste(dir_fa_db, "/",toupper(onto), sep = ""); 
                       system(paste("mkdir", dir_fa_dbcat),ignore.stderr=TRUE)
                       hg.res <- hypergeoAn(
                                            icaSet = icaSet,
                                            params = params,
                                            SlistSel = sel2entrez,
                                            hgCutoff = hgCutoff,
                                            cond = cond,
                                            path = (paste(dir_fa_dbcat,"/",sep = "")), 
                                            db = db,
                                            onto = onto,
                                            universe = universe,
                                            entrez2symbol = entrez2hugo)
                       resHg[[db]][[onto]] <- hg.res
                   }
                }
                else if (db == "KEGG") {
                       hg.res <- hypergeoAn(
                                            icaSet = icaSet,
                                            params = params,
                                            hgCutoff = hgCutoff,
                                            SlistSel = sel2entrez,
                                            cond = cond,
                                            path = (paste(dir_fa_db,"/", sep = "")), 
                                            db = db,
                                            universe = universe,
                                            entrez2symbol = entrez2hugo)
                       resHg[[db]] <- hg.res
                }
                
                }

            mergeGostatsResults(resPath = params["resPath"], GOstatsPath = "GOstatsEnrichAnalysis", rdata="hgres", cutoffs = params["selCutoff"], hgCutoff = hgCutoff, cond = cond, pathGenes =  genesPath(params))

    return(resHg)


}
