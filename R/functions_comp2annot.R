##' This function allows to cluster samples according to the results of an ICA decomposition.
##' One clustering is run independently for each component.
##'
##' @title Cluster samples from an IcaSet
##' @param icaSet An \code{IcaSet} object
##' @param params A \code{MineICAParams} object
##' @param funClus The function to be used for clustering, must be one of \code{c("Mclust","kmeans","pam","pamk","hclust","agnes")}
##' @param filename A file name to write the results of the clustering in
##' @param clusterOn Specifies the matrix used to apply clustering: \describe{
##' \item{\code{"A"}:}{the clustering is performed in one dimension, on the vector of sample contributions,}
##' \item{"S":}{the clustering is performed on the original data restricted to the contributing individuals.}}
##' @param level The level of projections to be used when \code{clusterOn="S"}, either \code{"features"} or \code{"genes"}.
##' @param nbClus The number of clusters to be computed, either a single number or a numeric vector whose length equals the number of components. If missing (only allowed if \code{funClus} is one of \code{c("Mclust","pamk")})
##' @param metric Metric used in \code{pam} and \code{hclust}, default is \code{"euclidean"}
##' @param method Method of hierarchical clustering, used in \code{hclust} and \code{agnes}
##' @param ... Additional parameters required by the clustering function \code{funClus}.res <- clusterSamplesByComp(icaSet=icaSetCarbayo, params=params, funClus="kmeans",
##' @return A list consisting of three elements \describe{\item{clus:}{a list specifying the sample clustering for each component,}\item{resClus:}{the complete output of the clustering function,}\item{funClus:}{the function used to perform the clustering.}}. When \code{clusterOn="S"}, if some components were not used because no contributing elements is selected using the cutoff, the icaSet with the corresponding component deleted is also returned. 
##' @author Anne Biton
##' @seealso \code{Mclust}, \code{kmeans}, \code{pam}, \code{pamk}, \code{hclust}, \code{agnes}, \code{cutree}
##' @export
##' @examples 
##' data(icaSetCarbayo)
##' params <- buildMineICAParams(resPath="carbayo/", selCutoff=4)
##' 
##' ## cluster samples according to their contributions
##' # using Mclust without a number of clusters
##' res <- clusterSamplesByComp(icaSet=icaSetCarbayo, params=params, funClus="Mclust",
##'                             clusterOn="A", filename="clusA")
##' 
##' # using kmeans
##' res <- clusterSamplesByComp(icaSet=icaSetCarbayo, params=params, funClus="kmeans",
##'                             clusterOn="A", nbClus=2, filename="clusA")
##' 
##' 
clusterSamplesByComp <- function (icaSet, params, funClus=c("Mclust","kmeans","pam","pamk","hclust","agnes"), filename, clusterOn=c("A","S"), level=c("genes","features"), nbClus, metric="euclidean", method="ward", ...) {


    level <- match.arg(level)
    clusterOn <- match.arg(clusterOn)
    #funClus <- match.arg(funClus)
    if (!missing(nbClus)) {
        if (length(nbClus)==1)
            nbClus <- rep(nbClus,nbComp(icaSet))
         else if (length(nbClus)<nbComp(icaSet))
             stop("Length of arg 'nbClus' must either 1 or equal to the number of components.")

    }
    
    emptysel <- c()
    
    G <- centers <- k <- uncert <- NULL
    ## build data the clustering is applied to
    switch(clusterOn,
           A={
              d <- Alist(icaSet)
           },
           S={
              samplesOrder <- lapply(Alist(icaSet),names)
              d <- build_sortHeatmap(icaSet = icaSet, selCutoff = selCutoff(params), level = level, samplesOrder=samplesOrder)
              d <- lapply(d,function(x) t(x))
              emptysel <- which(unlist(lapply(d,nrow))==0)
              if (length(emptysel)>0) {
                  message(paste("The component(s) ", paste(indComp(icaSet)[emptysel],collapse=" "), " has/have no contributing elements when using cutoff", if (length(selCutoff(params))==1) selCutoff(params) else selCutoff(params)[emptysel],".",sep=""))
                  icaSet <- icaSet[,,keepComp=indComp(icaSet)[-emptysel]]
                  d <- d[-emptysel]
              }
                  
           })

    switch(funClus,
           Mclust={ if (missing(nbClus))
                        resClus <- lapply(d,funClus, G=1:9)
                    else
                        resClus <- foreach(G=nbClus,dat=d) %do% {Mclust(dat,G=G, ...)}
                },
           kmeans={  resClus <- foreach(centers=nbClus,dat=d) %do% {kmeans(dat,centers=centers, ...)}},
           pam={ resClus <- foreach(k=nbClus,dat=d) %do% {pam(dat,k=k,metric=metric,...)}},
           pamk={
                 resClus <- foreach(dat=d) %do% {
                     res <- fpc:::pamk(dat)
                     res$pamobject
                 }
             },
           hclust={resClus <- foreach(dat=d, k=nbClus) %do% {
               dd <- dist(dat, method = metric) 
               res <- hclust(dd, method=method)
               resClus <-  list(cluster=cutree(res, k=k))
           }},
           agnes={resClus <- foreach(dat=d, k=nbClus) %do% {
               dd <- dist(dat, method = metric) 
               res <- agnes(dd, method=method)
               resClus <-  list(cluster=cutree(as.hclust(res), k=k))
           }}                   
       )

    nameClus <- c(kmeans="cluster",Mclust="classification",pam="clustering",pamk="clustering",hclust="cluster",agnes="cluster")

    ### extract clustering results
    onlyClus <- lapply(resClus,function(x,nameClus,funClus) x[[nameClus[funClus]]], nameClus=nameClus,funClus=funClus)
    names(onlyClus) <- indComp(icaSet)
    names(resClus) <- indComp(icaSet)

    
    ### write clustering results
    if (!missing(filename)) {
        onlyClusw <- onlyClus

        if (funClus=="Mclust") {
            uncertClus <- lapply(resClus,function(x) signif(x$uncertainty,3))
            names(uncertClus) <- paste(compNames(icaSet),"uncert.",sep="-")
            names(onlyClusw) <- paste(compNames(icaSet),"clus.",sep="-")
            
            onlyClusw <- 
                foreach(cluster=onlyClusw,uncert=uncertClus,.combine=cbind) %do% {
                    cbind(cluster,uncert)
                }
            colnames(onlyClusw) <- paste(as.vector(sapply(compNames(icaSet),rep,2)),colnames(onlyClusw),sep="-")
            
        }
        else {
            onlyClusw <- as.data.frame(onlyClusw)
            colnames(onlyClusw) <- compNames(icaSet)
        }
        
        
        filename <- paste(resPath(params),gsub(gsub(filename,pattern=resPath(params),replacement=""),pattern=".txt",replacement=""),".txt",sep="")
        write.table(onlyClusw,file=filename, sep = "\t", quote=FALSE)
        
    }

    ## if some components were not used because no contributing elements is selected using the cutoff, return icaSet
    if (length(emptysel)>0)
        return(list(clus=onlyClus,resClus=resClus,funClus=funClus, icaSet=icaSet))

    else
        return(list(clus=onlyClus,resClus=resClus,funClus=funClus))
 
}

clScore <- function(c1,c2){
  uc1 <- unique(c1)
  uc2 <- unique(c2)
  tmp <- 0
  for(c in uc1){
    for(cc in uc2){
      tmp <- tmp + (sum((c1 == c) & (c2 == cc))^2)/(sum(c1 == c)*sum(c2 == cc))
    }
  }
  return((length(uc1)+length(uc2))/2 - tmp)
}

##' This function allows to cluster samples according to the results of an ICA decomposition.
##' Several clustering functions and several levels of data for clustering can be performed by the function.
##'
##' One clustering is run independently for each component.
##' 
##' 
##' @title Cluster samples from an IcaSet
##' @param icaSet An \code{IcaSet} object
##' @param params A \code{MineICAParams} object
##' @param funClus The function to be used for clustering, must be several of \code{c("Mclust","kmeans","pam","pamk","hclust","agnes")}
##' @param filename A file name to write the results of the clustering in
##' @param clusterOn Specifies the matrix used to apply clustering, can be several of: \describe{
##' \item{\code{"A"}:}{the clustering is performed in one dimension, on the vector of sample contributions,}
##' \item{\code{"S"}:}{the clustering is performed on the original data restricted to the contributing individuals.}}
##' @param level The level of projections to be used when \code{clusterOn="S"}, either \code{"features"} or \code{"genes"}.
##' @param nbClus The number of clusters to be computed, either a single number or a numeric vector whose length equals the number of components. If missing (only allowed if \code{funClus} is one of \code{c("Mclust","pamk")})
##' @param metric Metric used in \code{pam} and \code{hclust}, default is \code{"euclidean"}
##' @param method Method of hierarchical clustering, used in \code{hclust} and \code{agnes}
##' @param ... Additional parameters required by the clustering function \code{funClus}.
##' @return A list consisting of three elements \describe{\item{clus:}{a data.frame specifying the sample clustering for each component using the different ways of clustering,}\item{resClus:}{the complete output of the clustering function(s),}\item{comparClus:}{the adjusted Rand indices, used to compare the clusterings obtained for a same component.}}
##' @author Anne
##' @seealso \code{Mclust}, \code{adjustedRandIndex}, \code{kmeans}, \code{pam}, \code{pamk}, \code{hclust}, \code{agnes}, \code{cutree}
##' @export
##' @examples
##' data(icaSetCarbayo)
##' params <- buildMineICAParams(resPath="carbayo/", selCutoff=3)
##' 
##' ## compare kmeans clustering applied to A and data restricted to the contributing genes
##' ## on components 1 to 3
##' res <- clusterSamplesByComp_multiple(icaSet=icaSetCarbayo[,,1:3], params=params, funClus="kmeans",
##'                                      nbClus=2, clusterOn=c("A","S"), level="features")
##' head(res$clus)
##' 
##' 
clusterSamplesByComp_multiple <- function (icaSet, params, funClus=c("Mclust","kmeans","pam","pamk","hclust","agnes"), filename, clusterOn=c("A","S"), level=c("genes","features"), nbClus, metric="euclidean", method="ward", ...) {


    level <- match.arg(level)
    clus <- list()
    i <- 1
    cl1 <- comp <- cl2 <- indc <- cl <- cluscomp <- NULL
    
    resclus <- list()
    for (func in funClus) {
        for (on in clusterOn) {
            resC <- clusterSamplesByComp(icaSet=icaSet, params=params, funClus=func, clusterOn=on, level=level, nbClus=nbClus, metric=metric, method=method, ...)
            clus[[i]] <- resC$clus
            resclus[[i]] <- resC$resClus
            names(clus)[i] <- names(resclus)[i] <- paste(func,"_on",on,sep="")
            names(clus[[i]]) <- paste(names(clus[[i]]),"_",func,"_on",on,sep="")
             i <- i+1
           
        }
    }

    if (length(funClus)>1 | length(clusterOn)>1) {

        comparClus <- 
        foreach (comp=1:nbComp(icaSet), indc=indComp(icaSet), .combine=rbind) %do% {
            comparcluscomp <- 
                foreach(cl1=clus, ncl1=names(clus), .combine=rbind) %do% {
                    rr <- foreach (cl2=clus, ncl2=names(clus)) %do% {
                        #if (ncl1 != ncl2)
                            signif(mclust:::adjustedRandIndex(cl1[[comp]],cl2[[comp]]),3)
                    }
                    rr <- unlist(rr)
                    names(rr) <- paste(indc,"_",names(clus),sep="")
                    return(rr)
                }
            rownames(comparcluscomp) <- paste(indc,"_",names(clus),sep="")
            return(comparcluscomp)
        }
        colnames(comparClus) <- names(clus)

 

    }
    


    ## reorder clus
    clusbis <- 
    foreach(cl=clus) %do%  {
        clbis <- 
        foreach(cluscomp=cl, comp=1:length(cl)) %do% {
            meanByClus <- aggregate(A(icaSet)[,comp],by=list(cluscomp), FUN=mean)
            clusOrder <- structure(order(meanByClus[,2]), .Names = meanByClus[,1])
            cluscompbis <- clusOrder[as.character(cluscomp)]
            names(cluscompbis) <- names(cluscomp)
            return(cluscompbis)
        }
        names(clbis) <- names(cl)
        return(clbis)
    }
    names(clusbis) <- names(clus)
    clus <- clusbis

    
    clusall <- lapply(clus,data.frame, check.names=FALSE)
    clusall <- foreach(clus=clusall, .combine=cbind) %do% {clus}

    if (!missing(filename)) {
        filename <- paste(resPath(params),gsub(filename,pattern=".txt",replacement=""),".txt",sep="")
        write.table(clusall,file=filename, sep = "\t", quote=FALSE)
        if (length(funClus)>1 | length(clusterOn)>1) {
            write.table("",file=filename, sep = "\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=FALSE)
            write.table("Comparison of clusters using adjusted Rand indices:",file=filename, sep = "\t", quote=FALSE, append=TRUE, row.names=FALSE, col.names=FALSE)
           cat("",colnames(comparClus),"\n",file=filename, sep = "\t", append=TRUE)
            write.table(comparClus,file=filename, sep = "\t", quote=FALSE, append=TRUE, col.names=FALSE, row.names=TRUE)
        }
        
        
    }
    
    if (length(funClus)>1 | length(clusterOn)>1) 
        return(list(clus=clusall,resClus=resclus, comparClus=comparClus))
    else
        return(list(clus=clusall,resClus=resclus))
        

}
    
#### Plot Gaussian(s) fitted to data by Mclust
#### mc: result of Mclust
#### data: values on which Mclust was applied
#### nbBreaks: number of breaks for the histogram
##' Given a result of function \code{Mclust} applied to a numeric vector, this function draws the fitted Gaussian on the histogram of the data values.  
##'
##' A shapiro test p-value is added to the plot title. 
##' This function can only deal with at the most three Gaussian.
##' @title Plots an histogram and Gaussian fitted  by \code{\link[mclust]{Mclust}} 
##' @param mc The result of Mclust function applied to argument \code{data}
##' @param data A vector of numeric values
##' @param nbBreaks The number of breaks for the histogram
##' @param traceDensity If TRUE (default) density are displayed on the y-axis, else if FALSE counts are displayed on the y-acis
##' @param title A title for the plot
##' @param xlim  x-axis limits to be used in the plot
##' @param ylim  y-axis limits to be used in the plot
##' @param ... additional arguments for hist
##' @return NULL
##' @seealso \code{\link{hist}}, \code{\link[mclust]{Mclust}}
##' @author Anne Biton
##' @export
##' @examples
##' ## create a mix of two Gaussian
##' v <-c(rnorm(80,mean=-0.5,sd=1),rnorm(80,mean=1,sd=0.2))
##' ## apply Mclust
##' mc <- Mclust(v)
##' ## plot fitted Gaussian on histogram of v
##' plotMix(mc=mc,data=v,nbBreaks=30)
plotMix <- function (mc, data, nbBreaks, traceDensity = TRUE, title = "", xlim, ylim, ...) {

         nbMix = mc$G
	 #check normality : shapiro test
	 if (length(data)<=5000)
             shapiro.pval = shapiro.test(data)$p.value

	 p<-seq(min(data), to=max(data), length=1000)
	 d=cdens(modelName=mc$modelName, data=p, parameters=mc$parameters)

	 h <- hist(as.numeric(data), breaks = nbBreaks, freq = FALSE,
                   xlim= if (missing(xlim) || is.null(xlim)) xlim = c(min(unlist(data)),max(unlist(data))) else xlim,
                   ylim= if (missing(ylim) || is.null(ylim)) ylim= c(min(unlist(d)),max(unlist(d))) else ylim,
                   main = (paste(title,"\n", if (length(data)<=5000) paste("shapiro_pval = ", signif(shapiro.pval,digits = 3), sep = "") else "" , sep = "")),
                   cex.main = 0.9,
                   xlab = "contributions",
                   ylab = if (traceDensity) "Density" else "count",
                   yaxt = if (!traceDensity) "n" , ...)

  	 if (nbMix == 1) {
		 points(p, d[,1], type="l", col="blue", lty=1)
	 }
	 else if (nbMix == 2) {
		 temp1 <- mc$parameters$pro[1]
		 temp2 <- mc$parameters$pro[2]
		 points(p, d[,1]*temp1, type="l", col="green3", lty=1)
		 points(p, d[,2]*temp2, type="l", col="blue", lty=1)
	 }
	 else if (nbMix == 3) {
		 temp1 <- mc$parameters$pro[1]
		 temp2 <- mc$parameters$pro[2]
		 temp3 <- mc$parameters$pro[3]
		 points(p, d[,1]*temp1, type="l", col="green3", lty=1)
		 points(p, d[,2]*temp2, type="l", col="blue", lty=1)
		 points(p, d[,3]*temp3, type="l", col="magenta", lty=1)
	 }
	 else {
		apply(d, MARGIN = 2, function (x,p) {
						points(p, x, type="l", col="green3", lty=1)
						
					}
					, p = p
		)
	 }


         if (!traceDensity) {
             samples <- rownames(A)
             x = lm(h$density~h$counts)
             lab = seq(0,length(data)/2,by=5) #c(0,5,10,15,20,25,30,35, 40, 45, 50) #lab = c(0,5,10,15,20,25,30,40)
             lab_coord = sapply(lab,calcCoord,x = x)
             axis(side = 2, at = lab_coord, labels = lab)
         }
         return(NULL)

}

##' Given a result of function \code{Mclust} applied on several numeric vectors, this function plots the fitted Gaussian on their histograms.  
##'
##' This function can only deal with at the most three Gaussian
##' 
##' @title Plots the Gaussian fitted  by \code{Mclust} on several numeric vectors
##' @param mc  A list consisting of outputs of function \code{Mclust} applied to each column of \code{A}, if this argument is missing \code{Mclust} is applied by the function. 
##' @param A A data.frame of dimensions 'samples x components'.
##' @param nbMix The number of Gaussian to be fitted.
##' @param nbBreaks The number of breaks for the histogram.
##' @param xlim x-axis limits to be used in the plot.
##' @param pdf A pdf file.
##' @return A list of \code{Mclust} results.
##' @seealso \code{\link{plotMix}}, \code{\link{hist}}, \code{\link[mclust]{Mclust}}
##' @author Anne Biton
##' @export
##' @examples
##' A <-matrix(c(c(rnorm(80,mean=-0.5,sd=1),rnorm(80,mean=1,sd=0.2)),rnorm(160,mean=0.5,sd=1),
##'            c(rnorm(80,mean=-1,sd=0.3),rnorm(80,mean=0,sd=0.2))),ncol=3)
##' ## apply function Mclust to each column of A
##' mc <- apply(A,2,Mclust)
##' ## plot the corresponding Gaussians on the histogram of each column 
##' plotAllMix(mc=mc,A=A)
##' ## apply function Mclust to each column of A, and impose the fit of two Gaussian (G=2)
##' mc <- apply(A,2,Mclust,G=2)
##' ## plot the corresponding Gaussians on the histogram of each column 
##' plotAllMix(mc=mc,A=A)
##' ## When arg 'mc' is missing, Mclust is applied by the function
##' plotAllMix(A=A)

plotAllMix <- function (mc, A, nbMix = NULL, pdf,  nbBreaks = 20, xlim=NULL) {

    
	if (missing(mc)) {
                Alist <- as.list(data.frame(A))
                Alist <- lapply(Alist, function(x,n) {names(x) <- n; x}, n = rownames(A))
                mc <- lapply(Alist, Mclust, G = nbMix)
	}
        else
            if (length(mc) != ncol(A))
                stop("The length of 'mc' (which contains the results of Mclust for each component) must be the same as the number of columns in A.")

        
        if (!missing(pdf))
            pdf(pdf, width = 9, height = 9, title=file)
        
	par(mfrow = c(3,3))

        sapply(1:ncol(A),
               function(indComp, A, mc, xlim, nbBreaks)  {
                   plotMix(mc = mc[[indComp]], data = A[,indComp], nbBreaks = nbBreaks, xlim = xlim, title = paste("Component",indComp))
               },
               A = A,
               mc = mc,
               xlim = xlim,
               nbBreaks = nbBreaks
               )
        

	if (!missing(pdf))
            dev.off()	

	return(mc)



}



calcCoord = function(xx,x) return((x$coefficients[2]*xx+x$coefficients[1])[[1]])


##' Given a result of function Mclust applied on a numeric vector, this function add the fitted Gaussian to a previous plot.
##' This is an internal function called by \code{plotPosSamplesInComp}.
##'
##' 
##' This function can only deal with at the most three Gaussian.
##' @title  Plots the Gaussian fitted  by \code{\link[mclust]{Mclust}} 
##' @param mc The result of Mclust function applied to argument \code{data}
##' @param data The vector of numeric values on which was applied Mclust
##' @return NULL
##' @author Anne Biton
##' @examples
##' ## create a mix of two Gaussian
##' v <-c(rnorm(80,mean=-0.5,sd=1),rnorm(80,mean=1,sd=0.2))
##' ## apply Mclust
##' mc <- Mclust(v)
##' ## plot fitted Gaussian on histogram of v
##' hist(v, freq=FALSE)
##' MineICA:::plotMclust(mc=mc,data=v)
##' 
##' @keywords internal
plotMclust <- function (mc, data) {
	
	 nb.mix = mc$G
	 
	 p<-seq(min(data), to=max(data), length=1000)
	 d=cdens(modelName=mc$modelName, data=p, parameters=mc$parameters)
  	 if (nb.mix == 1) {
		 points(p, d[,1], type="l", col="blue", lty=1)
	 }
	 else if (nb.mix == 2) {
		 temp1 <- mc$parameters$pro[1]
		 temp2 <- mc$parameters$pro[2]
		 points(p, d[,1]*temp1, type="l", col="green3", lty=1)
		 points(p, d[,2]*temp2, type="l", col="blue", lty=1)
	 }
	 else if (nb.mix == 3) {
		 temp1 <- mc$parameters$pro[1]
		 temp2 <- mc$parameters$pro[2]
		 temp3 <- mc$parameters$pro[3]
		 points(p, d[,1]*temp1, type="l", col="green3", lty=1)
		 points(p, d[,2]*temp2, type="l", col="blue", lty=1)
		 points(p, d[,3]*temp3, type="l", col="magenta", lty=1)
	 }
         return(NULL)
}



####---------------------------------------------------------------------------------------------
#### Plot global histogram of each component and superimpose the histograms of 2 sub-groups for which the user want to compare the distributions.
#### The goal is to see if one of the 2 classes is overrepresented in one of the group formed by the represented cutoffs.
#### Input :
####	ids_a : ids of the samples in the second group
####	ids_b : ids of the samples in the first group
####	cutoffs 
####	A : matrix (samples*components)
####	file : NULL
####	keepLev : labels of the 2 groups (appear in the legend)
####	labelsComp : labels of the components or NULL (default) 
####---------------------------------------------------------------------------------------------
plotPos2classInComp <- function (ids_a, ids_b, cutoffs = NULL, A, file = NULL, keepLev, labelsComp = NULL, typePlot = "hist", mc.list = NULL, col = NULL) {

    if(!is.null(file))
        pdf(file, height = 8.267717, width = 29.7/2.54, paper = "a4r", title=file)
    par(mfrow = c(3,3))

    if (typePlot == "hist") {
	    A.sub_a = A[as.character(ids_a),]
	    A.sub_b = A[as.character(ids_b),]
	    samples = c(ids_a, ids_b)
	    for (i in 1:ncol(A)) {
		h<-hist(as.numeric(A[,i]),  freq = F,
			   xlim = c(min(unlist(A)),max(unlist(A))),
			   breaks = 20,
			   main = if (is.null(labelsComp)) paste("Component", i) else labelsComp[i],    
			   xlab = "contributions", yaxt = "n",ylab = "Frequency")

		ha<-hist(as.numeric(A.sub_a[,i]), freq = FALSE, breaks = h$breaks, plot = FALSE)
		hb<-hist(as.numeric(A.sub_b[,i]), freq = FALSE, breaks = h$breaks, plot = FALSE)
                
		ha$density<-ha$density*length(ids_a)/nrow(A)
		hb$density<-hb$density*length(ids_b)/nrow(A)
                
		plot(ha, freq=F, add = TRUE, col = if (is.null(col)) "#00FF0080" else col[1], yaxt = "n", lty = 0) 
		plot(hb, freq=F, add = TRUE, col = if (is.null(col)) "#00640070" else col[2], yaxt = "n", lty = 0)

		da = density(as.numeric(A.sub_a[,i]))
		db = density(as.numeric(A.sub_b[,i]))

                plotMclust(mc = mc.list[[i]], data = A[,i])

		x = lm(h$density~h$counts)
		lab = c(0,5,10,15,20,25,30,40) #lab = c(0,5,10,15,20,25,30,40)
		lab_coord = sapply(lab,calcCoord,x = x)
		axis(side = 2, at = lab_coord, labels = lab)

		if (!is.null(cutoffs)) for (j in cutoffs) abline(v = j, col = "darkslateblue", lwd = 2, lty = 2)
		legend("topright", legend = keepLev,  col = c(if (is.null(col)) c("#00FF0080","#00640080") else col), pch = c(15,15))
	    }
   }

   if(!is.null(file)) dev.off()

}



##' This internal function is called by \code{\link{plotDensOneAnnotInAllComp}} and \code{\link{qualVarAnalysis}} and is dedicated to the plot of the densities or boxplots using the package \code{\link{ggplot2}}. 
##'
##' 
##' @title Plots the densities or boxplots of the component contributions using \code{\link{ggplot2}}. 
##' @param annot a data.frame of dimensions 'samples x annotations' with one column corresponding to the component to trait ("comp" column) and one column corresponding to the groups of interest  ("interest" column)
##' @param colAnnot  the name of a column of the argument \code{annot} with the groups of interest
##' @param global a vector with the global distribution, e.g the contribution values of all samples on the component
##' @param keepLev the groups of interest, i.e the levels of the annotation \code{colAnnot} to be considered
##' @param comp.label the label of the component
##' @param colours a vector of colours indexed by the names of the groups of interest 
##' @param legend.title the title of the legend, if NULL (default) colAnnot is used
##' @param pval the p-value of the test, will be written in the title
##' @param test name of test that gave the p-value 
##' @param title.add a title to add to the automatically generated title
##' @param data_ref a data.frame similar to the argument \code{annot} but restricted to a set of reference samples
##' @param geneExpr a vector of values representative of the component, e.g the expression of the witness gene of the component 
##' @param geneRef the ID of the feature/gene \code{geneExpr} corresponds to, e.g the name of the witness gene
##' @param ylab A label for the y-axis (character)
##' @param trace_globalExpression if TRUE, geneExpr is plotted below the graph as a set of points whose colour is representative of the amount of expression, default is FALSE
##' @param trace_groupExpression if TRUE (default), geneExpr is plotted below the graph, by group, as a set of points whose colour is representative of the amount of expression
##' @param typePlot The type of plot, either "density" or "boxplot"
##' @param addPoints If TRUE, points are superimposed on the boxplots
##' @return A
##' @seealso \code{\link[ggplot2]{geom_density}}, \code{\link[ggplot2]{geom_boxplot}}, \code{\link[ggplot2]{geom_point}}
##' @author Anne Biton
##' @keywords internal
plotDens2classInComp_plotOnly <- function (annot,
                                           colAnnot,
                                           global,
                                           keepLev,
                                           comp.label = NULL,
                                           colours,
                                           legend.title=NULL,
                                           pval,
                                           test,
                                           title.add = NULL,
                                           data_ref = NULL,
                                           geneExpr = NULL,
                                           geneRef = NULL,
                                           ylab = NULL,
                                           trace_globalExpression = FALSE,
                                           trace_groupExpression = TRUE,
                                           typePlot = c("density","boxplot"),
                                           addPoints=FALSE
                                           ) {

    typePlot <- match.arg(typePlot)
    interest <- witGene <- y <- comp <- x <- NULL

   
                annot$interest <- annot[[colAnnot]]
                
                if (missing(keepLev))
                    keepLev <- levels(as.factor(annot$interest))

                if (!is.factor(annot$interest))
                   annot$interest <- factor(annot$interest, levels = keepLev)
                annot <- subset(annot, !(is.na(interest)))
  
                ### global distribution 
                if (typePlot == "density")
                    g <- ggplot(annot) +  geom_density(aes(x=global),data = as.data.frame(global), colour = "red",alpha = 0.4, position = "identity")

                else
                    g <- ggplot(annot)
 
		if (!is.null(geneExpr)) {

			### Trace global gene expression (on all samples)
			if (trace_globalExpression) {
                            df  = data.frame(geneRef = geneExpr, witGene = rep(-0.1,length(geneExpr)), global = as.numeric(global))

                             if (typePlot == "density")
                                g <- g + geom_point( aes(x = global, y = witGene,colour = geneRef),size = 3,data = df, shape = 15) + scale_colour_gradientn(name = geneRef,colours = maPalette(low = "blue",high = "red", mid = "yellow", k=50)) 
                            else
                                g <- g + geom_point( aes(y = global, x= witGene,colour = geneRef),size = 5,data = df, shape = 15) + scale_colour_gradientn(name = geneRef,colours = maPalette(low = "blue",high = "red", mid = "yellow", k=50)) + scale_x_discrete(geneRef) 

                        }
                        
                        
			### Trace sub-groups gene expression
                        ## one line by sub group, in each line draw the position of the samples (colored by expression) that are members of the group
			if (trace_groupExpression) {
                            
				for (i in 1:length(keepLev)) {
					gr = keepLev[i]
					s = rownames(subset(annot,interest %in% gr)) 
					yord = -0.2*i- if(trace_globalExpression) 0.1 else 0
					df_g  = data.frame(geneRef = geneExpr[s], y = rep(yord ,length(s)), global = global[s])

					g <- g  + geom_hline(yintercept = yord, linetype = 3, colour = "black") + geom_point( aes(x = global, y= y,colour = geneRef),size = 1.55,data = df_g, shape = 15)  + scale_colour_gradientn(name = geneRef,colours = maPalette(low = "blue",high = "red", mid = "yellow", k=50))#, breaks = signif(c(min(geneExpr),q[1],mean(geneExpr),q[2],max(geneExpr)),4), legend = FALSE) 

					g <- g + annotate("text", x = max(global)+0.06, y = yord, label = gr, size = 2) 
				}
				
			}
				


		}

                ### add densities for classes to be compared
		hist = FALSE
		if (hist) g <- g + geom_histogram(aes(x=comp, group = interest, linetype = interest, fill = interest), alpha = 0.4, position = "identity", binwidth = 0.05) + scale_y_continuous("Density") + scale_x_continuous(if (is.null(ylab)) "Sample contributions" else ylab)
		else  {
                    if (typePlot == "density")
                        g <- g + geom_density(aes(x=comp, group = interest, linetype = interest, fill = interest),
                                              alpha = 0.4,
                                              position = "identity",
                                              colour = "black") + scale_y_continuous("Density") + scale_x_continuous(if (is.null(ylab)) "Sample contributions" else ylab)
                    else if (typePlot == "boxplot") {
                        g <- g + geom_boxplot(aes(x=interest, y = comp, fill = interest),
                                              position = "identity", data = annot,
                                              colour = "black", outlier.shape = if (addPoints) NA else 16) + scale_x_discrete(colAnnot) + scale_y_continuous(if (is.null(ylab)) "Sample contributions" else ylab) + theme_bw()                       
                        if (addPoints)
                            g <- g + geom_jitter(aes(x=interest, y = comp, fill = interest), data = annot, color="#1A1A1A99", size=1.9, position=position_jitter(width=.2))+ theme_bw()#  + coord_flip()+position="jitter", 
                    }

                }

                ### add positions of reference/normal samples if provided 
		if (!is.null(data_ref)) {

                    if (typePlot == "density") {
			data_ref$y <- rep(x = -0.02,times = nrow(data_ref))
			g <- g + geom_point(aes(x=comp,y=y),data = data_ref,  fill = "green", shape = 24, size = 1.7)
                    }
                    else {                        
                        if (!is.null(geneExpr))
                            data_ref$x <- rep(x = 0.1,times = nrow(data_ref))
                        else
                            data_ref$x <- rep(x = 0.5,times = nrow(data_ref))

			g <- g + geom_point(aes(y=comp,x=x),data = data_ref,  fill = "green", shape = 23, size = 3)
                        
                    }
		}
                
		bop<-theme(
                           plot.title =element_text(size=13,face="bold"),
			legend.text=element_text(size=12, colour = "black"), 
                          legend.title = element_text(size=12,face="bold", colour = "black",hjust = 0), 
			axis.text.x = element_text(size = 10, colour = "black"), 
			legend.position="right",
			axis.text.y = element_text(size = 10, colour = "black"),
			axis.ticks.margin = unit(0.3, "cm"), 
			axis.title.x = element_text(size = 11), 
			axis.title.y = element_text(size = 11, angle = 90),
                        axis.text = element_text(size = 10)	
		)


                ## if number of characters in colAnnot is too long, cut for legend
                if (nchar(colAnnot)>20)
                    colAnnotlegend <- substr(colAnnot, start = 1, stop = 20)
                else
                    colAnnotlegend <- colAnnot

                g <- (g+bop
                      + scale_linetype_manual(values =  c(rep(c(1,6,8,5,3,2,4,7),floor(length(keepLev))), c(1,6,8,5,3,2,4,7)[length(keepLev)%%8]),#deal with the aspect of the density lines according to the number of groups to plot
                                              name = if (is.null(legend.title)) colAnnotlegend else legend.title,
                                              labels = paste(keepLev, " (",table(annot$interest)[levels(annot$interest)],")",sep=""),
                                              breaks = levels(annot$interest))
                      + scale_fill_manual(values = c(colours[keepLev]),
                                            name = if (is.null(legend.title)) colAnnotlegend else legend.title,
                                            labels = paste(keepLev, " (",table(annot[[colAnnot]])[levels(annot$interest)],")",sep=""),
                                          breaks = levels(annot$interest))
                      + ggtitle(paste(comp.label, "\n", test," : ", signif(pval,4), title.add, sep = "")))

                return(g)
	
}



##' Given a variable of the phenotype data (i.e vector of sample annotations), this function tests if the groups of samples formed by this variable are differently distributed on the components, in terms of contribution values. The distribution of the groups on the components are represented using density plots. It is possible to restrict the tests and the plots to a subset of samples and/or components.
##'
##' Wilcoxon or Kruskal-Wallis tests are applied depending on the number of groups of interest from the considered annotation (argument \code{keepLev}).
##' The plots are saved in individual files (one file per component) in arg 'path' if specified or in the current directory if not specificied. Ech individual file is nameb 'index-of-component_colAnnot.png.'
##' Recall that the sample-contribution values are contained in \code{A(icaSet)}, and the sample annotations in \code{pData(icaSet)}. 
##' 
##' One png image is created by plot and located in \code{path}. Each image is named by 'index-of-component_keepVar.png'.
##' @title Tests if groups of samples are differently distributed on the components and do the corresponding plots.
##' @param icaSet an object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}
##' @param keepVar a variable label, i.e the label of a column of the pheno data of icaSet available in  (\code{varLabels(icaSet)}) wich contains the groups of interest 
##' @param path the directory where the plots will be located 
##' @param samples a subset of sample names available in \code{samplenames(icaSet)}, if NULL (default) all samples are used
##' @param keepComp a subset of components available in \code{indComp(icaSet)}, if NULL (default) all components are used
##' @param keepLev the groups of interest, i.e the levels of the annotation \code{keepVar} to be considered
##' @param colours A vector of colours indexed by the elements of \code{keepLev}, if NULL the colours are generated automatically using annot2Color
##' @param legend.title title of the legend
##' @param cutoff The threshold p-value for statistical significance
##' @param doPlot if TRUE (default), the plots are drawn, else if FALSE only test results are returned
##' @param onlySign if TRUE (default), only the significant results are plotted
##' @param resTests a vector of p-values per component, if NULL (default) the p-values are calculated using Wilcoxon or Kruskal-Wallis test
##' @return Returns a data.frame of dimensions 'components x 1' containing the results of the non-parametric tests (Wilcoxon or Kruskal-Wallis tests) that test if the groups of interest are differently distributed on the components 
##' @seealso \code{\link{wilcoxOrKruskalOnA}}, code{\link{writeHtmlResTestsByAnnot}}, code{\link{wilcox.test}}, code{\link{kruskal.test}}
##' @keywords internal
##' @author Anne Biton
##' @examples \dontrun{
##' ## load an example of IcaSet 
##' data(icaSetCarbayo)
##' 
##' ## have a look at the sample annotations which are available
##' varLabels(icaSetCarbayo)
##' 
##' ## with doPlot=TRUE trace the contributions of the samples according
##' ## to their grade on the components
##' restests <- plotDensOneAnnotInAllComp(icaSet=icaSetCarbayo, keepVar="GRADE",
##'                                                 doPlot=FALSE)
##' }
##' 
plotDensOneAnnotInAllComp <- function (
                                       icaSet,
                                       keepVar,
                                       path=NULL,
                                       samples,
                                       keepComp,
                                       keepLev = NULL,
                                       colours = NULL,
                                       legend.title=NULL,
                                       doPlot = TRUE, 
                                       cutoff = 0.05,
                                       onlySign = TRUE,
                                       resTests
                                       ) {

    if (!missing(keepComp)) 
        icaSet <- icaSet[,,keepComp]
    

    if (!missing(samples))
        icaSet <- icaSet[,samples] 
        
    print(paste("Plot Densities of sample contributions on each IC according to the ", keepVar, " annotation"))

    
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)
    
    annot <- pData(icaSet)
    A <- A(icaSet)
    refSamples <- refSamples(icaSet)
    
    if (length(witGenes(icaSet))>0)
        witData <- datByGene(icaSet)[witGenes(icaSet),]
    else
        witData <- data.frame()
    labelsComp <- compNames(icaSet)
    
        if (length(refSamples)>0) {
            annot_ref <- subset(annot, rownames(annot) %in% refSamples)
            annot <- annot[which(!(rownames(annot) %in% refSamples)),]
        }
        else
            annot_ref <- NULL
    
    if (missing(resTests)) {
        ## if keepLev is not null, we subset the annotation data.frame to these classes
        if (!is.null(keepLev)) {
            annot.bis = subset(annot, annot[[keepVar]] %in% keepLev)
            A.bis = A[rownames(annot.bis),]
        } else {
            keepLev = as.character(unique(annot[[keepVar]])[!is.na(unique(annot[[keepVar]]))])
            A.bis = A[rownames(annot),]
            annot.bis = annot
        }
        
        ### Delete annotation concerning less than 2 samples
        levRepr <- summary(as.factor(annot[[keepVar]]))
        levNoRepr <- names(levRepr)[which(levRepr < 2)]
        annot <- subset(annot, !(keepVar %in% levNoRepr))
        keepLev <- keepLev[!(keepLev %in% levNoRepr)]
        keepLev <- keepLev[!is.na(keepLev)]
        nbLev <- length(keepLev)
        
        if (length(levNoRepr)>0)
            warning(paste("The annotation levels '",paste(levNoRepr,collapse=", "), "' map to less than 4 samples, they have not been taken into account.",sep =""))
        # if only one group : return NULL
        if (length(keepLev[!is.na(keepLev)]) == 1) {
            warning(paste("'ALl samples are in the same group after removing annotation levels corresponding to less than 4 samples in column or NA",keepVar, ".", sep =""))
            resTests <-  rep(NA,ncol(A.bis)); return(resTests)
        }
        

        annot.bis$interest = factor(annot.bis[[keepVar]], levels = keepLev)
        annot.bis <- subset(annot.bis,!is.na(as.character(annot.bis$interest)))
        A.bis <- A.bis[rownames(annot.bis),]

        resTests <- wilcoxOrKruskalOnA (A = A.bis, annot = annot.bis, colAnnot = "interest")

    }

    whichSign <- which(resTests <= cutoff)

    
    if (doPlot) {
          ##  vector supplying one color for each annotation level
          if (is.null(colours) | length(colours)==0)
              colours = annot2Color(annot.bis)
          
         for (i in whichSign) {

            print(paste("Comp",indComp(icaSet)[i]))
            annot.bis$comp = A.bis[rownames(annot.bis),i]
            if (length(refSamples) > 0)
                annot_ref$comp = A[rownames(annot_ref),][[i]]
            global <- A[,i]
            names(global) <- rownames(A)
 
            g <- plotDens2classInComp_plotOnly(
                                               annot = annot.bis,
                                               colAnnot = keepVar,
                                               global = global,
                                               comp.label = if (!is.null(labelsComp)) labelsComp[i] else paste("Component",indComp(icaSet)[i]),
                                               colours = colours,
                                               legend.title=legend.title,
                                               pval = resTests[i],
                                               test = if (nbLev>2) "Kruskal-Wallis test" else "Wilcoxon test" ,
                                               title.add = NULL,
                                               data_ref = annot_ref,
                                               geneExpr = unlist(witData[i,]),
                                               geneRef = witGenes(icaSet)[i])
            ggsave(plot=g,filename=paste(if(!is.null(path)) paste(path,sep=""),indComp(icaSet)[i],"_",keepVar,".png",sep=""),height = 5, width = 6)
          }
        }

    resTests = signif(resTests,4)
    return(resTests)
}

##' This function tests if the groups of samples formed by the variables (i.e sample annotations) are differently distributed on the components, in terms of contribution value (i.e of values in matrix \code{A(icaSet)}). The distribution of the groups on the components are represented using density plots. It is possible to restrict the tests and the plots to a subset of samples and/or components.
##' 
##' This function writes an HTML file containing the results of the tests and links to the corresponding density plots.  One png image is created by plot and located in the sub-directory plots of \code{path}. Each image is named by index-of-component_var.png.
##' Wilcoxon or Kruskal-Wallis tests are applied depending on the number of groups of interest from the considered annotation (argument \code{keepLev}).
##' @title Tests if groups of samples are differently distributed on the components according and do the corresponding plots.
##' @param icaSet an object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}
##' @param params An object of the class \code{\link[MineICA:MineICAParams-class]{MineICAParams}} containing the parameters of the analysis
##' @param path the directory where the plots will be located 
##' @param keepVar The variable labels to be considered, i.e a subset of (\code{varLabels(icaSet)}) 
##' @param samples a subset of sample names available in \code{samplenames(icaSet)}, if NULL (default) all samples are used
##' @param keepComp a subset of components available in \code{indComp(icaSet)}, if NULL (default) all components are used
##' @param legend.title_list A list of titles for each component, indexed by elements of argument \code{keepVar}, default is NULL
##' @param colours A vector of colours indexed by the variable levels, if missing the colours are automatically generated using \code{annot2Color}
##' @param doPlot if TRUE (default), the plots are drawn, else if FALSE only the tests are performed
##' @param pval.cutoff The threshold p-value for statistical significance
##' @param typeImage The type of image file where each plot is saved
##' @param filename A file where the results will be displayed in format HTML, if NULL no file is created 
##' @param onlySign if TRUE (default), only the significant results are plotted
##' @return Returns a data.frame of dimensions 'components x variables' containing the p-values of the non-parametric tests (Wilcoxon or Kruskal-Wallis tests) wich test if the samples groups defined by each variable are differently distributed on the components 
##' @seealso \code{\link{wilcoxOrKruskalOnA}}, \code{\link{writeHtmlResTestsByAnnot}}, \code{\link{plotDensOneAnnotInAllComp}}
##' @author Anne Biton
##' @examples \dontrun{
##' ## load an example of IcaSet
##' data(icaSetCarbayo)
##' ## have a look at the sample annotations which are available
##' varLabels(icaSetCarbayo)
##' ## create parameters, specifying the result path
##' params <- buildMineICAParams(resPath="carbayo/")
##' 
##' ## trace the contributions of the samples according to their cancer stages and gender on the components
##' ## make sure the arg 'path' exists in the directory contained in resPath(params)!
##' restests <- plotDensAllAnnotInAllComp(icaSet=icaSetCarbayo, keepVar=c("stage","SEX"),
##'                                       params=params, path="testPlotDens")
##' }
##' @keywords internal
plotDensAllAnnotInAllComp <- function (
                                       icaSet,
                                       params,
                                       path,
                                       keepVar = NULL,
                                       keepComp,
                                       samples,
                                       legend.title_list = NULL,
                                       colours = params["annot2col"],
                                       doPlot = TRUE,
                                       pval.cutoff = params["pvalCutoff"],
                                       typeImage = "png",
                                       filename = NULL,
                                       onlySign = TRUE
                                       ) {

    colAnnot <- NULL

    if (!missing(keepComp)) 
        icaSet <- icaSet[,,keepComp]
    

    if (!missing(samples))
        icaSet <- icaSet[,samples]

  
    if (!is.null(path)) {
        if (path != "" & substr(path,start=nchar(path),stop=nchar(path)) != "/") path = paste(path,"/",sep="")
        path <- gsub(path,pattern=resPath(params),replacement="")
        path <- paste(resPath(params),"/",path,sep="")
        system(paste("rm -r ", path,"plots/",sep=""), ignore.stderr=TRUE)
        pathplot <- paste(path,"plots/",sep="")
    }
    else stop("argument 'path' is missing for density plots")

  if (!is.null(keepVar))
      pData(icaSet) <- pData(icaSet)[,intersect(keepVar,varLabels(icaSet))]

  system(paste("mkdir",pathplot), ignore.stderr=TRUE)

  resTestsByAnnot <- 
      foreach (colAnnot=varLabels(icaSet), .combine=cbind, .errorhandling = "stop") %dopar% {

           resTests <- plotDensOneAnnotInAllComp(
                                              icaSet = icaSet,                                
                                              path=pathplot,
                                              keepVar = colAnnot,
                                              colours = colours,
                                              legend.title=if (!is.null(legend.title_list)) legend.title_list[[colAnnot]] else NULL,
                                              doPlot = doPlot,
                                              cutoff = pval.cutoff,
                                              onlySign = onlySign) 
      }

  
  colnames(resTestsByAnnot) <- varLabels(icaSet)
  rownames(resTestsByAnnot) <- indComp(icaSet)

  if (!is.null(filename)) {
      writeHtmlResTestsByAnnot(icaSet = icaSet,
                               params = params,
                               res = resTestsByAnnot,
                               path = path,
                               keepVar = keepVar,
                               pathplot = pathplot,
                               filename = filename,
                               cutoff = pval.cutoff,
                               typeImage = typeImage,
                               caption = "<b>The contributions of pre-established groups of tumors were compared using Wilcoxon rank-sum test for 2-classes comparison and Kruskal-Wallis for the comparison of more than 2 classes.</b> <br> Click on the pvalue to visualize the densities of the different sub-groups.",
                               onlySign = onlySign)


  }

  return(resTestsByAnnot)
}



##' Given a data.frame consisting of sample annotations, this function returns a vector which gives a colour per annotation level.
##'
##' Arbitrary colours are attributed to some specific annotations met by the author, and for the remaining annotation levels, the colours are attributed using packages \code{RColorBrewer} and \code{rcolorspace}.
##' @title Association of a colour with each annotation level
##' @param annot a data.frame containing the sample annotations (of dimension 'samples x annotations').
##' @return A vector of colours indexed by the annotation levels.
##' @author Anne Biton
##' @export
annot2Color <- function (annot) {
    allLev <- unique(unlist(as.data.frame(apply(annot,2,as.character), stringsAsFactors=FALSE)))
    allLev <- allLev[!(allLev %in% c("", " "))]

	
        mutation = c(non = "palevioletred1", oui = "palevioletred4", yes = "palevioletred4", no = "palevioletred1", y = "palevioletred4", n = "palevioletred1", "RAS" = "burlywood1",FGFR3="brown", FGFR3_RAS_no = "palevioletred1", RAS_yes="brown", FGFR3_yes = "palevioletred4")
        
	stade = c( Ta_lowgrade = "#CCEBC5", Ta_highgrade =  "#A8DDB5", Ta =  "#CCEBC5", T1 = "#7BCCC4" , "T2+" = "#2B8CBE",   "T2" = "dodgerblue2", "T3" = "#0868AC", "T4" =  "#084081",  TaT1 =  "#A8DDB5", "Ta/T1"= "#A8DDB5", 
	"T2+CIS+" = "slateblue4", "T2+CIS-" = "slateblue1", "T1CIS+" = "slategray1", "T1CIS-" = "slategray4", "TaCIS+" = "slategray3", "TaCIS-" = "slategray1",
	"T2+_FGFR3oui" = "aquamarine4", "T2+_FGFR3non" = "aquamarine3", "T1_FGFR3oui" ="darkslategray4", "T1_FGFR3non" ="darkslategray3" , "Ta_FGFR3oui" = "aquamarine3", "Ta_FGFR3non" = "aquamarine"
	,T0 = "green", "TaG1G2" = "#CCEBC5", "TalowG" = "#CCEBC5", cellLine = "#ffafaf", "TaG3T1" = "dodgerblue2", "T1a" = "deepskyblue", "T1b" = "dodgerblue3", "TahighG" = "darkslategray3", "T1highG" = "dodgerblue2",
        IA = "bisque", IB="bisque3", IIA = "darkseagreen1" , IIB ="darkseagreen3" , IIC = "darkseagreen3", IIIA = "darkslategray1", IIIB = "darkslategray3", IIIC ="darkslategray3" , IV =" darkslategrey",
        StageIA = "bisque", StageIB="bisque3", StageIIA = "darkseagreen1" , StageIIB ="darkseagreen2" , StageIIC = "darkseagreen3", StageIIIA = "darkslategray1", StageIIIB = "darkslategray2", StageIIIC ="darkslategray3" , StageIV =" darkslategrey"
)

        ovarian = c("TUMOR_FREE" = "deepskyblue", "WITH_TUMOR" =  "dodgerblue4", "Tooearly" = "gray93",  "Missing"="gray97",   "Sensitive" = "mediumspringgreen", "Resistant" = "mediumpurple2")
    years <- structure(c("#CAB2D6", "#FFBFFFFF", "#FFFF33", "#D95F02", "#A6761D"), .Names = c("2004", "2006", "2005", "2010", "2007"))
    
    yearMonth <- structure(c("orange", "#DC9A51", "violetred4", "orchid", "#E41A1C", 
                "magenta1", "forestgreen", "#B2DF8A", "#999999", "#6A3D9A", "springgreen4", 
                "turquoise1"), .Names = c("2004-11", "2004-04", "2006-01", "2006-04", 
                            "2005-12", "2006-08", "2010-05", "2010-04", "2010-06", "2007-03", 
                            "2010-09", "2010-03"))
    
	grade98 =  c(G0= "green",bas_grade = "#ffafaf", basgrade ="#ffafaf", hautgrade = "#cf0000", haut_grade = "#cf0000", 'NA' = "#ffffff", lowGrade ="#ffafaf", highGrade = "#cf0000",  LG ="#ffafaf", HG = "#cf0000") 
	grade73 =  c(G1 = "#ffafaf", G2 = "#ff6767", G3 = "#cf0000") 

	cis = c( absent = "#afafaf", 'associe' = "#272727", 'associe' = "#272727" )
	cis_orntoft =  c( '-' = "#afafaf", '+' = "#272727", ND1 = "#ffffff", ND2 = "#ffffff", 'NA' = "#ffffff") 
	sexe = c(femme = "#FF33FF", homme = "#3366FF", FEMALE = "#FF33FF", MALE = "#3366FF")

	necrose = c(non_necrose = "#CCFFFF", necrose = "#3399FF")
	prel = c(RTUV = "#CCFFFF", cystectomie = "#3399FF")
	taille = c("<=30" = "#CCFFFF", ">30" = "#3399FF","<=20" = "#CCFFFF", ">20" = "#3399FF", "CLI" = "#CCFFFF", "CCI" = "#3399FF")
	statutP38 = c(O="#3399FF",faible_phosphorylation="#3399FF", N= "#CCFFFF")
	recidive = c("FALSE" = "#FF6347", nonrecidive = "#FF6347",nonrecidive_noninv = "#FF6347", nonrecidive_inv = "#FF6347", recidive = "#8B0000", recidive_noninv = "#8B0000", recidive_inv = "#8B0000", "TRUE" = "#8B0000")
	
	batchs = c(no = "deepskyblue4", yes = "hotpink", biopsy = "grey20", ND = "gray80", Limit = "deepskyblue", OK = "darkslateblue", "DenatBeforeBioan" = "chartreuse",   "RIN_ExpertSoft"=  "darkorange4", denaturation = "gray20", notgood = "darkseagreen1", inversed = "darkseagreen", "goodQuality" = "gray20"
		,"27/09/2006" = "gray80", "22/02/2007" = "gray50", "28/08/2008" = "grey10")

	other = c("1" = "#ffafaf", "2" = "#cf0000", "3" = "red4", Control = "blue", "Primary_bladder_cancer" = "gray80",  "Recurrent_bladder_tumor" = "grey20", "Surrounding" = "darkslateblue" , M="darkslategray1", F="indianred1",
		  NO = "gray80",  YES = "grey20", NED = "snow2", DOD = "snow4",
		  "0" = "bisque" , "-1" = "paleturquoise1", "-2" = "paleturquoise2" , "1" = "palevioletred1", "2" = "palevioletred3", "3" = "palevioletred4",  "4" = "rosybrown4" ,"deletion_homozygote"= "cyan", "deletion_heterozygote" = "blue",TPM0 = "#ffafaf", TPM1 = "#ff6767", metastasis = "#cf0000") 
	patients_tumult =  c("ROZ"= "#8DD3C7", "HOL"= "#FFFFB3", "TAZ"= "#BEBADA", "TAM"=  "#FB8072" , "DGI"= "#80B1D3", "DRE"= "#FDB462", "LEP"= "#B3DE69", "DEM"= "mistyrose4", "HEN"= "#FCCDE5", "ROS"= "#D9D9D9", "VHO"=  "#CCEBC5", "VAN"= "#BC80BD", "Muscle" = "green",
m0.5="#A6CEE3", m1="#1F78B4", m1.5="#B2DF8A", m2="#33A02C", m2.5="#FB9A99", m3="#E31A1C", m3.5="#FDBF6F", m4 = "#FF7F00", m5 = "#CAB2D6", inf18 = "rosybrown", sup18="purple4",
        "<1995"= "#A6CEE3", "1995-2000" = "#33A02C",  "2000-2005"="#E31A1C", ">2005" = "#FF7F00")
        
   patientsDoublePrel = c("260GIL"= "#8DD3C7","6FON" = "mistyrose4","318AKL" = "#FB8072","39VAS" = "#80B1D3", "74COU" =  "#FDB462",
  HM = "#8DD3C7", IGR = "#FB8072", HF = "#FDB462" )  

	sousTypes = c(papillaire_pediculee = "paleturquoise4", papillaire_sessile = "turquoise4", papillaire = "turquoise4",solide_vegetante = "chartreuse4", ulceree_necrosee = "chocolate", ulceree_vegetante = "89", ulceree = "89",  plane = "gray30", autre = "gray30",autres = "gray30", "ER+"="bisque","ER-"="palevioletred1",
	 classique = "paleturquoise4", epidermoede = "turquoise4", epidermoide = "turquoise4",micropapillaire = "chartreuse4", neuroendocrine = "violetred1", sarcomatoide = "violetred4",
        "endophytique" = "violetred1", "largesMassifs" = "chartreuse4", "plages" = "yellow", "travees_massifsPeuCohesifs" = "chocolate1", 
"travees_petitsMassifs" = "chocolate4")
        fgfr3mut <- c("WT" = "darkolivegreen1", "S249C" = "hotpink1", "S249C+Y375C" = "darkorchid1", "R248C"="cadetblue1", "Y375C" = "darkorchid4", "S249Cminor"= "hotpink1",
"K652E", "S249C+R248C" = "cadetblue4", "Y375Cminor" = "darkorchid", "WTminortobeconfirmed"= "darkolivegreen1"   )

        rbmut <- c(consMut = "red",   delLOH = "orchid",   methyl = "seagreen2",   oneDel = "royalblue1",   oneMut = "palegreen2", PTCslice = "peachpuff",   twoDel = "royalblue4",   twoMut = "palegreen4")
  
    stageAML <- c("M0Undifferentiated" = "orange", "M1" = "maroon1", "M2" = "magenta4", "M3" = "black", "M4" = "chartreuse2", "M5" = "cyan2", "M6" = "blue3", "M7" = "brown", NotClassified="gray60")

        surv <- c(LIVING = "#ffafaf", DECEASED = "#cf0000")
        cyto <- c("Normal"= "aquamarine", "t(8;21)" = "#377EB8", "del(5q)/5q-"= "#4DAF4A", "inv(16)"= "#984EA3", "del(7q)/7q-" = "#FF7F00", "Trisomy8" = "#FFFF33", "t(15;17)" = "#A65628", "complex" = "#8DD3C7" , "t(9;11)" = "#8DD3C7")##normal E41A1C
        tutu <- c("OtherCDNotTested" = "#E7E1EF", "OtherCDNegative" = "#C994C7" , "OtherCDPositive" = "#DD1C77")
        plate <- c("734"="blue4","735"="red","736"="darkred","748"="seagreen4","751"="darkmagenta","760"="chartreuse4","744"="darkorange","740"="deeppink1", "0734"="blue4","0735"="red","0736"="darkred","0748"="seagreen4","0751"="darkmagenta","0760"="chartreuse4","0744"="darkorange","0740"="deeppink1")
        amlothers <- c("Intermediate/Normal"= "aquamarine", "Favorable" = "aquamarine4", "Poor" = "chocolate1",  "NPMcNegative" = "#C994C7", "NPMcPositive" = "#DD1C77","PML-RAR"="#1B9E77", "BCR-ABL"="#D95F02", "CBF-Beta"="#7570B3", "8"="#E7298A",
"AML1-ETO"="#66A61E", "MLL"="#E6AB02", "TEL-AML1"="#A6761D", "Other"="#666666", 
"del(20q)"="#7FC97F", "-5ordel(5q)"="#BEAED4", "-7ordel(7q)"="#FDC086", 
"4"="#BF5B17", "5"="#666666") #"1"="#FFFF99", "2"="#386CB0", "3"="#F0027F", 
        aml <- c(stageAML, surv, cyto, tutu, plate, amlothers)

	merge.colours = c(mutation, stade, grade98, grade73, cis, cis_orntoft, sexe, necrose, taille, prel, statutP38, recidive, batchs, other, patients_tumult, sousTypes, patientsDoublePrel, fgfr3mut,rbmut, ovarian, aml, years, yearMonth)

    up2or <- allLev
    names(up2or) <- toupper(allLev)

    names(merge.colours) <- toupper(names(merge.colours))
    
    leftLev <- allLev[!(toupper(allLev) %in% toupper(names(merge.colours)))]
    qual <- c(brewer.pal(8,"Accent"),
              brewer.pal(8,"Dark2"),
              brewer.pal(12,"Paired"),    
               brewer.pal(9,"Set1"),     
              brewer.pal(8,"Set2"),       
              brewer.pal(12,"Set3"),
              (heat_hcl(5, c. = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))[1:4],
              (terrain_hcl(5, c. = c(65, 0), l = c(45, 95), power = c(1/3, 1.5)))[1:4],   
              (heat_hcl(4, h = c(0, -100), l = c(75, 40), c. = c(40, 80), power = 1)),
              (cm.colors(4)),
              rainbow_hcl(4, start = 60, end = 240),
              rainbow_hcl(4, start = 90, end = -30))
    leftLevCol <- sample(qual, size = length(leftLev), replace = length(leftLev)>length(qual))
    names(leftLevCol) <- leftLev

    avLevCol <- merge.colours[intersect(toupper(allLev),names(merge.colours))]
    names(avLevCol) <- up2or[intersect(toupper(allLev),names(merge.colours))]

    levCol <- c(avLevCol,leftLevCol)
    return(levCol)
    

        
}

cutoff.biMod = function (mixOnComp, A, ind.comp, proba.cutoff, use.classification) {

	mixOnComp2 = mixOnComp

	if (use.classification) {
		classification = mixOnComp2$classification
		names(classification) = row.names(A)
		mix1.sel = names(classification)[which(classification == 1)]
		mix2.sel = names(classification)[which(classification == 2)]

	}
	else {
            
            z.comp2 = mixOnComp2$z
            proba.mix1.comp2 = z.comp2[,1]
            proba.mix2.comp2 = z.comp2[,2]
            names(proba.mix1.comp2) = row.names(A)
            names(proba.mix2.comp2) = row.names(A)

            proba.mix1.comp.sel = proba.mix1.comp2[proba.mix1.comp2 >= proba.cutoff]
            proba.mix2.comp.sel = proba.mix2.comp2[proba.mix2.comp2 >= proba.cutoff]
            mix1.sel = names(proba.mix1.comp.sel)
	    mix2.sel = names(proba.mix2.comp.sel)

	}

	comp = A[[ind.comp]]
	names(comp) = row.names(A)

	mix1 = sort(comp[mix1.sel])
	mix2 = sort(comp[mix2.sel])
	cutoff.left = max(mix1)
	cutoff.right = min(mix2)
	return(c(left=cutoff.left, right=cutoff.right))


	    
}


##### ===============================
#####  With GGPLOT2
##### ===============================

##' Given a sample annotation (e.g a tumor specific stage), this function plots the positions of the corresponding samples (e.g the subset of samples having this tumor stage) within the histogram of the global sample contributions. This function is called by \code{\link{plotPosOneAnnotInComp_ggplot}} and is only dedicated to the plot of the histogram using the package \code{\link{ggplot2}}. 
##'
##' 
##' @title Plots the position of a subset of samples in the histogram of all samples using \code{\link{ggplot2}}. 
##' @param annot a data.frame of dimensions 'samples x annotations' with one column corresponding to the component to trait ("comp" column) and one column corresponding to the groups of interest  ("interest" column)
##' @param colAnnot  the name of a column of the argument \code{annot} with the groups of interest
##' @param selLev the name of the group of interest
##' @param comp a vector of sample contributions
##' @param colSel the colour of the histogram of the group of interest, default is "red"
##' @param colAll the colour of the global histogram 
##' @param geneExpr a vector of values representative of the component, e.g the expression of the witness gene of the component 
##' @param geneRef the ID of the feature/gene \code{geneExpr} corresponds to, e.g the name of the witness gene
##' @param title A title for the plot
##' @param binwidth set the width of the bins, see \code{\link[ggplot2]{geom_histogram}}
##' @param ... other parameters given to  \code{\link[ggplot2]{geom_histogram}}
##' @return An object of class ggplot2 containing the histogram
##' @seealso \code{\link[ggplot2]{geom_histogram}}
##' @keywords internal
##' @author Anne Biton
plotPosOneAnnotLevInComp_ggplot <- function(annot, colAnnot, selLev, comp, title = NULL, colSel = "red", colAll = "grey74", binwidth = 0.1, geneExpr = NULL, geneRef = NULL, ...) { 

    fake <- interest <- global <- y <- gene <- NULL

    annot$comp <- comp[rownames(annot)]
    annot$interest <- annot[[colAnnot]]
    annot$fake <- factor(c(rep("1",nrow(annot)),rep("2",0)), levels = c("1","2"))
    annotbis <- subset(annot, annot[[colAnnot]] == selLev & !is.na(annot[[colAnnot]]))
    annotbis$interest <- factor(annotbis$interest,levels=selLev)
    
    g <- ggplot(annot) +
          # Only to have adequat legend
          geom_histogram(aes(x=comp, y = (..count..), fill  = fake), colour = "black", binwidth = binwidth, ...)   + #, colour = "black" position = "identity", 
             scale_fill_manual (values = c("1"=colAll[[1]],"2"=colSel[[1]]), labels = c("All", selLev), breaks = c("1","2"), name = colAnnot) +
          # true needed histograms
          geom_histogram(aes(x=comp, y = (..count..)), colour = "gray20",fill = colAll, binwidth = binwidth, ...) + #fill = colAll, position = "identity", 
          geom_histogram(aes(x=comp, y = (..count..), groups = interest), colour = "black", fill = colSel, data = annotbis, position = "identity", binwidth = binwidth, ...) 

    bop<-theme(
               plot.title =element_text(size=13,face="bold"),
            legend.text=element_text(size=12, colour = "black"), 
              legend.title = element_text(size=12,face="bold", colour = "black",hjust = 0), 
            axis.text.x = element_text(size = 10, colour = "black"), 
            legend.position="right",
            #legend.background=theme_rect(colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.ticks.margin = unit(0.3, "cm"), 
            axis.title.x = element_text(size = 11), 
            axis.title.y = element_text(size = 11, angle = 90),
            axis.text = element_text(size = 10),
            title = paste(selLev, if (!is.null(title)) paste("in comp", "\n", title) else "") 
    )
   
    g <- g + bop + xlab("Sample contributions")  #geom_point(aes(x = comp,  y = y), shape = 22,  data = datafake) +

    if (!is.null(geneExpr)) {

            df  = data.frame("gene" = geneExpr, y = rep(-2,length(geneExpr)), global = comp)

            g <- g + geom_point( aes(x = global, y= y,colour = gene),size = 1.55,data = df, shape = 15) +
                     scale_colour_gradientn(name = geneRef,colours = maPalette(low = "blue",high = "red", mid = "yellow", k=50)) +
                     geom_hline(yintercept = -2, linetype = 3, colour = "black") 

            df_g  = data.frame(gene = geneExpr[rownames(annotbis)], y = rep(-4,nrow(annotbis)), global = annotbis$comp)

            g <- g  + geom_hline(yintercept = -4, linetype = 3, colour = "black") +
            geom_point( aes(x = global, y= y,colour = gene),size = 1.55,data = df_g, shape = 15)  +
            scale_colour_gradientn(name = geneRef,colours = maPalette(low = "blue",high = "red", mid = "yellow", k=50))#, breaks = signif(c(min(geneExpr),q[1],mean(geneExpr),q[2],max(geneExpr)),4), legend = FALSE) 

            g <- g + annotate("text", x = max(df$global)+0.1, y = -6, label = selLev, size = 3) 
            g <- g + annotate("text", x = max(df$global)+0.1, y = -3, label = "All", size = 3) 

        }

    
    return(g)

}

##' Given a variable of the phenoData, this function tests if the groups of samples formed by this variable are differently distributed, in terms of contribution value (i.e of values in matrix \code{A(icaSet)}), on the components. The distribution of the groups on the components are represented using density plots. It is possible to restrict the tests and the plots to a subset of samples and/or components.
##'
##' Wilcoxon or Kruskal-Wallis tests are applied depending on the number of groups of interest from the considered annotation (argument \code{keepLev}).
##' One png image is created by plot and located in \code{path}. Each image is named by component-of-component_colAnnot.png.
##' @title Tests if groups of samples are differently distributed on the components and do the corresponding plots.
##' @param icaSet An object of class \code{IcaSet}
##' @param params An object of the class \code{MineICAParams} containing the parameters of the analysis
##' @param colAnnot a variable label, i.e one of the variables available in  (\code{varLabels(icaSet)}) containing the groups of interest 
##' @param samples a subset of sample names available in \code{samplenames(icaSet)}, if NULL (default) all samples are used
##' @param keepComp a subset of components available in \code{indComp(icaSet)}, if NULL (default) all components are used
##' @param keepLev the groups of interest, i.e the levels of the variable \code{colAnnot} to be considered
##' @param colAll The colour of the global histogram, default is "grey74"
##' @param file the file where the histograms will be plotted
##' @param addExpr if TRUE (default) the expression profiles of the witness genes of each component are added below the plot
##' @param binwidth binwidth of the histogram (default is 0.1)
##' @param ... other parameters for geom_histogram function from ggplot2 package
##' @return NULL
##' @keywords internal
##' @seealso \code{\link{plotPosOneAnnotLevInComp_ggplot}}, \code{\link{geom_histogram}}
##' @author Anne Biton
plotPosOneAnnotInComp_ggplot<- function(icaSet, params, colAnnot, keepLev= NULL, keepComp, samples, colAll = "grey74", binwidth = 0.1, addExpr = TRUE, file = NULL,  ...) {
    
    if (!missing(keepComp))
        icaSet <- icaSet[,,keepComp]
    if (!missing(samples))
        icaSet <- icaSet[,samples]
    
    annot <- pData(icaSet)
    if (!is.null(keepLev)) {
        keepLev <- keepLev[keepLev %in% annot[[colAnnot]]]
        lev <- keepLev
    }
    else lev <- unique(annot[[colAnnot]])

    Alist <- Alist(icaSet)

    labelsComp <- compNames(icaSet)

    if (addExpr) {
        wit <- witGenes(icaSet)
        if (is.null(wit) | length(wit)==0) stop("Attribute witGenes is not available in IcaSet object and therefore the expression can't be plotted, please fill the corresponding attribute or put addExpr=FALSE to continue.")  
        witD <- datByGene(icaSet)[wit,rownames(annot)]
    }
    

    if(!is.null(file))
        pdf(file, height = 8.267717, width = 29.7/2.54, paper = "a4r", title=file)
    
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)


    indC <- labelComp <- NULL
    i <- 1
    foreach (comp = Alist, labelComp = labelsComp, indC = 1:length(Alist), .errorhandling='stop') %dopar% {
        
        comp <- comp[rownames(annot)]
        if (addExpr)
            wit <- witD[indC,]
        else
            wit <- NULL

        for (lev in keepLev) {    
            g <- plotPosOneAnnotLevInComp_ggplot(annot = annot, colAnnot = colAnnot, keepLev = keepLev, comp = comp, title = labelComp, colSel = params["annot2col"][keepLev], colAll = colAll, binwidth = 0.1, geneExpr = unlist(wit), geneRef = rownames(witD)[indC], ...)
            indPlot <- i
            ind = indPlot%%4
            if (ind == 1){
                grid.newpage();
                pushViewport(viewport(layout = grid.layout(2, 2)));
                vp = vplayout(1,1)
            }
            else if (ind == 2) vp = vplayout(1,2)
            else if (ind == 3) vp = vplayout(2,1)
            else if (ind == 0) vp = vplayout(2,2)

            print(g, vp = vp)
            i <- i+1

        
        }

    }
    
    if(!is.null(file))
        dev.off()


}

##### =================================
##### Input:
#####   icaSet: object of class IcaSet providing the results of an ICA decomposition
#####   params: object of class MineICAParams which provides the parameters of the ICA analysis
#####   samples: (character) ids of a subset of samples from icaSet whose contributions will be plotted on the ICs.
#####   keepComp: (numeric) index of the components that will be considered
#####   colAll: (character) color of the histogram of all samples
#####   addExpr: TRUE if expression levels of the witness genes of each component have to be plotted below the histograms
#####   file: (character) name of the file 
#####   ...: other parameters for geom_histogram function from ggplot2 package
##### =================================

plotPosSamplesInComp_ggplot <- function(icaSet, params, samples, keepComp,  colAll = "grey74", colSel = "red", labSel = NULL, binwidth = 0.1, addExpr = TRUE, file = NULL,  ...) {
    

    indC <- fake <- labelComp <- global <- gene <- NULL
    
    if (!missing(keepComp))
        icaSet <- icaSet[,,keepComp]

    #A <- getA(icaSet)
    Alist <- Alist(icaSet) #as.list(A); Alist <- lapply(Alist, function(x,n) {names(x) <- n; x}, n = sampleNames(icaSet))

    inter <- intersect(samples, sampleNames(icaSet))
    if (length(inter) == 0)
        stop ("The sample ids are not available in the IcaSet object")
    else
        if (length(inter) <= length(samples))
            warning(paste("The samples ", paste(setdiff(samples, sampleNames(icaSet)),collapse=", ",sep=""),"are not available in the IcaSet object"))
    samples <- inter
    
    labelsComp <- compNames(icaSet)

    if (addExpr) {
        witD <- witGenes(icaSet)
        if (is.null(witD) | length(witD)==0) stop("Attribute witGenes is not available in icaSet object and therefore the expression cannot be plotted, please fill the corresponding attribute or put addExpr=FALSE to continue.")  
        witD <- datByGene(icaSet)[witD,sampleNames(icaSet)]
    }
    
   

    
    if(!is.null(file))
        pdf(file, height = 8.267717, width = 29.7/2.54, paper = "a4r", title=file)
    
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)


    i <- 1
    foreach (comp = Alist, labelComp = labelsComp, indC = 1:length(Alist), .errorhandling='stop') %dopar% {
        
        comp <- comp[sampleNames(icaSet)]
        if (addExpr)
            geneExpr <- unlist(witD[indC,])
        else
            geneExpr <- NULL
        geneRef <- rownames(witD)[indC]
        
        annot <- data.frame(comp = comp, fake = as.factor(rep("a",nrow(A))), row.names = sampleNames(icaSet))
        annotbis <- annot[samples,]
    
        g <- ggplot(annot) +
             geom_histogram(aes(x=comp, y = (..count..)), colour = "gray20",fill = colAll, binwidth = binwidth, ...) + #fill = colAll, position = "identity", 
             geom_histogram(aes(x=comp, y = (..count..), groups = fake), colour = "black", fill = colSel, data = annotbis, position = "identity", binwidth = binwidth, ...) 

        bop<-theme(
                  plot.title =element_text(size=13,face="bold"),
                  legend.text=element_text(size=12, colour = "black"), 
                  legend.title = element_text(size=12,face="bold", colour = "black",hjust = 0), 
                  axis.text.x = element_text(size = 10, colour = "black"), 
                  legend.position="right",
                                        #legend.background=theme_rect(colour = "black"),
                  axis.text.y = element_text(size = 10, colour = "black"),
                  axis.ticks.margin = unit(0.3, "cm"), 
                  axis.title.x = element_text(size = 11), 
                  axis.title.y = element_text(size = 11, angle = 90),
                  axis.text = element_text(size = 10)                
                  )
      
   
        g <- g + bop + xlab("Sample contributions")  + ggtitle(if (!is.null(labelComp)) paste("Comp", "\n", labelComp) else paste("Comp",indComp(icaSet)[indC]))  

        if (!is.null(geneExpr)) {
            ### Plot global gene expression (on all samples)
            df  = data.frame("gene" = geneExpr, y = rep(-2,length(geneExpr)), global = comp)
            g <- g + geom_point( aes(x = global, y= y,colour = gene),size = 1.55,data = df, shape = 15) +
                     scale_colour_gradientn(name = geneRef,colours = maPalette(low = "blue",high = "red", mid = "yellow", k=50)) +
                     geom_hline(yintercept = -2, linetype = 3, colour = "black") 

            ### Trace expression of the samples belonging to the samples subset
            df_g  = data.frame(gene = geneExpr[rownames(annotbis)], y = rep(-4,nrow(annotbis)), global = annotbis$comp)
            g <- g  + geom_hline(yintercept = -4, linetype = 3, colour = "black") +
                      geom_point( aes(x = global, y= y,colour = gene),size = 1.55,data = df_g, shape = 15)  +
                      scale_colour_gradientn(name = geneRef,colours = maPalette(low = "blue",high = "red", mid = "yellow", k=50))#, breaks = signif(c(min(geneExpr),q[1],mean(geneExpr),q[2],max(geneExpr)),4), legend = FALSE) 

            ## add abline and 
            g <- g + annotate("text", x = max(df$global)+0.1, y = -6, label = if (is.null(labSel)) "subset" else labSel, size = 3) 
            g <- g + annotate("text", x = max(df$global)+0.1, y = -3, label = "All", size = 3) 

        }
        ### add positions of reference/normal samples if provided 
        if (refSamples(icaSet) != "") {
            
            data_ref <- annot[refSamples(icaSet),]
            y <- rep(x = -0.5,times = length(refSamples(icaSet)))
            data_ref$y = y
            g <- g + geom_point(aes(x=comp,y=y),data = data_ref,  fill = "green", shape = 24, size = 1.7) 
        }


        #### -----------------------

        
        indPlot <- i
        #choose viewport for a 2*2 frame
        ind = indPlot%%4
        if (ind == 1){
            grid.newpage();
            pushViewport(viewport(layout = grid.layout(2, 2)));
            vp = vplayout(1,1)
        }
        else if (ind == 2) vp = vplayout(1,2)
        else if (ind == 3) vp = vplayout(2,1)
        else if (ind == 0) vp = vplayout(2,2)

        print(g, vp = vp)
        i <- i+1

        
        

    }
    
    if(!is.null(file))
        dev.off()


}



##### ===============================
#####  With base plots
##### ===============================


##' This function plots the positions of several groups of samples across all the components of an \code{\link[MineICA:IcaSet-class]{icaSet}} object.
##' 
##' For each subgroup of samples this function plots their positions within the histogram of the global sample contributions.
##'
##' The values of interest are the sample contributions across the components, i.e across the columns \code{A(icaSet)}.
##' 
##' If argument \code{resClus} is not missing, the association between the clusters and the sub-groups of samples is tested using a chi-square test. The p-values of these tests are available in the title of each plot.
##' 
##' @title Histograms of sample subsets  
##' @param samplesByGroup A list whose elements are vector of sample names, these sample names must be available in \code{sampleNames(icaSet)}. The list should be indexed by the name of the corresponding groups. 
##' @param labGroups A vector of group names, will be used to add names to \code{sampleByGroup} if \code{names(samplesByGroup)} is NULL. 
##' @param icaSet An object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}
##' @param keepComp A subset of components available in \code{indComp(icaSet)}, if NULL (default) all components are used
##' @param file A pdf file
##' @param breaks The number of breaks to be used in the histograms
##' @param colSel The colour of the histogram of the group of interest, default is "red"
##' @param colAll The colour of the global histogram, default is "grey74"
##' @param resClus A list containing the outputs of function \code{\link{clusterSamplesByComp}}, which consists of results of clustering applied to matrix A of argument \code{icaSet}.
##' @param funClus Specifies the clustering method used, either \code{"Mclust"} or \code{"kmeans"}. If \code{resClus} is not missing, equals \code{resClus$funClus}.
##' @param titlesup Additional title for the histograms
##' @param ... Additional parameters for function \code{\link{hist}}
##' @return NULL
##' @seealso \code{\link{hist}}, \code{\link[MineICA:IcaSet-class]{IcaSet}}
##' @author Anne Biton
##' @export
##' @keywords internal
##' @examples \dontrun{
##' ## load an example of IcaSet
##' data(icaSetCarbayo)
##' 
##' ## selection of sample groups according to annotations STAGE 
##' samplesByGroup <- lapply(split(pData(icaSetCarbayo),pData(icaSetCarbayo)[c("STAGE")]), rownames)
##' # select groups including at least 2 samples
##' samplesByGroup <- samplesByGroup[which(unlist(lapply(samplesByGroup,length))>1)]
##' 
##' ## clustering of samples according to A using Mclust imposing two Gaussian
##' resClus <- clusterSamplesByComp(icaSet=icaSetCarbayo,funClus="Mclust", nbClus=2, clusterOn="A")
##'
##' ## Plot positions of the groups in 5th component
##' pdf(file="stageOnIC5.pdf", height = 8.267717, width = 29.7/2.54, paper = 'a4r', title="stageOnIC5")
##' plotPosSamplesInComp(samplesByGroup=samplesByGroup, icaSet=icaSetCarbayo, funClus="Mclust",
##'                      resClus = resClus, keepComp=5)
##' dev.off()
##' }
plotPosSamplesInComp <- function (samplesByGroup,
                                  labGroups = NULL,
                                  icaSet,
                                  keepComp=indComp(icaSet),
                                  file = NULL,
                                  breaks = 20,
                                  colAll = "grey74",
                                  colSel = "red",
                                  titlesup = "",
                                  resClus,
                                  funClus = c("Mclust","kmeans"),
                                  ...) {

   if (!missing(keepComp)) {
       icaSet <- icaSet[,,keepComp]
       if (!missing(resClus)) {
           resClus$clus <- resClus$clus[as.character(keepComp)]
           resClus$resClus <- resClus$resClus[as.character(keepComp)]
       }
   }
   else
       keepComp <- indComp(icaSet)

   if (!missing(resClus)) 
       funClus <- resClus$funClus

   
    if (length(witGenes(icaSet)) == nbComp(icaSet))
        addExpr <- TRUE
    else
        addExpr <- FALSE
   

   if (!is.list(samplesByGroup))
        stop("The samples have to be provided using a list. This list has to be indexed by the labels of the corresponding groups, if not, these labels have to be provided in 'labGroups'.")

    nullLabGroups <- is.null(labGroups)
    if (is.null(names(samplesByGroup)) & nullLabGroups)
        stop ("Labels of the groups are provided neither in names of 'samplesByGroup', nor in 'labGroups'.")
    if (nullLabGroups)
       labGroups <- names(samplesByGroup)
    else {
        if (length(labGroups) != length(samplesByGroup))
            stop("Lengths of 'labGroups' and 'samplesByGroup' have to be the same.")
        names(samplesByGroup) <- labGroups
    }

    if ((length(colSel) == 1) |
        (length(colSel)>1 & length(colSel)<length(samplesByGroup)) |
        (length(colSel) == length(samplesByGroup) & !all(labGroups %in% names(colSel)))) {
        warning("The colours provided in 'colSel' have to be indexed by the names of the groups. If not, only the first color is used")
        colSel <- rep(colSel[1],length(labGroups))
        names(colSel) <- labGroups
    }
    

    samplesByGroup <- 
    lapply(samplesByGroup,
           function(samples, icaSet) {
               inter <- intersect(samples, sampleNames(icaSet))
               if (length(inter) == 0)
                   stop ("The sample ids are not available in the IcaSet object")
               else
                   if (length(inter) < length(samples))
                       warning(paste("The samples ", paste(setdiff(samples, sampleNames(icaSet)),collapse=", ",sep=""),"are not available in the IcaSet object"))                   

               return(inter)
           }
           , icaSet = icaSet)

   if (!is.null(refSamples(icaSet)) & length(refSamples(icaSet))>0) {
       indallref <- lapply(samplesByGroup,
                           function(x,sref) {
                               if (all(x %in% sref))
                                   TRUE
                               else
                                   FALSE
                           },
                           sref=refSamples(icaSet))
       
       indallref <- which(unlist(indallref))
       if (length(indallref)>0) {
           samplesByGroup <- samplesByGroup[-indallref]
           labGroups <- names(samplesByGroup)
       }
       
       
   }

    
    
    labelsComp <- compNames(icaSet)

   if(!is.null(file))
        pdf(file, height = 8.267717, width = 29.7/2.54, paper = "a4r", title=file)
    
    witD <- NULL
    if (addExpr) {
        witD <- datByGene(icaSet)[witGenes(icaSet),]
        par(mfrow = c(2,2))
    }
    else
        par(mfrow = c(3,3))
        
   Al <- Alist(icaSet)
       
    ### "par" parameters according to frame
    parMfg <- c("1"="par(mfg=c(1,1))",
                "2"="par(mfg=c(1,2))",
                "3"="par(mfg=c(2,1))",
                "0"="par(mfg=c(2,2))")
    # for maColorBar
    parFig <- c("1"="par(fig = c(0.35,0.45,0.65,1), new = TRUE)",
                "2"="par(fig = c(0.88,0.98,0.65,1), new = TRUE)",
                "3"="par(fig = c(0.35,0.45,0.15,0.5), new = TRUE)",
                "0"="par(fig = c(0.88,0.98,0.15,0.5), new = TRUE)"
                )
   
    for (i in 1:length(Al)) {
    #foreach (comp = Alist, resCl = resClus, i = 1:length(Alist), .inorder = TRUE, .errorhandling = 'stop') %dopar% {
        comp <- Al[[i]]
        
        sapply(labGroups, function(labSel,
                                   samplesByGroup,
                                   icaSet,
                                   clus,
                                   resClus,
                                   addExpr,
                                   witExpr,
                                   witGene,
                                   breaks,
                                   colAll,
                                   colSel,
                                   titlesup,
                                   labComp,
                                   comp,
                                   indC,
                                   parFig,
                                   parMfg) {


            samples <- samplesByGroup[[labSel]]
            
            colSel <- colSel[labSel] 
            indGroup <- which(names(samplesByGroup) == labSel)
            indPlot <- (indC-1)*length(samplesByGroup)+indGroup
            
            #### -----------------------------------------------------------------------------------
            #### Test enrichment in provided clusters
            if (!missing(resClus)) {
                df <- data.frame(clus = clus)
                df$group <- rep(0,nrow(df))
                rownames(df) <- sampleNames(icaSet)
                df[samples,"group"] <- 1

                if (length(unique(df[,"group"]))>1 & length(unique(df[,"clus"]))>1) {
                    tt <- chisq.test(table(df$group, df$clus))
                    pval <- signif(tt$p.value,4)
                }
                else
                    pval <- NULL


            }

            #### -----------------------------------------------------------------------------------
            #mar1 <- par("mar")
            if (addExpr) {
                    par(mfrow = c(2,2))
                if (indPlot%%4 == 1) {
                    par(mfg=c(1,1))
                    par(mfrow = c(2,2), new = FALSE)

                }
                else
                    eval(parse(text=parMfg[as.character(indPlot%%4)]))

                par("mar"= c(5.1,4.1,5,2))
            }


            ddmax <- max(density(comp)$y)
            h<-hist(as.numeric(comp),  freq = FALSE, #xlim = c(-1,0.5),
                       breaks = breaks,
                       main = paste(if (length(keepComp)>1) paste("Component", labComp), if (!missing(resClus) & length(unique(df$clus))>1) paste("\n","Chi2 test p-val =",pval)),   
                       xlab = "contributions", yaxt = "n", if (addExpr) xaxt = "n", ylab = "Frequency",  cex.axis = 0.9, cex.lab = 1, col = colAll, ylim = if (addExpr) c(-mean(density(comp)$y)/2, ddmax+(max(density(comp)$y)/3)*2) else c(0, ddmax+(max(density(comp)$y)/3)*2)) 
            h2<-hist(as.numeric(comp[samples]), breaks = h$breaks, plot = FALSE)
            h2$density<-h2$density*length(samples)/length(comp)
            plot(h2, freq=F, add = TRUE, col = colSel, yaxt = "n", xaxt = "n")

            x = lm(h$density~h$counts)
            lab = seq(0,length(comp)/2,by=5)
            lab_coord = sapply(lab,calcCoord,x = x)
            axis(side = 2, at = lab_coord, labels = lab, cex.axis = 0.9)

            legend("topleft", legend = labSel, col = "black", pt.bg = colSel, pt.cex = 2, pch = 22, bty = "n",pt.lwd = 1.5, xjust = 1, cex = 1.1, inset = c(0.02,-0.05))

            if (funClus == "Mclust")
                plotMclust(mc = resClus, data = comp)
            else if (funClus == "kmeans") {
                cutoff <- (max(comp[rownames(subset(df,clus ==1))])+min(comp[rownames(subset(df,clus ==2))]))/2
                if (length(unique(df$clus) == 3))
                    cutoff <- c(cutoff,(max(comp[rownames(subset(df,clus ==2))])+min(comp[rownames(subset(df,clus ==3))]))/2)
                else if (length(unique(df$clus) == 4))
                    cutoff <- c(cutoff,(max(comp[rownames(subset(df,clus ==3))])+min(comp[rownames(subset(df,clus ==4))]))/2)
                abline(v=cutoff, lty = 2, col = "red")
            }

           if (!is.null(refSamples(icaSet)) | length(refSamples(icaSet))!=0) 
               points(comp[refSamples(icaSet)],y=rep(-0.04,length(refSamples(icaSet))),pch = 24, bg = "green", cex=1.35)

            if (addExpr) {

                df$comp <- comp[rownames(df)]
                df$expr <- witExpr[rownames(df)]
                df <- df[order(df$comp),]

                rm <- range(df$expr)

                
                qq1 <- quantile(as.numeric(df$expr),0.25)
                if (rm[1]==qq1)
                    qq1 <- quantile(as.numeric(df$expr),0.35)
                if (rm[1]==qq1)
                    qq1 <- quantile(as.numeric(df$expr),0.5)
                
                qq2 <- quantile(as.numeric(df$expr),0.75)
                breaksCol <- c(seq(from = rm[1], to = qq1, by = abs(qq1-rm[1])/10),
                            seq(from = qq1, to = qq2, by = (qq2-qq1)/22),
                            seq(from = qq2, to = rm[2], by = (rm[2]-qq2)/10))
                breaksCol <- unique(breaksCol)

                 expr <- df$expr
                names(expr) <- rownames(df)
                splitExpr <- split(expr,cut(expr,breaksCol))
                int2col <- maPalette(low = "blue",high = "red", mid = "yellow", k=length(splitExpr))
                if (length(int2col)<length(splitExpr))
                    int2col <- c(int2col,int2col[length(int2col)])
                names(int2col) <- names(splitExpr)

                splitExpr <- lapply(splitExpr, function(x) x <- names(x))
                splitExpr <- unlist(reverseSplit(splitExpr))
                colP <- int2col[splitExpr]
                names(colP) <- names(splitExpr)

                mm <- max(h$density)*0.60
                xax <- max(h$density)/7
                points(x = df$comp, y = rep(-xax/2,nrow(df)), pch = 15, col = colP[rownames(df)], cex = 1.5)

                eval(parse(text=parFig[as.character(indPlot%%4)]))

                ## breaks identiques a ceux dans heatmap.plus
                maColorBar(breaksCol[which(odd(1:length(breaksCol)))],
                           main = "",
                           horizontal = FALSE,
                           cex.axis = 0.8,
                           col = maPalette(low = "blue",high = "red", mid = "yellow", k=length(splitExpr)))

                # add gene names
                mtext(witGene,side=3,line=0.5,cex = 0.85) 

            }

            title(paste(if (length(keepComp)==1) paste("Component", labComp),  if (titlesup != "") paste("\n",titlesup)), outer = TRUE, line = if (length(keepComp)==1) -2 else -1) 
        },
          samplesByGroup = samplesByGroup,
          witExpr = if (!is.null(witD)) unlist(witD[i,]) else NULL,
          witGene = if (!is.null(witD)) rownames(witD)[i] else NULL,
          icaSet = icaSet,
          clus = resClus$clus[[i]],
          resClus = resClus$resClus[[i]],
          addExpr = addExpr,
          breaks = breaks,
          colAll = colAll,
          colSel = colSel,
          titlesup = titlesup,
          labComp = labelsComp[i],
          comp = comp,
          indC = i,
          parFig = parFig,
          parMfg = parMfg)       

    }

   if(!is.null(file))
    dev.off()

   return(NULL)

}

##' This function plots the positions of groups of samples formed by the variables (i.e the sample annotations) across all the components of an object of class \code{\link[MineICA:IcaSet-class]{icaSet}}.
##' For each variable level (e.g for each tumor stage) this function plots the positions of the corresponding samples (e.g the subset of samples having this tumor stage) within the histogram of the global sample contributions.
##' The plots are saved in pdf file, one file is created per variable. The pdf files are names 'variable.pdf' and save either in \code{pathPlot} if specified or the current directory.
##'
##' The plotted values are the sample contributions across the components, i.e across the columns of \code{A(icaSet)}.
##' 
##' 
##'  If argument \code{resClus} is missing, the function computes the clustering of the samples on each component (i.e on each column of \code{A(icaSet)}) using \code{funClus} and \code{nbClus}.
##' 
##' The association between the clusters and the considered sample group is tested using a chi-square test. The p-values of these tests are available in the title of each plot.
##' 
##' When \code{by="annot"} this function plots the histograms of each variable across all components, to plot the histograms for each component across variables, please use \code{by="component"}. 
##' 
##' @title Histograms of sample contributions for each annotation level
##' @param icaSet An object of class \code{IcaSet}
##' @param params A \code{MineICAParams} object
##' @param keepVar The variable labels to be considered, i.e a subset of the column labels of the pheno data of icaSet available in  (\code{varLabels(icaSet)}) 
##' @param keepComp A subset of components available in \code{indComp(icaSet)}; by default, all components are used
##' @param keepSamples A subset of samples, must be available in \code{sampleNames(icaSet)}; by default, all samples are used
##' @param pathPlot A character specifying the path where the plots will be saved
##' @param breaks The number of breaks to be used in the histograms
##' @param colSel The colour of the histogram of the group of interest, default is "red"
##' @param colAll The colour of the global histogram, default is "grey74"
##' @param resClus A list containing the outputs of function \code{clusterSamplesByComp}, which consists of sample clustering applied to matrix A of argument \code{icaSet}. If missing, the clustering is performed by the function.
##' @param funClus The clustering method to be used, either \code{"Mclust"} or \code{"kmeans"}. If \code{resClus} is not missing, equals \code{resClus$funClus}.
##' @param nbClus If \code{resClus} is missing, it provides the number of clusters to be computed by \code{funClus}, default is 2
##' @param by Either \code{"annot"} to plot the histograms of each variable across all components, or \code{"component"} to plot the histograms for each component across variables. When \code{by="annot"} one pdf file is created by variable name, while when \code{annot="component"}, one pdf file is created by component. 
##' @param typeImage The type of image to be created, either "pdf" (default) or "png". "png" is not recommended, unless there are at the most 4 histograms to be plotted, because it does not allow to deal with multiple pages of plots.
##' @param ... Additional parameters for function \code{\link{hist}}
##' @return NULL
##' @seealso \code{\link{plotPosSamplesInComp}}, \code{chisq.test}
##' @author Anne Biton
##' @export
##' @examples  \dontrun{
##' ## load an example of IcaSet
##' data(icaSetCarbayo)
##' 
##' ## Use icaSetCarbayo, look at the available annotations
##' varLabels(icaSetCarbayo)
##' 
##' ## Plot positions of samples in components according to annotations 'SEX' and 'STAGE'
##' # plots are saved in files SEX.pdf and STAGE.pdf created in the current directory 
##' plotPosAnnotInComp(icaSet=icaSetCarbayo, keepVar=c("SEX","STAGE"), keepComp=1:2,  funClus="Mclust")
##' # specifiy arg 'pathPlot' to save the pdf in another directory, but make sure it exists before
##' # specifiy arg 'by="comp"' to create one pdf file per component
##' }
##' 
plotPosAnnotInComp <- function(icaSet,
                               params,
                               keepVar = varLabels(icaSet),
                               keepComp = indComp(icaSet),
                               keepSamples = sampleNames(icaSet),
                               pathPlot = NULL,
                               breaks = 20,
                               colAll = "grey74",
                               colSel,
                               resClus,
                               funClus = c("Mclust","kmeans"),
                               nbClus = 2,
                               by = c("annot","component"),
                               typeImage=c("pdf","png", "none"),
                               ...) {


    
    icaSet <- icaSet[, keepSamples, keepComp]
    
    vardif <- setdiff(keepVar, varLabels(icaSet))
    if (length(vardif)>1)
        warning(paste("The corresponding annotations",vardif,"are not available in icaSet."))
    keepVar <- intersect(keepVar, varLabels(icaSet))

   if (missing(colSel))
       colSel = annot2Color(as.data.frame(pData(icaSet)[,keepVar]))
    
    funClus <- match.arg(funClus)
    if (missing(resClus)) {
        resClus <- clusterSamplesByComp(icaSet=icaSet, params=params, funClus=funClus, clusterOn="A", nbClus=nbClus)
    }

    by <- by[1]
    by <- match.arg(tolower(by),choices=c("annot","component"))
    typeImage <- tolower(typeImage[1])
    typeImage <- match.arg(typeImage, choices=c("pdf","png","none"))

    
    if (by=="annot") {
        sapply(keepVar,
           function(colAnnot,
                    icaSet,
                    breaks,
                    colAll,
                    colSel,
                    keepComp,
                    resClus,
                    funClus,
                    pathPlot,
                    typeImage,
                    ...){
               
               samplesByGroup <- lapply(split(pData(icaSet),pData(icaSet)[[colAnnot]]), rownames)
               samplesByGroup <- samplesByGroup[which(unlist(lapply(samplesByGroup,length))>=1)]

               if (typeImage != "none") {
                   file <- paste(pathPlot, colAnnot, ".", typeImage, sep="")
                   eval(parse(text=c(pdf="pdf(file, height = 8.267717, width = 29.7/2.54, paper = 'a4r', title=file)",png="png(file, height = 8.267717, width = 29.7/2.54, units='in', res=300)")[typeImage]))
               }

               plotPosSamplesInComp(samplesByGroup = samplesByGroup,
                                    icaSet = icaSet,
                                    breaks = breaks,
                                    colAll = colAll,
                                    colSel = colSel[names(samplesByGroup)],
                                    labGroups = names(samplesByGroup),
                                    titlesup = colAnnot,
                                    resClus = resClus,
                                    funClus = funClus,
                                    ...)
               if (typeImage != "none") 
                   graphics.off()

           },
           icaSet = icaSet,
           keepComp = keepComp,
           breaks = breaks,
           colAll = colAll,
           colSel = colSel,
           resClus = resClus,
           funClus = funClus,
           pathPlot = pathPlot,
           typeImage = typeImage,
           ...
           )

    }
    else {
        
        sapply(as.list(1:length(keepComp)),
               function(indC,
                        keepVar,
                        icaSet,
                        breaks,
                        colAll,
                        colSel,
                        resClus,
                        funClus,
                        pathPlot,
                        indCompAll,
                        typeImage,
                        ...){
                   
               if (typeImage != "none") {
                   filePlotComp <- paste(pathPlot,indCompAll[indC],".",typeImage,sep="")
                   eval(parse(text=c(pdf="pdf(filePlotComp, height = 8.267717, width = 29.7/2.54, paper = 'a4r', title=filePlotComp)",png="png(filePlotComp, height = 8.267717, width = 29.7/2.54, units='in', res=300)")[typeImage]))
               }

                   sapply(keepVar,
                          function(colAnnot,
                                   icaSet,
                                   breaks,
                                   colAll,
                                   colSel,
                                   indC,
                                   resClus,
                                   funClus,
                                   ...){

                              samplesByGroup <- lapply(split(pData(icaSet),pData(icaSet)[[colAnnot]]), rownames)
                              samplesByGroup <- samplesByGroup[which(unlist(lapply(samplesByGroup,length))>=1)]
                              
                              plotPosSamplesInComp(samplesByGroup = samplesByGroup,
                                                   icaSet = icaSet,
                                                   keepComp=indCompAll[indC],
                                                   breaks = breaks,
                                                   colAll = colAll,
                                                   colSel = colSel[names(samplesByGroup)],
                                                   labGroups = names(samplesByGroup),
                                                   titlesup = colAnnot,
                                                   resClus = resClus,
                                                   funClus = funClus,
                                                   ...)
                          },
                          icaSet = icaSet,
                          indC = indC,
                          breaks = breaks,
                          colAll = colAll,
                          colSel = colSel,
                          resClus = resClus,
                          funClus = funClus,
                          ...
                          )
               if (typeImage != "none") 
                   graphics.off()
                   
               },
               keepVar = keepVar,
               icaSet = icaSet,
               breaks = breaks,
               colAll = colAll,
               colSel = colSel,
               resClus = resClus,
               funClus = funClus,
               pathPlot = pathPlot,
               indCompAll = indComp(icaSet),
               typeImage = typeImage,
               ...
               )
    }

}


