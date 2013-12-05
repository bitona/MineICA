##' This function tests if the groups of samples formed by the variables are differently distributed on the components, in terms of contribution value (i.e of values in matrix \code{A(icaSet)}). The distribution of the samples  on the components are represented using either density plots of boxplots. It is possible to restrict the tests and the plots to a subset of samples and/or components.
##' 
##' This function writes an HTML file containing the results of the tests as a an array of dimensions 'variables * components' containing the p-values of the tests. When a p-value is considered as significant according to the threshold \code{cutoff}, it is written in bold and filled with a link pointing to the corresponding plot.
##' One image is created by plot and located into the sub-directory "plots/" of \code{path}. Each image is named by index-of-component_var.png.
##' Wilcoxon or Kruskal-Wallis tests are performed depending on the number of groups of interest in the considered variable (argument \code{keepLev}).
##' @title Tests association between qualitative variables and components.
##' @param params An object of class \code{\link[MineICA:MineICAParams-class]{MineICAParams}} providing the parameters of the analysis.
##' @param icaSet An object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}.
##' @param keepVar The variable labels to be considered, must be a subset of  \code{varLabels(icaSet)}.
##' @param keepComp A subset of components, must be included in \code{indComp(icaSet)}. By default, all components are used.
##' @param keepSamples A subset of samples, must be included in \code{sampleNames(icaSet)}. By default, all samples are used.
##' @param adjustBy The way the p-values of the Wilcoxon and Kruskal-Wallis tests should be corrected for multiple testing: \code{"none"} if no p-value correction has to be done, \code{"component"} if the p-values have to be corrected by component, \code{"variable"} if the p-values have to be corrected by variable
##' @param method The correction method, see \code{\link{p.adjust}} for details, default is \code{"BH"} for Benjamini & Hochberg.
##' @param doPlot If TRUE (default), the plots are done, else only tests are performed.
##' @param addPoints If TRUE, points are superimposed on the boxplot.
##' @param typePlot The type of plot, either \code{"density"} or \code{"boxplot"}.
##' @param onlySign If TRUE (default), only the significant results are plotted.
##' @param cutoff A threshold p-value for statistical significance.
##' @param colours A vector of colours indexed by the variable levels, if missing the colours are automatically generated using \code{\link{annot2Color}}.
##' @param path A directory _within resPath(params)_ where the files containing the plots and the p-value results will be located. Default is \code{"qualVarAnalysis/"}.
##' @param typeImage The type of image file to be used.
##' @param filename The name of the HTML file containing the p-values of the tests, if NULL no file is created.
##' @return Returns A data.frame of dimensions 'components x variables' containing the p-values of the non-parametric tests (Wilcoxon or Kruskal-Wallis tests) wich test if the samples groups defined by each variable are differently distributed on the components. 
##' @seealso , \code{\link{qualVarAnalysis}}, \code{\link{p.adjust}}, \code{link{writeHtmlResTestsByAnnot}}, \code{wilcox.test}, \code{kruskal.test}
##' @author Anne Biton
##' @export
##' @examples 
##' ## load an example of IcaSet 
##' data(icaSetCarbayo)
##' 
##' ## build MineICAParams object
##' params <- buildMineICAParams(resPath="carbayo/")
##' 
##' ## Define the directory containing the results
##' dir <- paste(resPath(params), "comp2annot/", sep="")
##' 
##' ## Run tests, make no adjustment of the p-values,
##' # for variable grade and components 1 and 2,
##' # and plot boxplots when 'doPlot=TRUE'.
##' qualVarAnalysis(params=params, icaSet=icaSetCarbayo, adjustBy="none", typePlot="boxplot",
##'                 keepVar="GRADE", keepComp=1:2, path=dir, doPlot=FALSE)
##' 
##' 
qualVarAnalysis <- function(params, icaSet, keepVar, keepComp=indComp(icaSet), keepSamples=sampleNames(icaSet), adjustBy = c("none", "component", "variable"), method = "BH", doPlot = TRUE, typePlot = "density", addPoints=FALSE, onlySign = TRUE, cutoff = params["pvalCutoff"], colours=annot2col(params), path="qualVarAnalysis/", filename = "qualVar", typeImage = "png") {

    res <- an <- NULL
    
    icaSet <- icaSet[, keepSamples, keepComp]
           
    adjustBy <- match.arg(adjustBy)

    if (!missing(keepVar)) {
        keepVar <- intersect(keepVar,varLabels(icaSet))
        if (length(keepVar)==0)
            stop("Arg 'keepVar' is not included in varLabels(icaSet).")
    }
    else
        keepVar <- varLabels(icaSet)
        
    if (length(keepVar)==1) {
        annot <- data.frame(pData(icaSet)[,keepVar], row.names=sampleNames(icaSet))
        colnames(annot) <- keepVar
    }
    else 
        annot <- pData(icaSet)[,keepVar]

    A <- A(icaSet)

    Alist <- Alist(icaSet)
        
    if (length(refSamples(icaSet))>0) {
        annot_ref <- subset(annot, rownames(annot) %in% refSamples(icaSet))
        annot <- subset(annot, !(rownames(annot) %in% refSamples(icaSet)))
    }
    else
        annot_ref <- NULL
    
        
    resTests <- lapply(colnames(annot), 
                            function(an, annot, Alist) {
                                
                                levRepr <- table(annot[[an]])
                                levNoRepr <- names(levRepr)[levRepr < 2]
                                if (length(levNoRepr)>0)
                                    warning(paste("The annotation levels '",paste(levNoRepr,collapse=", "), "' map to less than 2 samples, they have not been taken into account during the tests.",sep =""))
                                if (!is.factor(annot[[an]]))
                                    keepLev <- as.character(unique(annot[[an]]))
                                else
                                    keepLev <- levels(annot[[an]])
                                    
                                keepLev <- keepLev[!(keepLev %in% levNoRepr)]
                                keepLev <- keepLev[!is.na(keepLev)]

                                if (length(keepLev) <= 1)  {
                                    warning(paste("'All samples are in the same group after removing levels corresponding to less than 2 samples in or NA in variable ",keepVar, ".", sep =""))
                                    return(list(pvalByComp = rep(NA,length(Alist)), keepLev = list(keeplev = NA,annot.bis=NA)))
                                }
                                
                                annot.bis <- subset(annot,annot[[an]] %in% keepLev)
                                annot.bis[[an]] = factor(annot.bis[[an]], levels = keepLev)
                                annot.bis <- subset(annot.bis,!is.na(as.character(annot.bis[[an]])))
                                keepLevlist <- list(keepLev = keepLev, annot.bis = annot.bis)
                                pvalByComp <- lapply(Alist,
                                                     function(comp, annot.bis, keepLev, an) {
                                                         comp <- comp[rownames(annot.bis)]
                                                         if (length(keepLev) == 2) 
                                                             resTests <- wilcox.test(comp~eval(as.name(an)), data = annot.bis, na.action = "na.omit") #comp~eval(as.name(keepVar))
                                                         else 
                                                             resTests <- kruskal.test(comp~eval(as.name(an)), data = annot.bis, na.action = "na.omit")

                                                         return(resTests$p.value)

                                                     }
                                                     , annot.bis = annot.bis
                                                     , keepLev = keepLev
                                                     , an = an
                                                     )
                                pvalByComp <- unlist(pvalByComp)
                                names(pvalByComp) <- indComp(icaSet)
                                return(list(pvalByComp = pvalByComp, keepLev = keepLevlist))
                            }
                             , annot = annot
                             , Alist = Alist
                            )
    keepLevByAnnot <- lapply(resTests, function(x) return(x$keepLev))
    names(keepLevByAnnot) <- colnames(annot)
    resTests <- lapply(resTests, function(x) return(x$pvalByComp))
        
    resTests <- t(as.data.frame(resTests, check.names = FALSE))

    deleteAnnot <- which(unlist(apply(resTests,1,function(x) all(is.na(x)))))
    

    if(length(deleteAnnot)>0) {
        keepLevByAnnot <- keepLevByAnnot[-deleteAnnot]
        keepVar <- keepVar[-deleteAnnot]
        resTests <- resTests[-deleteAnnot,]

        if(length(keepVar)==1) {
            resTests <- t(data.frame(resTests))
            annot <- as.data.frame(annot[,-deleteAnnot])
            rownames(annot) <- sampleNames(icaSet)
            colnames(annot) <- rownames(resTests) <- keepVar
            
        }
        else
            annot <- annot[,-deleteAnnot]
            
        pData(icaSet) <- annot
    }
    
    if (adjustBy != "none") {

        if (adjustBy == "component") {
            resTests <- as.data.frame(apply(resTests, 2, p.adjust, method = method), stringsAsFactors = FALSE, check.names = FALSE)
        }
        else if (adjustBy == "variable") {
            resTests <- as.data.frame(t(apply(resTests, 1, p.adjust, method = method)), stringsAsFactors = FALSE, check.names = FALSE)         
        }
    }

    rownames(resTests) <- colnames(annot)
    colnames(resTests) <- indComp(icaSet)

    if (missing(path))
        path <- paste(resPath(params),"qualVarAnalysisOnA",if (adjustBy=="none") "_noAdjust" else paste("_adjustBy",adjustBy,sep=""),"/",sep="")
    else
        path <- gsub(paste(resPath(params),gsub(path,pattern=resPath(params),replacement=""),"/",sep=""),pattern="//",replacement="/")
    
    pathplot <- paste(path, "plots/",sep="")


    if (doPlot) {

        system(paste("rm -r ", path,"plots/",sep=""),ignore.stderr=TRUE)

        system(paste("mkdir", path, pathplot),ignore.stderr=TRUE)
        
        if (onlySign) {
            whichAnnotSign <- unlist(apply(resTests,1,function(x,cutoff) {x <- as.numeric(x); length(x[x<=cutoff & !is.na(x)])}, cutoff = cutoff))
            whichAnnotSign <- which(whichAnnotSign > 0)
        }
        else 
            whichAnnotSign <- 1:nrow(resTests)
        
        
        if (missing(colours) | length(colours)==0)
              colours <- annot2Color(annot)

        
        if (typePlot == "density") {
            trace_globalExpression <- FALSE
            trace_groupExpression <- TRUE
        }
        else if (typePlot == "boxplot") {
            trace_globalExpression <- TRUE
            trace_groupExpression <- FALSE
        }
        
        foreach(keepLev = keepLevByAnnot[whichAnnotSign], an = names(keepLevByAnnot)[whichAnnotSign]) %dopar% {

            keepLev <- keepLevByAnnot[[an]]
            print(paste("Plot distribution of samples on components according to variable",an))
            annot.bis <- keepLev$annot.bis

            if (onlySign) {
                whichCompSign <- which(resTests[an,] <= cutoff)
                compSign <- indComp(icaSet)[whichCompSign]
            }
            else {
                whichCompSign <- 1:ncol(resTests)
                compSign <- indComp(icaSet)
            }
            

            sapply(1:length(whichCompSign), function(i, an, keepLev, Alist, whichCompSign, compSign, icaSet, colours, annot_ref, annot.bis, typePlot, trace_globalExpression, trace_groupExpression) {
                indC <- whichCompSign[i]
                comp <- Alist[[indC]]
                numC <- compSign[i]
                labComp <- compNames(icaSet)[indC]
                
               print(paste("Comp",numC))
            
                if (length(refSamples(icaSet)) > 0)
                    annot_ref$comp = comp[rownames(annot_ref)]

                annot.bis$comp <- comp[rownames(annot.bis)]
                
 
                g <- plotDens2classInComp_plotOnly(
                                                   annot = annot.bis,
                                                   colAnnot = an,
                                                   global = comp,
                                                   comp.label = labComp,
                                                   if (!missing(colours)) colours = colours,
                                                   legend.title=NULL,
                                                   pval = resTests[an,indC],
                                                   test = if (length(keepLev$keepLev)>2) "Kruskal-Wallis test" else "Wilcoxon test" ,
                                                   title.add = NULL,
                                                   data_ref = annot_ref,
                                                   geneExpr = if (length(witGenes(icaSet))>0) unlist(datByGene(icaSet)[witGenes(icaSet)[indC],]) else NULL,
                                                   geneRef = if (length(witGenes(icaSet))>0) witGenes(icaSet)[indC] else NULL,
                                                   keepLev = keepLev$keepLev,
                                                   typePlot = typePlot,
                                                   addPoints = addPoints,
                                                   trace_globalExpression = trace_globalExpression,
                                                   trace_groupExpression = trace_groupExpression
                                                   )
                ggsave(plot=g,filename=paste(pathplot,numC,"_",an,".",typeImage,sep=""),height = 5, width = 6)

            }
            , an=an, keepLev=keepLev, Alist=Alist, whichCompSign=whichCompSign, compSign=compSign, icaSet=icaSet,colours=colours, annot_ref=annot_ref, annot.bis=annot.bis, typePlot=typePlot, trace_globalExpression=trace_globalExpression, trace_groupExpression=trace_groupExpression
            )

        }
    }
    
        if (!is.null(filename)) {
            writeHtmlResTestsByAnnot(icaSet = icaSet,
                               params = params,
                               res = as.data.frame(t(signif(resTests,5)), check.names=FALSE, stringsAsFactors=FALSE),
                               path = path,
                               pathplot = pathplot,
                               filename = filename,
                               cutoff = cutoff,
                               typeImage = typeImage,
                               caption = paste("<b>The contributions of pre-established groups of samples were compared using Wilcoxon rank-sum test for 2-classes comparison and Kruskal-Wallis for the comparison of more than 2 classes.</b> This matrix provides the ", if (adjustBy != "none") paste("adjusted p-values (by ", adjustBy, ", with method ", method, ") of these tests.",sep = "") else "p-values of these tests. <br> Click on the significant values to visualize the densities of the different sub-groups. "),
                               onlySign = onlySign,
                               keepVar=keepVar)
        }


    
    return(resTests)
}


##' This function tests if numeric variables are correlated with components.
##' 
##' This function writes an HTML file containing the correlation values and test p-values as a an array of dimensions 'variables * components' containing the p-values of the tests. When a p-value is considered as significant according to the threshold \code{cutoff}, it is written in bold and filled with a link pointing to the corresponding plot.
##' One image is created by plot and located into the sub-directory "plots/" of \code{path}. Each image is named by index-of-component_var.png.
##' @title Correlation between variables and components.
##' @param params An object of class \code{\link[MineICA:MineICAParams-class]{MineICAParams}} providing the parameters of the analysis.
##' @param icaSet An object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}.
##' @param keepVar The variable labels to be considered, must be a subset of  \code{varLabels(icaSet)}.
##' @param keepComp A subset of components, must be included in \code{indComp(icaSet)}. By default, all components are used.
##' @param keepSamples A subset of samples, must be included in \code{sampleNames(icaSet)}. By default, all samples are used.
##' @param adjustBy The way the p-values of the Wilcoxon and Kruskal-Wallis tests should be corrected for multiple testing: \code{"none"} if no p-value correction has to be done, \code{"component"} if the p-values have to be corrected by component, \code{"variable"} if the p-values have to be corrected by variable
##' @param method The correction method, see \code{\link{p.adjust}} for details, default is \code{"BH"} for Benjamini & Hochberg.
##' @param doPlot If TRUE (default), the plots are done, else only tests are performed.
##' @param onlySign If TRUE (default), only the significant results are plotted.
##' @param cutoff A threshold p-value for statistical significance.
##' @param cutoffOn The value the cutoff is applied to, either "cor" for correlation or "pval" for p-value
##' @param typeCor the type of correlation to be used, one of \code{c("pearson","spearman","kendall")}.
##' @param colours A vector of colours indexed by the variable levels, if missing the colours are automatically generated using \code{\link{annot2Color}}.
##' @param path A directory  _within resPath(params)_ where the files containing the plots and the p-value results will be located. Default is \code{"quantVarAnalysis/"}.
##' @param typeImage The type of image file to be used.
##' @param filename The name of the HTML file containing the p-values of the tests, if NULL no file is created.
##' @return Returns A data.frame of dimensions 'components x variables' containing the p-values of the non-parametric tests (Wilcoxon or Kruskal-Wallis tests) wich test if the samples groups defined by each variable are differently distributed on the components. 
##' @seealso \code{\link{qualVarAnalysis}}, \code{\link{p.adjust}}, \code{link{writeHtmlResTestsByAnnot}}, \code{code}
##' @author Anne Biton
##' @export
##' @examples 
##' ## load an example of IcaSet 
##' data(icaSetCarbayo)
##' 
##' # build MineICAParams object
##' params <- buildMineICAParams(resPath="carbayo/")
##' 
##' # Define the directory containing the results
##' dir <- paste(resPath(params), "comp2annottest/", sep="")
##'
##' # Check which variables are numeric looking at the pheno data, here only one  -> AGE
##' # pData(icaSetCarbayo)
##' 
##' ## Perform pearson correlation tests and plots association corresponding
##' # to correlation values larger than 0.2
##' quantVarAnalysis(params=params, icaSet=icaSetCarbayo, keepVar="AGE", keepComp=1:2,
##'                  adjustBy="none", path=dir, cutoff=0.2, cutoffOn="cor")
##'
##' \dontrun{
##' ## Perform Spearman correlation tests and do scatter plots for all pairs
##' quantVarAnalysis(params=params, icaSet=icaSetCarbayo, keepVar="AGE", adjustBy="none", path=dir,
##'                  cutoff=0.1, cutoffOn="cor", typeCor="spearman", onlySign=FALSE)
##' 
##' ## Perform pearson correlation tests and plots association corresponding
##' # to p-values lower than 0.05 when 'doPlot=TRUE'
##' quantVarAnalysis(params=params, icaSet=icaSetCarbayo, keepVar="AGE", adjustBy="none", path=dir,
##'                  cutoff=0.05, cutoffOn="pval", doPlot=FALSE)
##' }
##' 
quantVarAnalysis <- function(params, icaSet, keepVar, keepComp=indComp(icaSet), keepSamples=sampleNames(icaSet), adjustBy = c("none", "component", "variable"), method = "BH", typeCor = "pearson", doPlot = TRUE, onlySign = TRUE, cutoff = 0.4, cutoffOn = c("cor","pval"), colours, path="quantVarAnalysis/", filename = "quantVar", typeImage = "png") {

    res <- an <- NULL
    
    icaSet <- icaSet[,keepSamples,keepComp]

    cutoffOn <- match.arg(cutoffOn)
    cutoffDir <- unlist(c(pval="<=",cor=">=")[cutoffOn])
    adjustBy <- match.arg(adjustBy)

    if (!missing(keepVar)) {
        keepVar <- intersect(keepVar,varLabels(icaSet))
        if (length(keepVar)==0)
            stop("Arg 'keepVar' is not included in varLabels(icaSet).")
    }
    else
        keepVar <- varLabels(icaSet)

        
    if (length(refSamples(icaSet))>0) {
        annot_ref <- subset(pData(icaSet), sampleNames(icaSet) %in% refSamples(icaSet))
        annot <- subset(pData(icaSet), !(rownames(pData(icaSet)) %in% refSamples(icaSet)))
    }
    else
        annot_ref <- NULL

    if (length(keepVar)==1) {
        annot <- data.frame(pData(icaSet)[,keepVar], row.names=sampleNames(icaSet))
        colnames(annot) <- keepVar
    }
    else 
        annot <- pData(icaSet)[,keepVar]

    annot <- as.data.frame(apply(annot,2,as.numeric),stringsAsFactors=FALSE, row.names=sampleNames(icaSet))

    A <- A(icaSet)

    Alist <- Alist(icaSet)

        
    if (length(refSamples(icaSet))>0) {
        sampleRef <- refSamples(icaSet)
    }
    else
        sampleRef <- NULL
            
    resTests <- lapply(colnames(annot), 
                            function(an, annot, Alist, typeCor) {
                                
                                anval <- as.numeric(annot[[an]])
                                names(anval) <- rownames(annot)
                                if (all(is.na(anval)))
                                    stop(paste(an,"is not numeric (only includes NA after as.numeric)."))
                                    
                                pvalByComp <- lapply(Alist,
                                                     function(comp, annot.bis, an, typeCor) {
                                                         comp <- comp[names(annot.bis)]
                                                         testcor <- cor.test(comp,annot.bis, use = "complete.obs", method=typeCor)
                                                         return(list(cor=testcor$estimate,corpval=testcor$p.value))
                                                     }
                                                     , annot.bis = anval
                                                     , an = an
                                                     , typeCor = typeCor
                                                     )
                                corByComp <- unlist(lapply(pvalByComp,function(x) return(x$cor)))
                                pvalByComp <- unlist(lapply(pvalByComp,function(x) return(x$corpval)))
                                
                                names(pvalByComp) <- names(corByComp) <- indComp(icaSet)
                                
                                return(list(pvalByComp = pvalByComp, corByComp = corByComp))
                        }
                       , annot = annot
                       , Alist = Alist
                       , typeCor = typeCor
                       )
    resCor <- lapply(resTests, function(x) return(x$corByComp))        
    resCor <- t(as.data.frame(resCor, check.names = FALSE))
    
    resTests <- lapply(resTests, function(x) return(x$pvalByComp))        
    resTests <- t(as.data.frame(resTests, check.names = FALSE))
    
    if (adjustBy != "none") {

        ## adjust by column
        if (adjustBy == "component") {
            resTests <- resTests <- as.data.frame(matrix(apply(resTests, 2, p.adjust, method = method), nrow=ncol(annot)), stringsAsFactors = FALSE, check.names = FALSE)
        }
        else if (adjustBy == "variable") {
            resTests <- as.data.frame(matrix(t(apply(resTests, 1, p.adjust, method = method)),nrow=ncol(annot)), stringsAsFactors = FALSE, check.names = FALSE)

        }
    }

    rownames(resTests) <- rownames(resCor) <- colnames(annot)
    colnames(resTests) <- colnames(resCor) <- indComp(icaSet)

    ressave <- resTests
    ressavecor <- resCor
    if (cutoffOn == "cor") {
        resTests <- resCor
        resCor <- ressave
    }

    if (missing(path))
        path <- paste(resPath(params),"quantVarAnalysisOnA",if (adjustBy=="none") "_noAdjust" else paste("_adjustBy",adjustBy,sep=""),"/",sep="")
    else
        path <- gsub(paste(resPath(params),gsub(path,pattern=resPath(params),replacement=""),"/",sep=""),pattern="//",replacement="/")
    
    pathplot <- paste(path, "plots/",sep="")


    if (doPlot) {

        system(paste("rm -r ", path,"plots/",sep=""),ignore.stderr=TRUE)
        system(paste("mkdir", path, pathplot),ignore.stderr=TRUE)
        
        ## select rows that have at least one significant result
        if (onlySign) {
                whichAnnotSign <- unlist(apply(resTests,1,function(x,cutoff) {x <- as.numeric(x); length(x[eval(parse(text=paste("abs(x)",cutoffDir,cutoff,sep=""))) & !is.na(x)])}, cutoff = cutoff))
                whichAnnotSign <- which(whichAnnotSign > 0)
        }
        else {
            whichAnnotSign <- 1:nrow(resTests)
        }
        

        if (missing(colours))
              colours <- annot2Color(annot)

                
        foreach(an = keepVar[whichAnnotSign]) %dopar% {
            print(paste("Scatter plot of samples contributions vs variable",an))
            annot.bis <- annot[[an]]

            if (onlySign)
                whichCompSign <- which(eval(parse(text=paste("abs(resTests[an,])",cutoffDir,cutoff,sep=""))))
            else
                whichCompSign <- 1:ncol(resTests)

            compSign <- indComp(icaSet)[whichCompSign]
            
            sapply(1:length(whichCompSign), function(i, an, Alist, whichCompSign, compSign, icaSet, colours, annot_ref, annot.bis, typeCor) {
                indC <- whichCompSign[i]
                comp <- Alist[[indC]]
                numC <- compSign[i]
                labComp <- compNames(icaSet)[indC]
                
               print(paste("Comp",numC))
                            
                g <- plotNumVarComp(annot = annot.bis,
                                    comp = comp,
                                    title=if (cutoffOn=="cor") paste(an, " vs comp ", labComp, "\n", typeCor, " cor=",signif(resTests[an,indC],3),", pval=",signif(resCor[an,indC],4), sep="") else paste(an, " vs comp ", labComp, "\n", typeCor, " cor=",signif(resCor[an,indC],3),", pval=",signif(resTests[an,indC],4), sep=""),
                                    sampleRef = sampleRef,
                                    geneExpr = if (length(witGenes(icaSet))>0) unlist(datByGene(icaSet)[witGenes(icaSet)[indC],]) else NULL,
                                    geneRef = if (length(witGenes(icaSet))>0) witGenes(icaSet)[indC] else NULL,
                                    ylab = an,
                                    xlab = paste("comp", labComp)
                                 )
                ggsave(plot=g,filename=paste(pathplot,numC,"_",an,".",typeImage,sep=""),height = 5, width = 6)

            }
            , an=an, Alist=Alist, whichCompSign=whichCompSign, compSign=compSign, icaSet=icaSet,colours=colours, annot_ref=annot_ref, annot.bis=annot.bis, typeCor=typeCor
            )

        }
    }
    
        if (!is.null(filename)) {
            
            writeHtmlResTestsByAnnot(icaSet = icaSet,
                               params = params,
                               res = as.data.frame(t(signif(resTests,3))),
                               res2 = as.data.frame(t(signif(resCor,3))) ,
                               nameres = cutoffOn,
                               nameres2 = c("cor","pval")[!c("cor","pval") %in% cutoffOn],
                               path = path,
                               pathplot = pathplot,
                               filename = filename,
                               cutoff = cutoff,
                               cutoffDir = cutoffDir,
                               typeImage = typeImage,
                               caption = paste("<b>Association between numeric variables and sample contributions.</b> This matrix contains the <b>", typeCor, " </b> correlations and ", if (adjustBy != "none") paste("adjusted p-values (by ", adjustBy, ", with method ", method, ") of the corresponding tests.",sep = "") else "p-values of the tests. <br>", "Click on the", if (cutoffOn=="cor") "correlations" else "p-values", if (onlySign) "in bold", "to visualize the scatter plots."),
                               onlySign = onlySign,
                               keepVar=keepVar)
        }


    
    return(list(pval=ressave,cor=ressavecor))
}


plotNumVarComp <- function (annot,
                         comp,
                         title,
                         sampleRef=NULL,
                         geneExpr,
                         geneRef,
                         file,
                         xlab,
                         ylab) {
    

    xvar <- yvar <- ref <- expr <- NULL
    dat <- data.frame(xvar=comp, yvar=annot, stringsAsFactors=FALSE, row.names=names(comp))

    if (!is.null(sampleRef)) {
        dat$ref <- as.character(rep("nonRef",nrow(dat)))
        dat[sampleRef,"ref"] <- "ref"
        dat$ref <- as.factor(dat$ref)
        g <- ggplot(dat)  +
            geom_point(aes(x=xvar, y=yvar, shape=ref), size=2.5) +    
                geom_smooth(aes(x=xvar, y=yvar),method=lm,se=FALSE) +
                    scale_colour_manual(values=c("nonRef"="black","ref"="green")) +
                        scale_shape_manual(values=c("nonRef"=20,"ref"=17)) 
    }
    else {
    g <- ggplot(dat) +
                geom_point(aes(x=xvar, y=yvar), size=2.5, shape=20) 
                
    }
    g <- g + geom_smooth(aes(x=xvar, y=yvar),method=lm,se=FALSE) + xlab(xlab) + ylab(ylab) + ggtitle(title)

    
    if (!is.null(geneExpr)) {
        yy <-  min(dat$yvar)-(abs(max(dat$yvar)-min(dat$yvar)))/15
        df <- data.frame(expr = geneExpr,  xvar = as.numeric(comp), yvar=rep(yy,length(comp)))
        genename <- geneRef
        g <- g + geom_point(data=df, aes(x = xvar, y = yvar, colour = expr),size = 3, shape = 15) + scale_colour_gradientn(name = genename,colours = maPalette(low = "blue",high = "red", mid = "yellow", k=50))

    }
        

    return(g)
    
                    
    

}


##' Compare the sample contributions according to their annotation level across the components.
##' 
##' Wilcoxon or Kruskal-Wallis tests are performed depending on the number of levels in the considered annotation.
##' @title Comparison of distributions of sample groups 
##' @param A A matrix of dimensions 'samples x components' containing the sample contributions
##' @param annot A matrix of dimensions 'samples x variables' containing the sample annotations
##' @param colAnnot The name of the column of \code{annot} to be considered
##' @return A vector of p-values
##' @author Anne Biton
##' @seealso \code{wilcox.test}, \code{kruskal.test}
##' @keywords internal
wilcoxOrKruskalOnA <-
    function (A,
              colAnnot,
              annot) {

        comp <- NULL
        A <- A[rownames(annot),]
        
	res_tests <- foreach(comp=as.list(A),.combine = c, .errorhandling = "stop") %dopar% {
 			annotComp <- data.frame(comp=comp) #interest=as.factor(annot[[colAnnot]]),
                        annotComp[[colAnnot]] <- as.factor(annot[[colAnnot]])
 			if (length(levels(annotComp[[colAnnot]])) == 2) 
				res.test <- wilcox.test(comp~eval(as.name(colAnnot)), data = annotComp, na.action = "na.omit") #comp~eval(as.name(colAnnot))
			else 
				res.test <- kruskal.test(comp~eval(as.name(colAnnot)), data = annotComp, na.action = "na.omit")
			return(res.test$p.value)			
	}

	return(unlist(res_tests))

}

##' Computes the relative path between two imbricated paths
##'
##' \code{path1} and \code{path2} must be imbricated.
##' @title Relative path 
##' @param path1 The first path
##' @param path2 The second path
##' @return  The relative path between path1 and path2
##' @author Anne Biton
##' @examples
##' path1 <- "home/lulu/res/gene2comp/"
##' path2 <- "home/lulu/res/comp2annot/invasive/"
##' relativePath(path1,path2)
##' @export
relativePath <- function(path1, path2) {
    path1 <- gsub(path1, pattern = "//", replacement = "/")
    path2 <- gsub(path2, pattern = "//", replacement = "/")
    if (substr(path1,start=nchar(path1),stop=nchar(path1)) != "/")
        path1 = paste(path1,"/",sep="")
    if (substr(path2,start=nchar(path2),stop=nchar(path2)) != "/")
        path2 = paste(path2,"/",sep="")
    if (nchar(path1) > nchar(path2)) {
        path <- path1
        path1 <- path2
        path2 <- path
    }
    commonpath <- lcPrefix(c(path1,path2), ignore.case=FALSE)
    path2 <- gsub(path2, pattern = commonpath, replacement = "")
    relpath <- paste(rep("../",length(gregexpr(text = path2, pattern= "/")[[1]])), collapse = "")

    return(relpath)
}

##' This internal function creates an HTML file containing a table of dimensions 'variables x components' with p-values. When a p-value is considered as significant according to the threshold \code{cutoff}, it is written in bold and filled with a link pointing to the corresponding plot.
##' These plots are contained in images located in the path \code{pathplot}. To be identified by the function, the file syntax of each image file must be "index-of-component_colAnnot.typeImage". 
##' 
##' If argument \code{onlySign} is TRUE, then only links to plots that are significant according to the given threshold are provided. 
##'
##' When \code{res2} is not missing, the values contained in \code{res2} are pasted to the values contained in res in the output array. \code{nameres} and \code{nameres2} are used such as every element in the ouput array contains two indexed values: \code{nameres}=x, \code{nameres2}=y.  
##' 
##' @title Tests if groups of samples are differently distributed on the components according and do the corresponding plots.
##' @param params An object of class \code{\link[MineICA:MineICAParams-class]{MineICAParams}} containing the parameters of the analysis
##' @param icaSet An object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}
##' @param res A matrix or data.frame of dimension 'components x variables' containing numeric values that quantify the association of the components with sample variables (e.g p-values, FDR, correlation values). This is the matrix used to select the significant results according to \code{cutoff} and \code{cutoffDir}.
##' @param res2 A matrix or data.frame of dimension 'components x variables' containing numeric values that quantify the association of the components with sample annotations (e.g p-values, FDR, correlation values). It is only used as an additional result displayed in the output.
##' @param nameres Name of the values contained in  \code{res}, default is "p"
##' @param nameres2 Name of the values contained in  \code{res2}, default is "cor"
##' @param onlySign If TRUE (default), only the significant results are plotted
##' @param cutoff The threshold p-value for statistical significance
##' @param path A directory for the HTML file containing the p-value results  
##' @param pathplot A directory for the plots  
##' @param filename The name of the file where the results will be displayed in format HTML, if NULL no file is created 
##' @param typeImage The type of image file where each plot is saved
##' @param caption The title of the HTML table
##' @param cutoffDir The direction to be used with the cutoff: \code{"inf"} for "<=" and \code{"sup"} for ">="
##' @param keepVar The variable labels to be considered, i.e a subset of the variables of icaSet available in \code{varLabels(icaSet)}.
##' @return Returns a data.frame of dimensions 'components x variables' containing the p-values of the non-parametric tests (Wilcoxon or Kruskal-Wallis tests) wich test if the samples groups defined by each variable are differently distributed on the components 
##' @seealso \code{\link{p.adjust}}, \code{\link{qualVarAnalysis}}, \code{\link{quantVarAnalysis}}
##' @author Anne Biton
##' @export
##' @keywords internal
writeHtmlResTestsByAnnot <- function(params, icaSet, res, res2, nameres="p", nameres2="cor", onlySign = TRUE, cutoff = params["pvalCutoff"], cutoffDir = c("<=",">="), path, pathplot = "plots/", filename = NULL, typeImage = "png", caption = "", keepVar) {

    typeImage <- gsub(typeImage,pattern=".",replacement="",fixed=TRUE)
    
    if (!(is.null(path)) && path != "" && substr(path,start=nchar(path),stop=nchar(path)) != "/")
        path <- paste(path,"/",sep="")

    A <- A(icaSet)
    pathComp <- gsub(params["genesPath"],pattern = params["resPath"],replacement = "")

    if (!missing(keepVar)) {
        keepVar <- intersect(keepVar, varLabels(icaSet))
        if (length(keepVar)==0)
            stop("Arg 'keepVar' is not included in varLabels(icaSet).")
    }
    else
        keepVar <- varLabels(icaSet)
    
    
    cutoffDir <- cutoffDir[1]
    cutoffDir <- match.arg(cutoffDir)
    
    ## relative path
    relPath <- relativePath(path1 = path, path2 = genesPath(params))
    system(paste("mkdir", path, pathplot),ignore.stderr=TRUE)
    allimages <- list.files(path = paste(path,gsub(pattern=path,x=pathplot,replacement=""),sep = ""), pattern = paste(".",typeImage,sep=""), all.files = TRUE,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = FALSE)

    comp <- rownames(res)

    if (length(keepVar)==1) {
        annot <- data.frame(pData(icaSet)[,keepVar], row.names=sampleNames(icaSet))
        colnames(annot) <- keepVar
    }
    else 
        annot <- pData(icaSet)[,keepVar]

    saverownames <- colnames(res)
    resnum <- t(res)
    res <- sapply(1:nrow(res),
                     function(comp,res,pathplot,cutoff,typeImage,onlySign,indC,cutoffDir) {
                         x <- unlist(res[comp,]);
                         names(x) <- colnames(res)
                         newx <- x
                         xsign <- which(eval(parse(text=paste("abs(x)",cutoffDir,cutoff,sep=""))) & !is.na(x))
                         
                         newx[xsign] <- paste("<a class='sign' href='",basename(pathplot),"/",indC[comp],"_",names(x[xsign]),".",typeImage,"'>",x[xsign]," </a>",sep="");
                         
                         if (!onlySign) {
                             xnosign <- which(!is.na(x))
                             xnosign <- intersect(xnosign,setdiff(1:length(x),xsign))
                             newx[xnosign] <- paste("<a class='normal' href='",basename(pathplot),"/",indC[comp],"_",names(x[xnosign]),".",typeImage,"'>",x[xnosign]," </a>",sep="")
                         }

               
                         return(newx)
                     },res=res,pathplot=pathplot,cutoff=cutoff,typeImage=typeImage, onlySign=onlySign, indC = indComp(icaSet), cutoffDir=cutoffDir)
    res <- if (length(keepVar)==1) as.data.frame(t(as.data.frame(res))) else (as.data.frame(res))
    colnames(res) <- comp
    rownames(res) <- saverownames

    if (onlySign) {
            whichRowSign <- unlist(apply(resnum,1,function(x,cutoff) {x <- as.numeric(x); length(x[eval(parse(text=paste("abs(x)",cutoffDir,cutoff,sep=""))) & !is.na(x)])}, cutoff = cutoff))
            whichCompSign <- unlist(apply(resnum,2,function(x,cutoff) {x <- as.numeric(x); length(x[eval(parse(text=paste("abs(x)",cutoffDir,cutoff,sep=""))) & !is.na(x)])}, cutoff = cutoff))
        whichRowSign <- which(whichRowSign > 0)
        whichCompSign <- which(whichCompSign > 0)
        compSign <- indComp(icaSet)[whichCompSign]

    }
    else {
        whichCompSign <- 1:ncol(res)
        whichRowSign <- 1:nrow(res)
        compSign <- indComp(icaSet)
    }
    
    if (!missing(res2)) {
        rn <- rownames(res)
        cn <- colnames(res)
        res <- as.data.frame(matrix(paste(nameres,"=",as.matrix(res),"\n",nameres2,"=",t(res2),sep=""),ncol=ncol(res)), stringsAsFactors=FALSE, check.names=FALSE)
        rownames(res) <- rn
        colnames(res) <- cn
    }

    ###list densities graphs by variable column and write the references in an html file for each variable
    graphByAnnot <- llply(keepVar[whichRowSign],#colnames(annot)
                          function(colAnnot,path,pathplot,allimages,typeImage) {
                              annotImages <- allimages[grep(x=allimages, pattern =paste("_",colAnnot,".",typeImage,sep=""),fixed = TRUE)]
                              #extract component index and sort by index
                              names(annotImages) <- sapply(basename(annotImages),function(x) strsplit(x,split="_")[[1]][[1]])
                              annotImages <- annotImages[as.character(sort(as.numeric(names(annotImages))))]
                              #write(paste("<IMG SRC='",basename(pathplot),"/",basename(annotImages),"' ALT='Error when loading image' TITLE=''>",sep="\n"), file = paste(path,colAnnot,".htm",sep=""))
                              file <- paste(path,basename(pathplot),"/",colAnnot,".htm",sep="")
                              write(paste("<IMG SRC='",basename(annotImages),"' ALT='Error when loading image' TITLE=''>",sep="\n"), file = file)
                              return(file)

                          },
                          path=path,pathplot=pathplot,allimages=allimages,typeImage=typeImage)
                                    
    graphByAnnot <- sapply(unlist(graphByAnnot),function(x, pathplot) paste(basename(pathplot),basename(x),sep="/"), pathplot = pathplot)
    names(graphByAnnot) <- keepVar[whichRowSign]

    ###list densities graphs by Comp and write the references in an html file for each component
    graphByComp <- llply(compSign,
                          function(indComp,path,pathplot,allimages) {
                              compImages <- allimages[grep(x=basename(allimages), pattern =paste("^",indComp,"_",sep=""),fixed=FALSE)]
                              #write(paste("<IMG SRC='",basename(pathplot),"/",basename(compImages),"' ALT='Error when loading image' TITLE=''>",sep="\n"),file = paste(path,indComp,".htm",sep=""))
                              file <- paste(path,basename(pathplot),"/",indComp,".htm",sep="")
                              write(paste("<IMG SRC='",basename(compImages),"' ALT='Error when loading image' TITLE=''>",sep="\n"),file = file)
                              return(file)
                              
                          },path=path,pathplot=pathplot,allimages=allimages)
    
    graphByComp <- sapply(unlist(graphByComp),basename)
    
    ## link to the genes projections on the component
    ## 3 directories to join the main path with results
    colnames(res)[whichCompSign] <- paste("<a  class='comp' href='",relPath, basename(pathComp),"/",colnames(res)[whichCompSign],".htm'>",colnames(res)[whichCompSign],"</a>","\n","<a  class='comp' href='",basename(pathplot), "/",graphByComp,"'>","(graphs)","</a>",sep="")
    # index of no significant components
    whichNoSign <- which(!(1:ncol(res) %in% whichCompSign))
    colnames(res)[whichNoSign] <-  paste("<a  class='comp' href='",relPath, basename(pathComp),"/",colnames(res)[whichNoSign],".htm'>",colnames(res)[whichNoSign],"</a>",sep="")
    
    ### link 
    res$Annotation <-  rownames(res)

    ## In case the Annotation column does not provide only annotation names but also annotation levels (annotation.level),
    ## I assume that there are only one plot by annotation column, so I have to map each row to the corresponding
    ## annotation file plot.
    ## pour chaque annot je sais le nom annot.level correspondant dans res$Annotation
    ## donc pour chaque fichier d'annot (graphByAnnot) je sais les indices correspondants dans res$Annotation
    ## plus qua retourner la situation

    ## if one result by annotation level
    if (length(res$Annotation)>length(keepVar))
        annot2ind <- unlist(sapply(names(graphByAnnot), function(x, y) which(regexpr(x, y) != -1), y = res$Annotation))
    else ## if one result by annotation
        annot2ind <- unlist(sapply(names(graphByAnnot), function(x, y) match(x, y), y = res$Annotation))
        
    ## if not only one plot by annotation, names of the unlist result will be automatically merged with index to avoid name repetition,
    ## here we want to retrieve the annotation column corresponding to each plot so we delete this number
    ## TODO : if number >= 10 ???
    if (length(annot2ind) > nrow(res)) {
         names(annot2ind) <- gsub(names(annot2ind), pattern = "[0-9]$", replacement = "", perl = TRUE)
         names(annot2ind) <- gsub(names(annot2ind), pattern = "[0-9]$", replacement = "", perl = TRUE)
     }

    ## if some significant results
    if (length(annot2ind)>0) {
        
        names(annot2ind) <- graphByAnnot[names(annot2ind)]
        res$Annotation[intersect(annot2ind,whichRowSign)] <-  paste("<a  class='sign' href='",names(annot2ind[annot2ind %in% intersect(annot2ind,whichRowSign)]),"'>",res$Annotation[intersect(annot2ind,whichRowSign)],"</a>","\n", sep = "")

    }
    res <- res[,c("Annotation",colnames(res)[colnames(res) != "Annotation"])]
    
	style <- 
	"
	<head>
	 <style type='text/css'>

	     a.comp:link{text-decoration:none;color:white;}
	     a.comp:visited{text-decoration:none;color:white;}
	     a.comp:hover{text-decoration:none;font-weight:bold;color:black;}

	     a.sign:link{text-decoration:none;color:black;font-weight:bold;}
	     a.sign:visited{text-decoration:none;color:black;font-weight:bold;}
	     a.sign:hover{text-decoration:none;font-weight:bold;color:#B22222;}

	     a.normal:link{text-decoration:none;color:black;}
	     a.normal:visited{text-decoration:none;color:black;}
	     a.normal:hover{text-decoration:none;font-weight:bold;color:black;}

	     TH.geneTable {
	        padding: 2px;
	            border:2px solid;
	            color: white;
	            border-color: #8B1A1A;
	            text-align: center;
	            font-weight: bold;
	            background-color: #B22222
	     }

	 </style>
	</head>"

      x <- xtable(res, caption = caption)
      x <- capture.output(print(x,
                 type = "html",
                 sanitize.text.function = force,
                 include.rownames=FALSE,
                 caption.placement = "top"))#,
  
      x <- gsub(x,pattern="<TABLE",replacement = "<TABLE style='font-family:Helvetica; font-size:12px; border-top:solid thin black; border=2px;' ", ignore.case = TRUE)
      x <- gsub(x,pattern="<CAPTION",replacement="<CAPTION  style='color:black; text-align:left; font-size:14px;' ", ignore.case = TRUE)
      x <- gsub(x,pattern="<TH>",replacement = "<TH  class='geneTable'>", ignore.case = TRUE)
      x <- paste(style,x,sep="")                
      write(x, file = paste(path,filename,".htm",sep=""))


}



##' From a clustering of samples performed according to their contribution to each component, this function computes the chi-squared test of association between each variable level and the cluster, and summarizes the results in an HTML file. 
##'
##' When \code{doPlot=TRUE}, this function writes an HTML file containing the results of the tests as a table of dimension 'variable levels x components' which contains the p-values of the tests. When a p-value is considered as significant according to the threshold \code{cutoff}, it is written in bold and filled with a link pointing to the corresponding barplot displaying the distribution of the clusters across the levels of the variables.
##'
##' One image is created by plot and located into the sub-directory "plots/" of \code{path}. Each image is named by index-of-component_var.png
##'
##' @title Tests association between clusters of samples and variables
##' @param icaSet An object of class \code{\link[MineICA:IcaSet-class]{IcaSet}}
##' @param params An object of class \code{\link[MineICA:MineICAParams-class]{MineICAParams}} providing the parameters of the analysis
##' @param resClus A list of numeric vectors indexed by sample IDs, which specifies the sample clusters. There must be one clustering by component of \code{icaSet}. The names of the list must correspond to the component indices.
##' @param keepVar The variable labels to be considered, i.e a subset of the variables of icaSet available in \code{varLabels(icaSet)}.
##' @param keepComp A subset of components available in \code{indComp(icaSet)} to be considered, if missing all components are used.
##' @param funClus The name of the function used to perform the clustering (just for text in written files).
##' @param adjustBy The way the p-values of the Wilcoxon and Kruskal-Wallis tests should be corrected for multiple testing: \code{"none"} if no p-value correction has to be done, \code{"component"} if the p-values have to be corrected by component, \code{"variable"} if the p-values have to be corrected by variable.
##' @param testBy Chi-square tests of association can be performed either by \code{"variable"} (one test by variable, default) or by variable \code{"level"} (as many tests as there are annotation levels).
##' @param method The correction method, see \code{\link{p.adjust}} for details, default if \code{"BH"} for Benjamini & Hochberg.
##' @param doPlot If TRUE, the barplots showing the distribution of the annotation levels among the clusters are plotted and the results are provided in an HTML file 'cluster2annot.htm', else no plot is created.
##' @param cutoff The threshold for statistical significance.
##' @param filename File name for test results, if \code{doPlot=TRUE} will be an HTML file else will be a 'txt' file. If missing when \code{doPlot=TRUE},  will be "clusVar".
##' @param path A directory _within resPath(params)_ where the outputs are saved if \code{doPlot=TRUE}, default is \code{'cluster2annot/'}.
##' @param onlySign If TRUE (default), only the significant results are plotted.
##' @param typeImage The type of image file where each plot is saved.
##' @return This function returns a list whose each element gives, for each component,  the results of the association chi-squared tests between the clusters and the annotation levels.  
##' @author Anne Biton
##' @export
##' @seealso \code{clusterSamplesByComp}
##' @examples 
##' ## load an example of IcaSet
##' data(icaSetCarbayo)
##' ## build object of class MineICAParams 
##' params <- buildMineICAParams(resPath="carbayo/")
##'
##' ## cluster samples according to the columns of the mixing matrix A with kmeans in 2 groups
##' resClus <- clusterSamplesByComp(icaSet=icaSetCarbayo, params=params, funClus="kmeans",
##'                                clusterOn="A", nbClus=2)$clus
##'
##' ## specify directory for the function outputs (here same directory as the default one)
##' ## this directory will be created by the function in resPath(params)
##' dir <- "clus2var/"
##'
##' ## compute chi-square tests of association, p-value are not adjusted (adjustBy="none"),
##' # test results are written in txt format (doPlot=FALSE and filename not missing)
##' resChi <- clusVarAnalysis(icaSet=icaSetCarbayo, params=params, resClus=resClus, funClus="kmeans",
##'                           adjustBy="none", doPlot=FALSE, path=dir, filename="clusVarTests")
##'
##' \dontrun{
##' ## compute chi-square tests of association, p-value are not adjusted (adjustBy="none"),
##' # write results and plots in HTML files (doPlot=TRUE)
##' resChi <- clusVarAnalysis(icaSet=icaSetCarbayo, params=params, resClus=resClus, funClus="kmeans",
##'                           path=dir, adjustBy="none", doPlot=TRUE, filename="clusVarTests")
##'
##' ## compute chi-square tests of association by only considering a subset of components and variables,
##' # adjust p-values by component  (adjustBy="component"),
##' # do not write results (doPlot=FALSE and filename is missing).
##' resChi <- clusVarAnalysis(icaSet=icaSetCarbayo, params=params, resClus=resClus, keepComp = 1:10,
##'                           keepVar=c("GENDER","STAGE"), funClus="kmeans", adjustBy="none", 
##'                           doPlot=FALSE)
##' }
##' 
clusVarAnalysis <- function (icaSet, params, resClus, keepVar, keepComp, funClus = "", adjustBy = c("none", "component", "variable"), method = "BH", doPlot = FALSE, cutoff = params["pvalCutoff"], path = paste(resPath(params),"clus2var/",sep=""), onlySign = TRUE, typeImage = "png", testBy=c("variable","level"), filename) {

    
    group <- G <- k <- centers <- uncert <- NULL
    testBy <- match.arg(testBy)
    adjustBy <- match.arg(adjustBy)            
    
    nbClus <- unlist(lapply(resClus,function(x) length(unique(x))))
    if (length(which(nbClus==1))>0) {
        if (missing(keepComp))
            keepComp <- indComp(icaSet)[-which(nbClus==1)]
        else
             keepComp <- c(keepComp,setdiff(keepComp,indComp(icaSet)[which(nbClus==1)]))
    }
    
    if (!missing(keepComp)) {
        icaSet <- icaSet[,,keepComp=keepComp]
        resClus <- resClus[as.character(keepComp)]
    }

    typeImage <- gsub(typeImage,pattern=".",replacement="",fixed=TRUE)
        
    if (missing(keepVar))
        keepVar <- varLabels(icaSet)
    else {
        keepVar <- intersect(varLabels(icaSet),keepVar)
        dif <- setdiff(keepVar,varLabels(icaSet))
        if (length(dif)>0)
            warning(paste("Element",dif,"of arg 'keepVar' is not included in 'varLabels(icaSet)'"))
        if (length(keepVar)==0)
            stop('keepVar is not included in varLabels(icaSet)')
    }

    switch(testBy,
           variable={
               resByComp <- 
                   lapply(as.list(1:ncol(A(icaSet))),
                          function (indC, icaSet, resClus) {
                              
                              clus <- resClus[[indC]][sampleNames(icaSet)]
                              resByAnnot <-
                                  lapply(keepVar, function(keepVar, icaSet, clus) {
                                      
                                      ann <- icaSet[[keepVar]]
                                      tt <- chisq.test(table(ann, clus))
                                      pval <- signif(tt$p.value,4)
                                      
                                  }, icaSet = icaSet, clus = clus)
                              resByAnnot <- unlist(resByAnnot)
                              names(resByAnnot) <- keepVar
                              return(resByAnnot)
                          }
                          , icaSet = icaSet
                          , resClus = resClus
                          )
               resByComp <- data.frame(resByComp, stringsAsFactors = FALSE, check.names = FALSE)               
           },
           level={
               resByComp <- 
                   lapply(as.list(1:ncol(A(icaSet))),
                          function (indC, icaSet, resClus) {
                              clus <- resClus[[indC]]
                              resByAnnot <-
                                  lapply(keepVar, function(keepVar, icaSet, clus) {
                                      
                                      samplesByGroup <- lapply(split(pData(icaSet),pData(icaSet)[[keepVar]]), rownames)
                                      lapply(samplesByGroup,
                                             function(samples, icaSet, clus) {
                                                 
                                                 df <- data.frame(clus = clus)
                                                 df$group <- rep(0,nrow(df))
                                                 rownames(df) <- sampleNames(icaSet)
                                                 df[samples,"group"] <- 1
                                                 
                                                 tt <- chisq.test(table(df$group, df$clus))
                                                 pval <- signif(tt$p.value,4)
                                             }
                                             , icaSet = icaSet
                                             , clus = clus
                                             )
                                      
                                  }, icaSet = icaSet, clus = clus)
                              names(resByAnnot) <- keepVar
                              return(resByAnnot)
                          }
                          , icaSet = icaSet
                          , resClus = resClus
                          )
               
               resByComp <- as.data.frame(lapply(lapply(resByComp,function(x) lapply(x, unlist)),unlist), stringsAsFactors = FALSE, check.names = FALSE)
               
           }
           )

    
    colnames(resByComp) <- indComp(icaSet) 


    if (adjustBy != "none") {

        ## adjust by column
        if (adjustBy == "component") {
            resByComp <- as.data.frame(apply(resByComp, 2, p.adjust, method = method), stringsAsFactors = FALSE, check.names = FALSE)
        }
        else if (adjustBy == "variable") {
            resByComp <- as.data.frame(t(apply(resByComp, 1, p.adjust, method = method)), stringsAsFactors = FALSE, check.names = FALSE)            
        }
    }
    
    if (doPlot) {
        if (missing(path))
            path <- paste(resPath(params),"cluster2annot/",sep="")
        else
            path <- gsub(paste(resPath(params),gsub(path,pattern=resPath(params),replacement=""),"/",sep=""),pattern="//",replacement="/")

        if (missing(filename))
            filename <- "clusVar"
        
        pathplot <- paste(path, "plots/",sep="")
        system(paste("rm -r", pathplot),ignore.stderr=TRUE)
        system(paste("mkdir", path, pathplot),ignore.stderr=TRUE)

        if (length(annot2col(params))==0)
            annot2col(params) <- annot2Color(pData(icaSet))
        
        
        comp2sign <- lapply(as.list(resByComp), function(x) rownames(resByComp)[which(x<=cutoff)])

        signAnnot <- clus <- indC <- NULL
        foreach(signAnnot = comp2sign, clus = resClus, indC = 1:length(comp2sign)) %do% {
            signColAnnot <- c()
            if (length(signAnnot)>0) {
                if (testBy=="variable")
                    colannots <- keepVar
                else
                    colannots <- colnames(pData(icaSet))[which(!is.na(charmatch(signAnnot, x = colnames(pData(icaSet)))))]
                
                if (length(colannots)>0) {
                    for(colannot in colannots) { 
                        if (!(colannot %in% signColAnnot)) {
                            signColAnnot <- c(signColAnnot, colannot)
                            df <- data.frame(cluster = as.factor(clus[rownames(pData(icaSet))]))
                            df$group <- as.factor(pData(icaSet)[[colannot]])
                            rownames(df) <- sampleNames(icaSet)
                            q <- ggplot(df)
                            q <- q + geom_bar(aes(x=cluster,fill = group), position="dodge") + xlab("clusters") + ggtitle(paste("Component ",indComp(icaSet)[indC], ": " , colannot,sep="")) 
                            
                            if (nchar(colannot)>20)
                                colAnnotlegend <- substr(colannot, start = 1, stop = 20)
                            else
                                colAnnotlegend <- colannot
                            
                            q <- q + scale_fill_manual(values = c(annot2col(params)[levels(df$group)]),
                                                       name = colAnnotlegend,
                                                       labels = paste(levels(df$group), " (",table(df$group)[levels(df$group)],")",sep=""),
                                                       breaks = levels(df$group))

                            
                            filen <- paste(pathplot,indComp(icaSet)[indC],"_",colannot,".",typeImage,sep="")
                            ggsave(plot=q,filename=filen,height = 5, width = 6)

                            if (testBy=="level") {
                                filen <- paste(indComp(icaSet)[indC],"_",colannot,".",typeImage,sep="") 
                                levelfiles <- (paste(pathplot,indComp(icaSet)[indC],"_", paste(colannot, unique(pData(icaSet)[[colannot]]), typeImage, sep = "."),sep=""))
                                ## symbolic link to link the levels annotation (rows of resbycomp) to the plot made for the whole annotation column
                                sapply(levelfiles, function(x) system(paste("ln -s", filen, x),ignore.stderr=TRUE))
                            }
                        }
                    }
                }
            }
        }
        
        caption <- paste("<b>Given the clustering obtained with ",funClus, " the association of each ", if (testBy=="variable") "variable" else "variable level" , " with the clusters hosas been tested using chi-squared tests of independence.  This matrix provides the ", if (adjustBy != "none") paste("adjusted p-values (by ", adjustBy, ", with method ", method, ") of these tests.",sep = "") else "p-values of these tests. </b> <br> Click on the significant values to visualize corresponding barplots describing the distribution of the clusters within an annotation group.", sep="")
        writeHtmlResTestsByAnnot (res = signif(t(resByComp), 5), icaSet = icaSet, params = params, path = path, pathplot = "plots/", filename = gsub(filename,pattern=path,replacement=""), typeImage = typeImage, caption = caption, onlySign = onlySign)
    }
    else {
        if (!missing(filename)) {
            if (missing(path))
                path <- paste(resPath(params),"cluster2annot/",sep="")
            else
                path <- gsub(paste(resPath(params),gsub(path,pattern=resPath(params),replacement=""),"/",sep=""),pattern="//",replacement="/")
       
            system(paste("mkdir", path),ignore.stderr=TRUE)

            resByComp <- signif(resByComp, 5)
            write.table(resByComp, file=paste(path,gsub(filename,pattern=path,replacement=""),sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
        }
    }


       return(resByComp)
 

}


