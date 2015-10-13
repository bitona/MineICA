##' This function returns the matrices that will be used to plot the heatmaps of each component. It restricts the data matrix of the \code{\link[MineICA:IcaSet-class]{icaSet}} object to the contributing genes/features, and order the features/genes and samples.
##'
##' This function is called by function \code{\link{plot_heatmapsOnSel}} and is not likely to be called alone.
##' @title Build the heatmap matrices
##' @param icaSet The IcaSet object
##' @param selCutoff The threshold used to select the contributing features/genes based on their projection values. Must be either of length 1 and the same treshold is applied to all components, or of length equal to the number of components and one specific threshold is used for each component.
##' @param selectionByComp The list of gene projections per components already restricted to the contributing genes
##' @param level A character indicating which data level is used to plot the heatmaps: 'features' to plot measured feature levels (e.g probe sets expression values), 'genes' to plot measured gene values (e.g gene expression values).
##' @param samplesOrder A list providing the order of the samples, per component, to be used in the heatmaps. If NULL, the contribution values of the samples are used to rank the columns of the heatmaps.
##' @param featuresOrder A list providing the features or genes order, per component, to be used in the heatmaps. If NULL, the projection values of the genes are used to rank the rows of the heatmaps.
##' @return A list of matrices 
##' @author Anne Biton
##' @keywords internal
build_sortHeatmap <- function (icaSet,
                               selCutoff,
                               selectionByComp,
                               level = c("features","genes"),
                               samplesOrder,
                               featuresOrder
                               ) {

        A <- A(icaSet) 
        level <- match.arg(level)
        switch(level,
               features={S=Slist(icaSet)},
               genes={S=SlistByGene(icaSet)}
               )
        
        if (missing(samplesOrder))
            samplesOrder <- NULL
        
        if (missing(featuresOrder))
            featuresOrder <- NULL
        
 	if (missing(selectionByComp))             
           selectionByComp <- selectContrib(S, cutoff = selCutoff)
       

        if (level == "features") {
            data.expr <- assayDataElement(icaSet,"dat")
        }
        else if (level == "genes") { #if selection was made on genes
            data.expr <- datByGene(icaSet) 
            
        }


        mlist = list()


        for (i in 1:ncol(A)) {

            if (length(selectionByComp[[i]])>0) {

                if (length(selectionByComp[[i]])>1)  
                    m <- data.expr[names(selectionByComp[[i]]), ]
                else {
                    m <- as.data.frame(t(data.expr[names(selectionByComp[[i]]), ]), row.names=names(selectionByComp[[i]]))
                    colnames(m) <- colnames(data.expr)
                }
                    
                if (is.null(samplesOrder)) {
                    proj <- A[,i]
                    names(proj) <- rownames(A)
                    samplesOrderComp <- names(sort(proj, decreasing = TRUE))
                }
                else
                    samplesOrderComp <- samplesOrder[[i]]
         
                if (is.null(featuresOrder)) {
                    maxabs <- which.max(abs(selectionByComp[[i]]))
                    signMax <- sign(selectionByComp[[i]][maxabs])
                    featuresOrder<- names(sort(selectionByComp[[i]], decreasing = if (signMax==1) FALSE else TRUE))
                }
                m <- m[,samplesOrderComp]
                m <- m[featuresOrder,]
                featuresOrder <- NULL

            }
            else
                m <- data.frame()
            
            mlist[[i]] <- m


    }
    names(mlist) = indComp(icaSet)
    return(mlist)


}




##' This function plots the heatmaps representing the measured values of the contributing features/genes on each component. It also plots the sample annotations above each heatmap using colours.
##' 
##' This function restricts the data matrix of an \code{\link{IcaSet}} object to the contributing genes/features, and order features/genes and samples either as asked by the user or according to their values in the ICA decomposition.
##'
##' The heatmap is plotted using a slightly modified version of the function \code{heatmap.plus} from the package of the same name.
##' By default in this function, the hierarchical clustering is calculated using the function \code{\link[cluster]{agnes}} with euclidean metric and Ward's method.
##'  
##' 
##' @title Plot heatmap associated with each component
##' @param icaSet The IcaSet object
##' @param selCutoff A numeric threshold used to select the contributing genes based on their projection values. Must be either of length 1 and the same treshold is applied to all components, or of length equal to the number of components and one specific threshold is used for each component.
##' @param samplesOrder A list providing the order of the samples, per component, to be used in the heatmaps. If missing, the contribution values of the samples are used to rank the columns of the heatmaps.
##' @param featuresOrder A list providing the order of the genes, per component, to be used in the heatmaps. If missing, the projection values of the genes are used to rank the rows of the heatmaps.
##' @param selectionByComp A list of gene projections per component already restricted to the contributing genes, if missing is computed by the function.
##' @param level A character indicating which data level is used to plot the heatmaps: either \code{'features'} to represent the data at the feature levels (e.g expression profiles of probe sets), or \code{'genes'} to represent the data at the annotated-features level (e.g gene expression profiles).
##' @param keepVar The variable labels to be considered, i.e a subset of the column labels of the pheno data of icaSet available in  (\code{varLabels(icaSet)}) 
##' @param keepComp A subset of components, must be included in \code{indComp(icaSet)}. By default, all components are used.
##' @param doSamplesDendro A logical indicating whether a hierarchical clustering has to be performed on the data matrix restricted to the contributing features/genes, and whether the corresponding dendrogram has to be plotted, default is TRUE.
##' @param doGenesDendro A logical indicating if the dendrogram of features/genes has to be plotted, default is FALSE.
##' @param heatmapCol A list of colors used to for heatmap coloring (see argument \code{col} of the function \code{image}).
##' @param file A character to add to each pdf file name. This function creates one file by component named "index-of-component_\code{file}.pdf" .
##' @param path A directory for the output pdf files, must end with "/". Default is current directory.
##' @param annot2col A vector of colours indexed by the levels of the variables of \code{icaSet} (i.e all the annotation values available in \code{pData(icaSet)}). If missing the colours are generated automatically using the function \code{annot2Color}
##' @param ... Additional parameters for function \code{heatmap.plus}
##' @return A list with one element per component, each of them being a list consisting of three elements: \describe{\item{x}{the matrix represented by the heatmap},\item{breaks}{the breaks used for the colours of the heatmap},\item{dendro}{the dendrogram}.}
##' @seealso \code{heatmap.plus}, \code{\link{image}}, \code{\link{annot2Color}}, \code{\link{build_sortHeatmap}}
##' @author Anne Biton
##' @export
##' @examples \dontrun{
##' ## load an example of IcaSet object
##' data(icaSetCarbayo)
##'
##' ## check which variables you would like to use in the heatmap 
##' varLabels(icaSetCarbayo)
##' keepVar <- c("STAGE","SEX")
##' ## Use only component 1
##' keepComp <- 1
##' 
##' ## For each component, select contributing *genes* using a threshold of 2 on the absolute projection values,
##' ## and plot heatmaps of these contributing genes by ordering genes and samples according to their contribution values
##' plot_heatmapsOnSel(icaSet = icaSetCarbayo, selCutoff = 2, level = "genes", keepVar = keepVar,
##'                    keepComp=1, doSamplesDendro = TRUE, doGenesDendro = TRUE,
##'                    heatmapCol = maPalette(low = "blue",high = "red", mid = "yellow", k=44),
##'                    file = "heatmapWithoutDendro_zval3.pdf")
##' 
##' ## For each considered component, select contributing *features* using a threshold of 2 on the absolute projection values,
##' ## and plot heatmaps of these contributing genes with dendrograms
##' plot_heatmapsOnSel(icaSet = icaSetCarbayo, selCutoff = 2, level = "features", keepVar = keepVar,
##'                    keepComp=1, doSamplesDendro = TRUE, doGenesDendro = TRUE,
##'                    heatmapCol = maPalette(low = "blue",high = "red", mid = "yellow", k=44),
##'                    file = "heatmapWithDendro_zval3.pdf")
##'
##' 
##'
###' }
plot_heatmapsOnSel <- function (icaSet,
                                selCutoff = 4,
                                level = c("features","genes"),
                                samplesOrder,
                                featuresOrder,
                                selectionByComp,
                                keepVar,
                                keepComp=indComp(icaSet),
                                doSamplesDendro = TRUE,
                                doGenesDendro = TRUE,
                                heatmapCol=maPalette(low = "blue",high = "red", mid = "yellow", k=44),
                                file = "",
                                path = "",
                                annot2col, ...
                               


) {


    level <- match.arg(level)
    
    icaSet <- icaSet[,, keepComp]

    annot <- pData(icaSet)
    if (!missing(keepVar)) {
        keepVar <- intersect(keepVar,varLabels(icaSet))
        
        if (length(keepVar)==0)
            annot <- NULL
        else {
            if (length(keepVar)==1) {
                annot <- data.frame(annot[,keepVar],row.names=sampleNames(icaSet))
                colnames(annot) <- keepVar
            }
            else
                annot <- annot[,keepVar]
        }
               
        
    }
    else if (ncol(annot)==0)
        annot <- NULL

    if (length(selCutoff)==1)
        selCutoff <- rep(selCutoff[1],nbComp(icaSet))
    else
        if ((length(selCutoff)>1) & (length(selCutoff)<nbComp(icaSet))) {
            cutoff <- rep(selCutoff[1],nbComp(icaSet))
            warning("The length of arg 'selCutoff' must be either 1 or equal to the number of components")    
        }


	### Selection of contributing genes by comp, threshold = zval.cutoff
	### Building of sorted expression matrix based on these genes and with selected samples sorted by values on components
	exprSub_sort <- build_sortHeatmap(icaSet = icaSet,
                                          selCutoff = selCutoff,
                                          selectionByComp=selectionByComp,
                                          level = level,
                                          samplesOrder = samplesOrder,
                                          featuresOrder = featuresOrder
                                          )

        if (length(selCutoff) == 1)
            selCutoff <- rep(selCutoff,nbComp(icaSet))

        if (!is.null(annot)) {

            if (missing(annot2col))
                    annot2col <- annot2Color(annot)

            annot2color <- annot2col
            names(annot2col) <- tolower(names(annot2col))
            mcol <- apply(annot, MARGIN = 2, 
                       function (x,annot2col) {
                           x = tolower(x)
                           x = annot2col[gsub(x = as.character(x), pattern = " ", replacement = "")]
                       }
                       , annot2col = annot2col
                       
                       )
            rownames(mcol) = rownames(annot)
            annot2col_m <- mcol
            annot2col_v <- annot2col
            

        }
        else {
            annot2col.res <-  NULL
            annot2col_m<- NULL
            annot2col_v <- NULL
        }
        
  
    if (!is.null(annot)) {
        colors <- colornames <- borders <- c()
        foreach(varlev=annot, var=colnames(annot)) %do% {
            colors <- c(colors,0,annot2col_v[tolower(unique(varlev))],0)
        }
        
        ncol <- if (ceiling(length(colors)/26)<4) ceiling(length(colors)/26) else 3
        if (length(colors)<ncol*26)
            nbcolbycol <- ceiling(length(colors)/ncol)
        else
            nbcolbycol <- 26
        
        pdf(file = paste(path,dirname(file),"/variableLegend.pdf",sep=""), width = 12, height = 11, paper = "a4r", title=paste(path,dirname(file),"/variableLegend.pdf",sep=""))
        plot.new()

        colors <- colornames <- borders <- c()
        indvar <- trueIndVar <- 0
        foreach(varlev=annot, var=colnames(annot), trueIndVar=1:ncol(annot)) %do% {
            if (indvar>0) {
                if (nbcolbycol-(length(colors)%%nbcolbycol)<(length(unique(varlev))+2))
                    repwhite <- nbcolbycol-(length(colors)%%nbcolbycol)
                else
                    repwhite <- 1
            }
            else
                repwhite <- 0


            if ((length(colors)+length(unique(varlev))+repwhite+1)>ncol*nbcolbycol) {
                legend(x=0, y=1, legend=colornames,fill = colors, border=borders,title=expression(bold("Variable colours")), ncol = ceiling(length(colors)/nbcolbycol), bty = "n", y.intersp=1)
                plot.new()
                colors <- colornames <- borders <- c()
                repwhite <- 0
                indvar <- 0
            }
            else {    
                varlev <- as.character(varlev)
                varlev[is.na(as.character(varlev))] <- "NA"
                colornames <- c(colornames,rep("",repwhite),eval(parse(text=paste("expression(bold('",var,"'))", sep=""))),sort(unique(varlev)))
                colors <- c(colors,rep(0,repwhite+1),annot2col_v[tolower(sort(unique(varlev)))])
                borders <- c(borders,rep("white",repwhite+1),rep("black",length(unique(varlev))))
                indvar <- indvar+1

                if (trueIndVar==ncol(annot)) {
                    legend(x=0, y=1, legend=colornames,fill = colors, border=borders,title=expression(bold("Variable colours")), ncol = ceiling(length(colors)/nbcolbycol), bty = "n", y.intersp=1)
                }
                    
            }
        }

                        
        
        dev.off()
    }

    
    i <- resH <- NULL
    res <- 
         foreach(m=exprSub_sort, i=1:length(exprSub_sort), .inorder=TRUE) %dopar% {
                    

                    # to solve problems with a row with only NA
                    if ("NA" %in% rownames(m))
                        m <- m[-which(rownames(m)=="NA"),]
                    
                    if(ncol(m) != 0) {
                        ff <- paste(path,dirname(file),"/",indComp(icaSet)[i],"_",gsub(x=basename(file),pattern=".pdf",replacement="",fixed=TRUE),"zval",selCutoff[i],".pdf",sep="")
                        pdf(file = ff, title=ff, width = 12, height = 11, paper = "a4r")
                        ## at least two samples are required to plot heatmap
                        if (!is.null(annot)) 
                            annot2col_m <- data.frame(annot2col_m[colnames(m),],stringsAsFactors=FALSE,check.names=FALSE)
                        
                        if (nrow(m) >2) {	

                                Colv <- if (doSamplesDendro) NULL else NA
                                Rowv <- if (doGenesDendro) NULL else NA

                                ## rm <- range(m)
                                ## breaks <- c(seq(from = 0, to = 0.5, by = (0+0.5)/10),
                                ##              seq(from = 0.5, to = 1, by = (1-0.5)/11)
                                ##              )
                                resH <- heatmap.plus(as.matrix(m[,ncol(m):1]),
                                                     Rowv = Rowv,
                                                     Colv = Colv,
                                                     scale = "none", 
                                                     ColSideColors = if (!is.null(annot)) as.matrix(annot2col_m[ncol(m):1,]) else NULL,
                                                     margins = c(5,8), 
                                                     main = paste("Component ",if (!is.null(compNames(icaSet))) compNames(icaSet)[i] else i, ", cutoff = ", selCutoff[i], sep = ""),
                                                     cex.main = 0.8,
                                                     heatmapCol = heatmapCol,
                                                     legend.args = NULL,
                                                     ...
                                                     )
                                breaks <- resH$breaks

                          }
                          else {
                              plot(1, type = "n", axes = "F",xlab = "", ylab = "", main = paste("Component ",compNames(icaSet)[i],": Less than 2 genes have scaled projection higher than ", selCutoff[i], ".",sep = ""))
                              breaks <- NULL
                              resH <- NULL
                          }


                          if (nrow(m) >2 & !is.null(breaks)) {
                                op1 = par(fig = c(0.87,1,0,0.4), new=TRUE)
                                op2 = par(cex.axis = 0.65)
                                ## breaks identical to the ones in heatmap.plus
                                 maColorBar(breaks[which(odd(1:length(breaks)))],
                                           main = "",
                                           horizontal = FALSE,
                                           cex.axis = 0.8,
                                           col = heatmapCol)
                                par(op1)
                                par(op2)
                          }
                        dev.off()
                    }



                    return(resH)
                    
                }


	
    return(res)
}

