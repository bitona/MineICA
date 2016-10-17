heatmap.plus = function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
    distfun = dist, metric = "euclidean", distres = NULL, hclustfun = agnes, method = "ward", reorderfun = function(d, 
        w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
        "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
    margins = c(5, 5), ColSideColors = NULL, RowSideColors = NULL, cexRow = 0.2 + 
        1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, heatmapCol = NULL,
    verbose = getOption("verbose"), legend.args = NULL,breaks=NULL, ...) 
{

    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x, method = metric), method = method)
            ddr <- as.dendrogram(as.hclust(hcr))
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            distRes <- if (is.null(distres)) distfun(if (symm) x else t(x), method = metric) else distres
            hcc <- hclustfun(distRes, method = method)
            ddc <- as.dendrogram(as.hclust(hcc), method = method)

            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 2.5 else 0.05) + if (!is.null(main)) 0.4 else 0, 4)
    #lhei <- c((if (doCdend) 1 else ) + if (!is.null(main)) 0.2 else 0, 4)

    if (!is.null(ColSideColors)) {
        if (!is.matrix(ColSideColors)) 
            stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || dim(ColSideColors)[1] != 
            nc) 
            stop("'ColSideColors' dim()[2] must be of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
	lhei <- c(lhei[1], if (doCdend) 1.4 else 1.3, lhei[2])   #lhei <- c(lhei[1], if (doCdend) 1.4 else 1, lhei[2])
	#lhei <- c(lhei[1], if (doCdend) 0.6 else 0.4, lhei[2])
	
    }
    if (!is.null(RowSideColors)) {
        if (!is.matrix(RowSideColors)) 
            stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || dim(RowSideColors)[1] != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, 
            "; lmat=\n")
        #print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!is.null(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = RowSideColors[rowInd, ]
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
            rsc.colors[rsc.i] = rsc.name
            rsc[rsc == rsc.name] = rsc.i
            rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (ncol(RowSideColors) > 0) {
            axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), 
                las = 2, tick = FALSE)
        }
    }
    if (!is.null(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, ] 
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
            csc.colors[csc.i] = csc.name
            csc[csc == csc.name] = csc.i
            csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (ncol(ColSideColors) > 0) {
            axis(4, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(ColSideColors), 
                las = 2, tick = FALSE, cex.axis = cexCol)
        }
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
        x <- t(x)
    }
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr

    if (is.null(breaks)) {
        ## add color breaks, extend the number of colorsd between the 25th and 75th quantiles
        rm <- range(x)
        qq1 <- quantile(unlist(x),0.25)
        if (qq1 == rm[1])
            qq1 <- quantile(unlist(x),0.35)
        
        qq2 <- quantile(unlist(x),0.75)
        nbCol <- length(heatmapCol)
        extrNbCol <- floor(nbCol/4)-1
        midNbCol <- floor(nbCol/2) 
        breaks <- c(seq(from = rm[1], to = qq1, by = abs(qq1-rm[1])/extrNbCol),
                    seq(from = qq1, to = qq2, by = (qq2-qq1)/midNbCol),
                    seq(from = qq2, to = rm[2], by = (rm[2]-qq2)/extrNbCol))

        image(x = 1:nc, y = 1:nr, z = x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
              c(0, nr), axes = FALSE, xlab = "", ylab = "",
              col = if (!is.null(heatmapCol)) heatmapCol else heat.colors(12),
              breaks = breaks, ...)
    }
    else {
            image(x = 1:nc, y = 1:nr, z = x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
              c(0, nr), axes = FALSE, xlab = "", ylab = "",
              col = if (!is.null(heatmapCol)) heatmapCol else heat.colors(12),
              breaks = breaks, ...) ## add colors
    }
    
    
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)

    #to have legend with bold font
    opbold = par(font = 2)
    if (!is.null(legend.args)) do.call(what = "legend", legend.args)
    par(opbold)

    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend) 
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1.3 else 0, margins[2]))
    if (doCdend)
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
        frame()

    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))

        return(list(x=x,breaks=breaks))

}

