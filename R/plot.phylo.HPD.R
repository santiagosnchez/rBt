#' plot.phylo.HPD
#'
#' This function plots a \code{phylo} BEAST tree read by
#' \code{read.annot.beast} with its 95% HPD node height
#' intervals and adds a x-axis scale. 
#'
#' ** Now includes support for PhyloBayes chronograms **
#'
#' @param x annotated \code{phylo} object
#' @param nodes a vector with node numbers to 
#'              plot. Default is \code{NULL}
#' @param bar.width fraction of 1 for bar thickness
#' @param bar.col bar color
#' @param vline add a vertical lines (default: TRUE)
#' @param pb TRUE if the tree is a chronogram from PhyloBayes with HPD
#' @param border if bar borders should be drawn. See ?rect.
#' @param ... further arguments passed by \code{phylo.plot}
#' @return plot
#' @seealso \code{\link{plot.phylo}} \code{\link{rect}}
#' @importFrom scales alpha
#' @export
#' @examples
#' file <- system.file("data/mcc.tre", package="rBt")
#' tr <- read.annot.beast(file)
#' # plot HPD heights on all nodes:
#' plot.phylo.HPD(tr, cex=0.5, bar.width=0.4, bar.col="red", border=NA)
#' # with library(scales) for transparency
#' # plot.phylo.HPD(tr, cex=0.5, bar.width=0.4, bar.col=alpha("red",0.5), border=NA)
#' # plot HPD heights on specific nodes:
#' tips <- tr$tip.label
#' tipsfac <- factor(sapply(tips, function(x) strsplit(x, "__")[[1]][1] ))
#' names(tipsfac) <- NULL
#' tipslist <- split(tips, tipsfac)
#' allanc <- c(unlist(sapply(tipslist, function(x, y=tr) Ancestors(tr, which.node(y, x)) )), 
#'             sapply(tipslist, function(x) which.node(tr, x)))
#' allanc <- sort(allanc)
#' allanc <- allanc[ !duplicated(allanc) ]
#' names(allanc) <- NULL
#' plot.phylo.HPD(tr, cex=0.5, bar.width=0.4, bar.col="red", border=NA, nodes=allanc)

plot.phylo.HPD <- function(x, nodes=NULL, pb=FALSE, bar.width=0.3, bar.col=NA, border=NULL, 
			   at = NULL, minor=NULL, vline=TRUE, ...){
	op <- par(no.readonly = TRUE)
	plot(x, plot=F, ...)
	ppenv <- get("last_plot.phylo",envir=.PlotPhyloEnv)
	N <- length(x$tip.label)
	yycrds <- ppenv$yy[(N+1):length(ppenv$yy)]
	yycrdsu <- yycrds+bar.width
	yycrdsl <- yycrds-bar.width
	maxxlim <- max(ppenv$x.lim)
	xxmax <- max(ppenv$xx)
	if (pb){
		hpd <- sub("_",",",x$node.label)
		hpd <- matrix(as.numeric(unlist(strsplit(hpd,","))), ncol=2, byrow=T)
		hpd <- rbind(matrix(0,ncol=2, nrow=length(x$tip.label)), hpd)
		hpd <- hpd[,c(2,1)]
	} else {
		if (is.null(x$metadata))
			stop("phylo object was not read by read.annot.beast; try with pb=T, perhaps?")
		hpd <- x$metadata[,"height_95%_HPD"]
		hpd <- matrix(as.numeric(unlist(strsplit(hpd,","))), ncol=2, byrow=T)
	}
	hpdl <- hpd[(N+1):dim(hpd)[1],1]
	hpdu <- hpd[(N+1):dim(hpd)[1],2]
	xxu <- -(hpdu - xxmax)
	xxl <- -(hpdl - xxmax)
	par(op)
	plot(x, x.lim=c(min(xxl, na.rm=TRUE),maxxlim), ...)
	if (!is.null(nodes)){
		rect(xxu[nodes-N], yycrdsl[nodes-N], xxl[nodes-N], yycrdsu[nodes-N], border=border, col=bar.col)
	} else {
		rect(xxu, yycrdsl, xxl, yycrdsu, border=border, col=bar.col)
	}
	ax <- simpleAxisPhylo(at=at, minor=minor)
	if (vline)
		abline(v=ax, lty=2, lwd=0.5)
}
