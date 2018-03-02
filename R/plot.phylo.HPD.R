#' plot.phylo.HPD
#'
#' This function plots a \code{phylo} BEAST tree read by
#' \code{read.beast.annot} with its 95% HPD node height
#' intervals and adds a x-axis scale. 
#'
#' @param x annotated \code{phylo} object
#' @param nodes a vector with node numbers to 
#'              plot. Default is \code{NULL}
#' @param bar.width fraction of 1 for bar thickness
#' @param bar.col bar color
#' @param border if bar borders should be drawn. See ?rect.
#' @param ... further arguments passed by \code{phylo.plot}
#' @return plot
#' @seealso \code{\link{plot.phylo}} \code{\link{rect}}
#' @export
#' @examples
#' file <- system.file("data/mcc.tre", package="rBt")
#' tr <- read.beast.annot(file)
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

plot.phylo.HPD <- function(x, nodes=NULL, bar.width=0.3, bar.col=NA, border=NULL, ...){
	plot(x, plot=F, ...)
	ppenv <- get("last_plot.phylo",envir=.PlotPhyloEnv)
	N <- length(x$tip.label)
	yycrds <- ppenv$yy[(N+1):length(ppenv$yy)]
	yycrdsu <- yycrds+bar.width
	yycrdsl <- yycrds-bar.width
	maxxlim <- max(ppenv$x.lim)
	xxmax <- max(ppenv$xx)
	hpdu <- x$`height_95%_HPD_UPPER`
	hpdl <- x$`height_95%_HPD_LOWER`
	xxu <- -(hpdu - xxmax)
	xxl <- -(hpdl - xxmax)
	plot(x, x.lim=c(min(xxu, na.rm=TRUE),maxxlim), ...)
	if (!is.null(nodes)){
		rect(xxu[nodes-N], yycrdsl[nodes-N], xxl[nodes-N], yycrdsu[nodes-N], border=border, col=bar.col)
	} else {
		rect(xxu, yycrdsl, xxl, yycrdsu, border=border, col=bar.col)
	}
	ax <- axisPhylo()
	abline(v=ax, lty=2, lwd=0.5)
}
