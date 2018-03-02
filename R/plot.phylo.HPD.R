


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
