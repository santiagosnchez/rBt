#' RplotEBS
#'
#' Plots one or multiple EBSP given a CSV ouput from BEAST.
#'
#' @param path       a \code{character} string with the path to 
#'                   CSV files (\code{character}). Default is ".".
#' @param pattern    a \code{character} string pattern to search 
#'                   for CSV files in the \code{path}. Default is ".csv".  
#' @param trim.x     a \code{numeric} value indicating the lower limit
#'                   to the x-axis.
#' @param trim.y     a \code{numeric} value indicating the upper limit
#'                   to the y-axis. 
#' @return a \code{plot}
#' @export
#' @examples
#' file <- system.file("data/EBSP.csv", package="rBt")
#' RplotEBS(file)

RplotEBS <- function(file=NULL, path=".", pattern=".csv", trim.x=NULL, trim.y=NULL, ...){
	if (is.null(file)){
		files <- dir(paste(path), pattern=paste(pattern))
		if (length(files) == 0)
			stop('the path to the csv files is incorrect')
		names <- sub(".csv", "", files)
	} else {
		names <- sub(".*\\/","",file)
		names <- sub(".csv", "", file)
	}
	file.numb <- length(names)
	data = list("data.frame", file.numb)
	min.l = vector()
	max.u = vector()
	min.t = vector()
	for (i in 1:file.numb){
		if (is.null(file)){
			data[[i]] = read.csv(paste(path, "/", names[i], ".csv", sep=""))
		} else {
			data[[i]] = read.csv(file)
		}
		time = -data[[i]]$time
		data[[i]]$time = time
		min.l[i] = min(data[[i]]$hpd.lower.95)
		max.u[i] = max(data[[i]]$hpd.upper.95)
		min.t[i] = min(time)
	}
	for (i in 1:file.numb){
		if (!is.null(trim.x) && !is.null(y.eq)){
			plot(data[[i]]$hpd.lower.95~data[[i]]$time, cex=1, type="n", ylim=c(min(min.l),max(max.u)), xlim=c(trim.x,0), ylab=expression(paste("Population (", theta, ")")), xlab="Time", main=names[i], ...)
		} else if (is.null(trim.x) && is.null(y.eq) && is.null(trim.y)) {
			plot(data[[i]]$hpd.lower.95~data[[i]]$time, cex=1, type="n", xlim=c(min(min.t),0), ylim=c(min(data[[i]]$hpd.lower.95),max(data[[i]]$hpd.upper.95)), ylab=expression(paste("Population (", theta, ")")), xlab="Time", main=names[i], ...)
		} else if (!is.null(trim.x) && is.null(y.eq) && is.null(trim.y)){
			plot(data[[i]]$hpd.lower.95~data[[i]]$time, cex=1, type="n", xlim=c(trim.x,0), ylim=c(min(data[[i]]$hpd.lower.95),max(data[[i]]$hpd.upper.95)), ylab=expression(paste("Population (", theta, ")")), xlab="Time", main=names[i], ...)
		} else if (is.null(trim.x) && !is.null(y.eq) && is.null(trim.y)){
			plot(data[[i]]$hpd.lower.95~data[[i]]$time, cex=1, type="n", ylim=c(min(min.l),max(max.u)), xlim=c(min(min.t),0), ylab=expression(paste("Population (", theta, ")")), xlab="Time", main=names[i], ...)
		} else if (is.null(trim.x) && !is.null(trim.y) && is.null(y.eq)){
			plot(data[[i]]$hpd.lower.95~data[[i]]$time, cex=1, type="n", ylim=c(0,trim.y), xlim=c(min(min.t),0), ylab=expression(paste("Population (", theta, ")")), xlab="Time", main=names[i], ...)
		} else if (!is.null(trim.x) && is.null(trim.y) && is.null(y.eq)){
			plot(data[[i]]$hpd.lower.95~data[[i]]$time, cex=1, type="n", xlim=c(trim.x,0), ylim=c(min(data[[i]]$hpd.lower.95),max(data[[i]]$hpd.upper.95)), ylab=expression(paste("Population (", theta, ")")), xlab="Time", main=names[i], ...)
		} else if (!is.null(trim.x) && !is.null(trim.y) && is.null(y.eq)){
			plot(data[[i]]$hpd.lower.95~data[[i]]$time, cex=1, type="n", xlim=c(trim.x,0), ylim=c(0,trim.y), ylab=expression(paste("Population (", theta, ")")), xlab="Time", main=names[i], ...)
		}
		x = c(data[[i]]$time, rev(data[[i]]$time))
		y = c(data[[i]]$hpd.lower.95, rev(data[[i]]$hpd.upper.95))
		polygon(x = x, y= y, col="grey60", border=NA)
		lines(data[[i]]$mean~data[[i]]$time)
		lines(data[[i]]$median~data[[i]]$time, lty=2)
		axis(4, labels=F)
	}
}
