#!/usr/bin/Rscript
# ©Santiago Sanchez-Ramirez

args <- commandArgs(trailingOnly=TRUE)

if (length(grep("help", args)) != 0){
	stop("\n\nTry:\nRscript RplotEBS.R path=PATH/TO/CSV/FILES pattern=.csv trim.x=-2.0 trim.y=0.6 y.eq=TRUE
Defaults: path=. pattern=csv trim.x=NULL trim.y=NULL y.eq=FALSE\n\n")
}

if (length(grep("path=", args)) != 0){
	path <- strsplit(args[grep("path=", args)], "=")[[1]][2]
} else {
	path = '.'
}

if (length(grep("pattern=", args)) != 0){
	pattern <- strsplit(args[grep("pattern=", args)], "=")[[1]][2]
} else {
	pattern = "csv"
}

if (length(grep("trim.x=", args)) != 0){
	trim.x <- as.numeric(strsplit(args[grep("trim.x=", args)], "=")[[1]][2])
} else {
	trim.x = NULL
}

if (length(grep("trim.y=", args)) != 0){
	trim.y <- as.numeric(strsplit(args[grep("trim.y=", args)], "=")[[1]][2])
} else {
	trim.y = NULL
}

if (length(grep("y.eq=", args)) != 0){
	y.eq <- as.logical(strsplit(args[grep("y.eq=", args)], "=")[[1]][2])
	if (y.eq != TRUE)
		y.eq == NULL
} else {
	y.eq = NULL
}

if (length(grep("panel=", args)) != 0){
	panel <- strsplit(args[grep("panel=", args)], "=")[[1]][2]
	n <- as.numeric(strsplit(panel, ",")[[1]][1])
	t <- as.numeric(strsplit(panel, ",")[[1]][2])
} else {
	panel = NULL
}

file.numb <- length(dir(paste(path), pattern=pattern))

if (file.numb == 0){
	stop('Pattern not found in directory')
}

RplotEBS <- function(path, pattern, file.numb, trim.x, trim.y, ...){
	files = dir(paste(path), pattern=paste(pattern))
	if (length(files) == 0)
		stop('the path to the csv files is incorrect')
	names = sub(".csv", "", files)
	data = list("data.frame", file.numb)
	min.l = vector()
	max.u = vector()
	min.t = vector()
	for (i in 1:file.numb){
		data[[i]] = read.csv(paste(path, "/", names[i], ".csv", sep=""))
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
quartz()

if (is.null(panel)){
	if (file.numb==2){
		#pdf(file=paste(path,'/',"EBPS.pdf", sep=''), height=3, width=3*file.numb)
		par(mfrow=c(1,2))
	} else if (file.numb==3){
		#pdf(file=paste(path,'/',"EBPS.pdf", sep=''), height=7, width=7*file.numb)
		par(mfrow=c(1,3))
	} else if (file.numb==4){
		#pdf(file=paste(path,'/',"EBPS.pdf", sep=''), height=7*(file.numb/2), width=7*(file.numb/2))
		par(mfrow=c(2,2))
	} else if (file.numb>4){
		stop("\n\nWith more than 4 files panel configuration needs to be specified through the argument \"panel\"\ne.g. panel=2,3\n\n")
	}
} else {
	par(mfrow=c(n,t))
}
RplotEBS(path,pattern,file.numb,trim.x,trim.y,las=1)
cat('Press <control><c> when ready to close', "\n")
Sys.sleep(Inf)
dev.off()
