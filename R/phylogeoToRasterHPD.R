#' gphylogeoToRasterHPD
#'
#' @description
#' Extracts coordinate data from a full posterior BEAST trees file
#' read by \code{read.annot.beast} and generates \code{raster} objects
#' with 2d kernel HPD densities.
#' 
#' @details
#' Thee trees object must be \code{multiPhylo} or \code{list}.
#' It requires two additional arguments: \code{coordnames}, which is 
#' character vector with the column names pointing to the coordinate 
#' data in the BEAST tree (see examples); and \code{intervals}, which is
#' a numeric vector with the time intervals or slices that will be used 
#' to summarize and claculate kernel densities.
#'
#' @param file The path to the MCC file from BEAST
#' @return \code{raster} object
#' @seealso \code{\link{raster}}
#' @export
#' @importFrom MASS kde2d
#' @importFrom scales rescale
#' @examples
#' file <- system.file("data/HBV.trees", package="rBt")
#' trees <- read.annot.beast(file)
#' # remove burnin
#' trees <- trees[1000:10000]
#' # look for column names
#' head(trees[[1]]$metadata)
#' # node            locationsnewTrait
#' #1    1 33.3431887795,134.9307404236
#' #2    2 -2.3045941592,118.5275669247
#' #3    3  15.880672531,106.7863444469
#' #4    4  15.2357758333,97.1170427803
#' #5    5  10.826174473,100.0321841486
#' #6    6 11.8647366524,104.4581350327
#' coordnames <- "locationsnewTrait"
#' # set intervals
#' intervals <- c(0,250,500,750)
#' # get 95% HPD raster
#' st <- phylogeoToRasterHPD(trees, coordnames, intervals)
#' plot(st)
#' # add a map:
#' library(maptools)
#' data(wrld_simpl)
#' par(mfrow=c(2,2))
#' for (i in 1:dim(st)[3]){
#'    image(st[[i]], col=rev(terrain.colors(20)))
#'    plot(wrld_simpl, add=T)
#' }

phylogeoToRasterHPD <- function(trees, coordnames=NULL, intervals=NULL, resolution=0.1){
	if (length(grep("metadata", names(trees[[1]]))) == 0){
		stop("trees need to be read by read.annot.beast")
	}
	if (is.null(coordnames)){
		stop("coordnames needs to be give; try with \"locationsnewTait\" or c(\"lat\",\"long\")")
	}
	coordsNnodes <- matrix(NA,ncol=3)
	colnames(coordsNnodes) <- c("x","y","time")
	cat("Extracting metadata ...\r")
	for (i in 1:length(trees)){
		if (length(coordnames) == 1){
			coords <- matrix(as.numeric(unlist(strsplit(trees[[i]]$metadata[,coordnames],","))), byrow=T, ncol=2)
			coords <- coords[,c(2,1)]
		} 
		else if (length(coordnames) == 2){
			coords <- matrix(c(as.numeric(trees[[i]]$metadata[,grep("lo",coordnames)]),as.numeric(trees[[i]]$metadata[,grep("la",coordnames)])), ncol=2)
		}
		coords <- coords[(length(trees[[i]]$tip.label)+1):dim(coords)[1],]
		coords <- cbind(coords, branching.times(trees[[i]]))
		colnames(coords) <- c("x","y","time")
		coordsNnodes <- rbind(coordsNnodes, coords)
	}
	coordsNnodes <- coordsNnodes[-1,]
	rownames(coordsNnodes) <- 1:dim(coordsNnodes)[1]
	coordsNnodes <- as.data.frame(coordsNnodes)
	coordinates(coordsNnodes) <- ~x+y
	cat("\nChecking intervals ...\r")
	if (is.null(intervals)){
		int <- seq(floor(range(coordsNnodes$time)[1]), ceiling(range(coordsNnodes$time)[2]), length.out=10)
		int <- cbind(int[-10], int[-1])
	} else {
		int <- cbind(intervals[-length(intervals)], intervals[-1])
	}
	all_coords_int <- apply(int, 1, function(x, y=coordsNnodes) y@coords[y$time > x[1] & y$time < x[2], ] )
	names(all_coords_int) <- round(int[,2],0)
	empty <- unlist(lapply(all_coords_int, function(x) dim(x)[1] == 0 ))
	all_coords_int <- all_coords_int[ !empty ]
	cat("\nGenerating 2D kernel density ...\r")
	st <- list()
	for (i in 1:length(all_coords_int)){
		f1 <- kde2d(all_coords_int[[i]][,1], all_coords_int[[i]][,2], n=100)
		f1$z <- rescale(f1$z, c(0,1))
		r <- raster(f1)
		r[ r[] < 0.05] <- NA
		st <- c(st, r)
	}
	cat("\nGetting global extent ...\r")
	d <- matrix(unlist(lapply(st, function(x) as.matrix(extent(x)) )), ncol=2, byrow=T)
	extm <- extent(range(d[,1])[1], range(d[,1])[2], range(d[,2])[1], range(d[,2])[2])
	r <- raster(ext=extm)
	res(r) <- resolution
	cat("\nResampling raster ...\r")
	st <- lapply(st, resample, r)
	st <- stack(st)
	st <- subset(st, c(rev(1:dim(st)[3])))
	names(st) <- rev(names(all_coords_int))
	cat("\nDone\n")
	return(st)
}