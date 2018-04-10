#' mcc2
#'
#' Adds a circular time slice to a ploted \code{phylo} tree
#' of \code{type = "fan"}
#'
#' @param phy \code{phylo} tre object
#' @param start start of the time slice (lower)
#' @param end end of the time slice (upper)
#' @param col choose color (default = "grey")
#' @seealso \code{\link{plot.phylo}} \code{\link{sp}}
#' @export circular.time.slice
#' @examples
#' path <- system.file("data/mcc.tre", package="rBt")
#' tr <- read.beast.annot(path)
#' timerange <- range(branching.times(tr))
#' # equally spaced time slices
#' timeslices <- seq(timerange[1], timerange[2], length.out=5)
#' ts <- cbind(timeslices[-5], timeslices[-5]+((timeslices[3]-timeslices[2])/2))
#' # plot tree
#' plot(tr, type="fan", show.tip.label=F)
#' # add slice
#' for (i in 1:dim(ts)[1])
#'    circular.time.slice(tr, ts[i,1], ts[i,2])

circular.time.slice <- function(phy,start,end,col="grey"){
    # see http://r.789695.n4.nabble.com/how-to-draw-a-circle-td798480.html
    t <- seq(0,2*pi,length=100)
    root <- rev(sort(branching.times(phy)))[1]
    st <- root-end
    en <- st+(end-start)
    cs <- t(rbind( 0+sin(t)*st, 0+cos(t)*st)) 
    ce <- t(rbind( 0+sin(t)*en, 0+cos(t)*en)) 
    # http://gis.stackexchange.com/questions/64537/clip-polygon-and-retain-data
    ps <- Polygon(cs, hole=T)
    pe <- Polygon(ce)
    ps <- Polygons(list(ps), ID = "ps")
    pe <- Polygons(list(pe), ID = "pe")
    both <- SpatialPolygons(list(pe, ps))
    spdf = SpatialPolygonsDataFrame(both, data.frame(variable1 = c(NA,NA),
    variable2 = c(NA, NA), row.names = c("pe", "ps")))
    AddHoleToPolygon <-function(poly,hole){
        # http://stackoverflow.com/questions/29624895/how-to-add-a-hole-to-a-polygon-within-a-spatialpolygonsdataframe
        # invert the coordinates for Polygons to flag it as a hole
        coordsHole <-  hole@polygons[[1]]@Polygons[[1]]@coords
        newHole <- Polygon(coordsHole,hole=TRUE)
        # punch the hole in the main poly
        listPol <- poly@polygons[[1]]@Polygons
        listPol[[length(listPol)+1]] <- newHole
        punch <- Polygons(listPol,poly@polygons[[1]]@ID)
        # make the polygon a SpatialPolygonsDataFrame as the entry
        new <- SpatialPolygons(list(punch),proj4string=poly@proj4string)
        new <- SpatialPolygonsDataFrame(new,data=as(poly,"data.frame"))
        return(new)
    }
    ring <-AddHoleToPolygon(spdf[1,],spdf[2,])
    plot(ring, add=T, col=adjustcolor(col,alpha.f=0.5), border=F)
}

