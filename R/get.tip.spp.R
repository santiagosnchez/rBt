#' get.tip.spp
#'
#' Returns a vector is species names based on a \code{phylo}  
#' object's tip names.
#'
#' @param phy a \code{phylo} object
#' @param del a \code{character} string with the delimitator.  
#'             Default is \code{"_"}.               
#' @return a \code{character} vector
#' @export
#' @examples
#' file <- system.file("data/mcc.tre", package="rBt")
#' tr <- read.nexus(file)
#' spp <- get.tip.spp(tr, del="__")
#' tipcols <- rainbow(8)[ factor(spp) ]

get.tip.spp <- function(phy,del="_"){
    if (length(grep(paste(del), phy$tip.label)) == 0)
        stop(paste(del,"was not found in tip.label"))
    else if (length(grep(paste(del), phy$tip.label)) != length(phy$tip.label))
        stop(paste(del,"was not found in all tip.label"))
    re <- unlist(lapply(strsplit(phy$tip.label, paste(del)), function(x) x[1] ))
    return(re)
}
