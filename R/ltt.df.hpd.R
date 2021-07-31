#' ltt.df.HPD
#'
#' Produces a list with two data frames; one with
#' mean node values that can be represented as points
#' and or lines, and a another with 95% HPD values that
#' that can be represented as a polygon.
#'
#' It's inteded use is to generate a format that is
#' easily compatible with ggplot2 (see example).
#'
#' @param phy a \code{phylo} object read by \code{read.annot.beast}
#' @return \code{list} object with 2 data frames, one for points and the other for polygons
#' @export
#' @examples
#' file <- system.file("data/mcc.tre", package="rBt")
#' tr <- read.annot.beast(file)
#' # find the root node:
#' my_ltt <- ltt.df.hpd(tr)
#' library(ggplot2)
#' ggplot(my_ltt$points, aes(x=mean, y=N)) +
#'     geom_poly(data=my_ltt$poly, aes(x=x, y=y), fill="grey", alpha=0.5) +
#'     geom_point() + geom_line()

ltt.df.hpd <- function(phy){
    if (is.null(phy$metadata))
        stop("phy was not read by read.annot.beast")
    root = length(phy$tip.label)+1
    mean_height = sort(as.numeric(phy$metadata$height[root:dim(phy$metadata)[1]]), decreasing=T)
    hpd_height = phy$metadata$`height_95%_HPD`[root:dim(phy$metadata)[1]]
    hpd_height = matrix(as.numeric(unlist(strsplit(hpd_height,","))), ncol=2, byrow=T)
    hpd_height_min = sort(hpd_height[,1], decreasing=T)
    hpd_height_max = sort(hpd_height[,2], decreasing=T)
    node_height = -data.frame(mean=mean_height, lower=hpd_height_min, upper=hpd_height_max)
    node_height$N = 1:dim(node_height)[1]
    res = list()
    res$points = node_height[,c("mean","N")]
    res$poly = data.frame(x=c(node_height$lower,rev(node_height$upper)),y=c(node_height$N,rev(node_height$N)))
    return(res)
}
