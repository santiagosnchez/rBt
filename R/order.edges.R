#' order.edges
#'
#' Gets the edge index ordered by node (e.g.
#' order by node index). This function is practical
#' for plotting node values on edges instead of nodes.
#'
#' @param phy tree, \code{phylo} object
#' @return \code{numeric} vector with edge indices
#' @seealso \code{\link{edgelabels}} \code{\link{read.beast.annot}}
#' @note tree files read by \code{read.beast.annot} will include an
#'       element with ordered edges: phy$edge.ordered
#' @export
#' @examples
#' file <- system.file("data/trees/mcc.tre", package="rBt")
#' tr <- read.nexus(tr)
#' nodes <- (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)
#' edges.ordered <- order.edges(tr)
#' plot(tr)
#' edgelabels(edge=edges.ordered, text=nodes[-1]) # -1 is needed to exclude the root node

order.edges <- function(phy){
    pp <- prop.part(phy)
    edg <- sapply(pp, function(x,y=phy) min(which.edge(y,x))-1 )
    edg[1] <- NA
    return(edg)
}
