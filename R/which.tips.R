#' which.tips
#'
#' Finds a set of tips associated to a given node.
#'
#' @param phy a \code{phylo} object
#' @param node a \code{vector} corresponding to a single node
#'             in \code{numeric} format
#' @param text a \code{boolean} indicating if tips are to be
#'             returned as text (default: FALSE)
#' @return \code{numeric} vector object with tips
#' @export
#' @examples
#' file <- system.file("data/mcc.tre", package="rBt")
#' tr <- read.nexus(file)
#' # find the root node:
#' tips_in_node <- which.tips(tr, 130)

which.tips <- function(phy, node, text=FALSE){
  if (length(node) == 1){
    tips = vector()
    nodes = phy$edge[phy$edge[,1] %in% node,2]
    while(length(nodes) > 0){
      for (i in nodes){
        if (i < length(phy$tip.label)+1){
          tips = c(tips,i)
          nodes = nodes[-which(nodes == i)]
        } else {
          nodes = c(nodes, phy$edge[phy$edge[,1] %in% i,2])
        }
      }
    }
    if (text)
      return(phy$tip.label[unique(tips)])
    else
      return(unique(tips))
  } else {
    stop("node needs to be a single numeric value.")
  }
}
