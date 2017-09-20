#' which.node
#'
#' Finds the a node given a vector of tip labels or 
#' tip indices. 
#'
#' @param phy a \code{phylo} object
#' @param tips a \code{vector} corresponding to a set tips. 
#'             It can be either type \code{character} or \code{numeric}               
#' @return \code{numeric} object with a node number
#' @export
#' @examples
#' file <- system.file("data/mcc.tre", package="rBt")
#' tr <- read.nexus(file)
#' # find the root node:
#' rnode <- which.node(tr, tr$tip.label)

which.node <- function(phy, tips){
   if (is.character(tips) == TRUE){
       tips_n = vector("numeric",length(tips))
       for (i in 1:length(tips))
           tips_n[i] = which(phy$tip.label == tips[i])
       rows = vector("numeric",length(tips_n))
       for (i in 1:length(tips_n))
           rows[i] = which(phy$edge[,2] == tips_n[i])
       nodes = sort(phy$edge[rows[which.min(rows)]:rows[which.max(rows)],1])
       return(nodes[1])
   } else if (is.numeric(tips) == TRUE){
       rows = vector("numeric",length(tips))
       for (i in 1:length(tips))
           rows[i] = which(phy$edge[,2] == tips[i])
       nodes = sort(phy$edge[rows[which.min(rows)]:rows[which.max(rows)],1])
       return(nodes[1])
   } else if (is.null(tips) == TRUE){
       stop('tips vector is empty')
   }
}
