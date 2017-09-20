#' mcc
#'
#' This function is based on the \code{\link{maxCladeCred}}
#' in \code{phangorn}, and it finds the maximum clade 
#' crebility tree given a list of trees. 
#' 
#' 
#'
#' @param phy List of \code{multiPhylo} trees
#' @return mcc tree \code{phylo} object
#' @seealso \code{\link{maxCladeCred}} \code{\link{prop.part}}
#' @importFrom fastmatch fmatch
#' @importFrom phangorn checkLabels
#' @export mcc
#' @examples
#' path <- system.file("data/trees/", package="rBt")
#' trs <- mcc.trees2multi(path)
#' mcctr <- mcc(trs)
#' 
#' 
#' 
#' 
#' 

mcc <- function(phy){
    pp <- prop.part(phy)
    pplabel <- attr(pp, "labels")
    m <- max(attr(pp, "number"))
    nb <- log(attr(pp, "number")/m)
    L <- length(phy)
    res <- numeric(L)
    for (i in 1:L){
        tmp <- phangorn:::checkLabels(phy[[i]], pplabel)
        ppi <- prop.part(tmp)
        indi <- fmatch(ppi, pp)
        if (any(is.na(indi))) 
            res[i] <- -Inf
        else res[i] <- sum(nb[indi])
    }
    k <- which.max(res)
    cat("Clade credibility (log):",res[k],"\n")
    tr <- phy[[k]]
    tr$clade.credibility <- res[k]
    ppk <- prop.part(tr)
    pmt <- matrix(NA,ncol=L, nrow=length(ppk))
    for (i in 1:L){
        ppi <- prop.part(phy[[i]])
        indi <- fmatch(ppk, ppi)
        pmt[,i] <- phy[[i]]$posterior[indi]
    }
    mcp <- rowSums(pmt, na.rm=TRUE)/L
    tr$MCposterior <- mcp
    return(tr)
}
