#' gsi
#'
#' Computes genealogical sorting index (Cummings et al. 2008)
#'
#' @param phy        a \code{phylo} object
#' @param groups     a \code{list} where each element is a character
#'                   vector with groups of tips
#' @return a \code{vector} with gsi values
#' @export
#' @references Cummings, M.P., M.C. Neal, K L. Shaw (2008) A genealogical approach to 
#'             quantifying lineage divergence. Evolution 62:2411-2422
#'             \url{doi:10.1111/j.1558-5646.2008.00442.x}
#' @examples
#' set.seed(1)
#' tr = rtree(10)
#' groups = split(tr$tip.label, factor(rep(c(1,2,3), c(3,3,4))))
#' tr_gsi = gsi(tr, groups)
#' tr_gsi
#' #     1     2     3 
#' # 1.000 1.000 0.625

gsi <- function(phy, groups){
    res = vector()
    for (group in groups){
        n = length(group)-1
        anc = get_ancestors(phy, group)
        D = sum(sapply(anc, function(x) length(Children(phy, x))+1-2 ))
        gs = n/D
        allnodes = (length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode)
        I = sum(sapply(allnodes, function(x) length(Children(phy, x))+1-2 ))
        gs.min = n/I
        gs.max = 1
        gsi = (gs-gs.min)/(gs.max-gs.min)
        res = append(res, gsi)
    }
    if (is.null(names(groups))){
        return(res)
    }
    else {
        names(res) = names(groups)
        return(res)
    }
}
