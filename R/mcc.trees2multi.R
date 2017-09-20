#' mcc.trees2multi
#'
#' This function reads separate MCC files from BEAST
#' an stores them as a \code{multiPhylo} object.
#'
#' @param path The path to the directory with the MCC files.
#'             By default it will list files in the current path(.).
#' @return list of trees of class \code{multiPhylo}
#' @seealso \code{\link{read.beast.annot}} \code{\link{read.nexus}}
#' @export
#' @examples
#' path <- system.file("data/trees/", package="rBt")
#' trs <- mcc.trees2multi(path)
#' class(trs)
#' # [1] "multiPhylo"

mcc.trees2multi <- function(path="."){
    tnames <- dir(path)
    files <- paste0(path,"/",tnames)
    trs <- lapply(files, read.beast.annot)
    msg <- paste("Read",length(trs),"trees")
    message(msg)
    class(trs) <- "multiPhylo"
    names(trs) <- tnames
    return(trs)
}


