#' simpleAxisPhylo
#'
#' Adds an x-axis to the latest plotted tree with minor
#' with custom tickmarks and minor tickmarks.
#'
#' @param at a vector of points where to place major tickmarks
#' @param minor NULL or numeric factor by which to divide intervals
#' @param ... further arguments passed to \code{axis}
#' @return x-axis (side 1)
#' @seealso \code{\link{axis}} \code{\link{axisPhylo}}
#' @export
#' @examples
#' 
#' # simulate coalescent tree and plot
#' set.seed(1)
#' tr <- rcoal(10)
#' plot(tr)
#'
#' # add default axis
#' simpleAxisPhylo()
#' 
#' # add custom scale
#' plot(tr)
#' trange <- c(0, max(branching.times(tr)))
#' trange <- pretty(trange)
#' x <- seq(trange[1], rev(trange)[1], (diff(trange)/2)[1])
#' simpleAxisPhylo(at=x, minor=2)

simpleAxisPhylo <- function (at = NULL, minor=NULL, ...) 
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    type <- lastPP$type
    if (type %in% c("phylogram", "cladogram")) {
        xscale <- range(lastPP$xx)
        tscale <- c(0, xscale[2] - xscale[1])
        tscale <- tscale[2:1]
        beta <- diff(xscale)/diff(tscale)
        alpha <- xscale[1] - beta * tscale[1]
        if (is.null(at)){
            lab <- pretty(tscale)
            if (is.null(minor)){
                minor.tick <- seq(lab[1],rev(lab)[1],diff(lab)[1]/5)
                minor.tick <- minor.tick[ !minor.tick %in% lab ]
            } else {
                minor.tick <- seq(lab[1],rev(lab)[1],diff(lab)[1]/minor)
                minor.tick <- minor.tick[ !minor.tick %in% lab ]
            }
        }
        else {
            lab <- at
            if (is.null(minor)){
                minor.tick <- seq(lab[1],rev(lab)[1],diff(lab)[1]/5)
                minor.tick <- minor.tick[ !minor.tick %in% lab ]
            } else {
                minor.tick <- seq(lab[1],rev(lab)[1],diff(lab)[1]/minor)
                minor.tick <- minor.tick[ !minor.tick %in% lab ]
            }
        }
        x <- beta * lab + alpha
        xx <- beta * minor.tick + alpha
        axis(side = 1, at = xx, labels = FALSE, lwd = 0, lwd.ticks = 0.5, tck = -.01)
        axis(side = 1, at = x, labels = lab, ...)
    } else {
        stop("type needs to by phylogram or cladogram")
    }
}


