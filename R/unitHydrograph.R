#' Unit Hydrograph
#' 
#' Calculate continuous unit hydrograph with given n and k (in the framework of the linear storage cascade)
#' 
#' @return Vector with the unit hydrograph along t
#' @note The sum under the UH should always be 1 (if t is long enough). This needs yet to be checked...
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2013
#' @seealso \code{\link{lsc}} on how to estimate n and k for a given discharge dataset. \code{deconvolution.uh} in the package hydromad, \url{http://hydromad.catchment.org}
#' @keywords hplot ts
#' @export
#' @examples
#' 
#' Time <- 0:100
#' plot(Time, unitHydrograph(n=2,  k=3, t=Time), type="l", las=1,
#'      main="Unit Hydrograph - linear storage cascade")
#' lines(Time, unitHydrograph(n=2,  k=8, t=Time), col=2)
#' lines(Time, unitHydrograph(n=5.5,k=8, t=Time), col=4)
#' text(c(12, 20, 50), c(0.1, 0.04, 0.025), c("n=2, k=3","n=2, k=8","n=5.5, k=8"),
#'      col=c(1,2,4), adj=0)
#' 
#' @param n Numeric. Number of storages in cascade.
#' @param k Numeric. Storage coefficient [1/s] (resistance to let water run out). High damping = slowly reacting landscape = high soil water absorbtion = high k.
#' @param t Numeric, possibly a vector. Time [s].
#' @param force Logical: Force the integral of the hydrograph to be 1? DEFAULT: FALSE
#'
unitHydrograph <- function(
n,
k,
t,
force=FALSE
)
{
if(length(n)>1 | length(k)>1) stop("n and k can only have one single value!
For vectorization, use sapply (see documentation examples).")
UH <- t^(n-1) / k^n / gamma(n) * exp(-t/k)  # some say /k^(n-1) for the second term!
if(force) UH <- UH/sum(UH)
UH
}

