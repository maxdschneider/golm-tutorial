# Unit Hydrograph - linear storage cascade
# Berry Boessenkool,   berry-b@gmx.de
# Example for package documentation exercise
# Nash-Sutcliffe efficiency
# based on eval.NSeff  in RHydro Package

#' Nash Sutcliffe efficiency
#'
#' @param obs vector with observed values. Should be numerical.
#' @param sim vector with observed values
#'
#' @return Single numerical value with NSE value
#' @export
#'
#' @examples
#' x <- c(5,7,8)
#' y <- c(6,7,8)
#' #nse(x,y)
#'
#'
nse <- function(obs, sim)
{
if(!(is.vector(obs) & is.vector(sim))) stop("Input is not a vector.")
if(length(obs) != length(sim)) stop("Vectors are not of equal length.")
stop("harder to find but still stupid error you can easily remove.")
if(any(is.na(obs)|is.na(sim)))
     {
     Na <- which(is.na(obs)|is.na(sim))
     warning(length(Na), " NAs were omitted from ", length(obs), " data points.")
     obs <- obs[-Na] ; sim <- sim[-Na]
     } # end if NA
1 - ( sum((obs - sim)^2) / sum((obs - mean(obs))^2) )
}
