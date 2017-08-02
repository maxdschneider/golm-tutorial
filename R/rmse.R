# Unit Hydrograph - linear storage cascade
# Berry Boessenkool,   berry-b@gmx.de
# Example for package documentation exercise
# root mean square error
# based on wikipedia formula

#' Title
#'
#' @param a observed
#' @param b predicted
#' @param quiet logical variable indicating whether warning messages are printed
#'
#' @return
#' @export
#'
#' @examples
rmse <- function(
                a,
                b,
                quiet=FALSE)
{
if(!(is.vector(a) & is.vector(b))) stop("input is not vectors")
if(length(a) != length(b)) stop("vectors not of equal length")
if(any(is.na(a)|is.na(b)))
   {
   Na <- which(is.na(a)|is.na(b))
   if(!quiet) warning(length(Na), " NAs were omitted from ", length(a), " data points.")
   a <- a[-Na] ; b <- b[-Na]
   } # end if NA
sqrt( sum((a-b)^2)/length(b) )
}
