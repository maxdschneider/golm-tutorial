#' Linear storage cascade, unit hydrograph
#'
#' Optimize the parameters for unit hydrograph as in the framework of the
#' linear storage cascade. Plot observed & simulated data
#'
#' @return \emph{Either} vector with optimized n and k and the Nash-Sutcliffe Index,
#'         \emph{or} simulated discharge, depending on the value of \code{returnsim}
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2013
#' @seealso \code{\link{unitHydrograph}}, \code{\link{superPos}}, \code{\link{nse}}, \code{\link{rmse}}.
#'          \code{deconvolution.uh} in the package hydromad, \url{http://hydromad.catchment.org}
#' @references \url{http://ponce.sdsu.edu/onlineuhcascade.php}\cr
#'    Skript 'Abflusskonzentration' zur Vorlesungsreihe Abwasserentsorgung I von Prof. Krebs an der TU Dresden\cr
#'   \url{http://tu-dresden.de/die_tu_dresden/fakultaeten/fakultaet_forst_geo_und_hydrowissenschaften/fachrichtung_wasserwesen/isiw/sww/lehre/dateien/abwasserbehandlung/uebung_ws09_10/uebung_awe_1_abflusskonzentration.pdf}\cr
#'   \url{http://www.uni-potsdam.de/fs-g3/file.php?fileserver=klausuren&file=\%2FMaster_of_Science\%2FHydroII_Lernzettel.pdf}
#' @keywords hplot ts optimize
#' @importFrom graphics axis legend lines mtext par plot text title
#' @importFrom stats optim
#' @export
#' @examples
#' # Simplified for this course. True documentation is in berryFunctions::lsc
#' \donttest{ # this will not run directly, as the debugging exercise requires to fix the bugs
#' data(calib)
#' lsc(calib$P, calib$Q, area=1.6, quietNA=T)
#' }
#'
#' @param P       Vector with precipitation values \bold{in mm in hourly spacing}
#' @param Q       Vector with observed discharge (runoff) \bold{in m^3/s}
#'                with the same length as precipitation.
#' @param area    Single numeric. Catchment area \bold{in km^2}
#' @param Qbase   Baseflow that is added to UH-induced simulated Q
#'                DEFAULT: Q[1] thus cutting off baseflow in a very simple manner.
#' @param n       Numeric. Initial number of storages in cascade.
#'                Not necessarily integer. DEFAULT: 2
#' @param k       Numeric. Initial storage coefficient (resistance to let water run out).
#'                High damping, slowly reacting landscape, high k. DEFAULT: 3
#' @param x       Vector for the x-axis of the plot. DEFAULT: sequence along P
#' @param fit     Integer vector. Indices for a subset of Q that Qsim is fitted to.
#'                DEFAULT: all of Q
#' @param plot    Logical. plot input data? DEFAULT: TRUE
#' @param main    Character string. DEFAULT: "Precipitation and discharge"
#' @param plotsim Logical. add best fit to plot? DEFAULT: TRUE
#' @param returnsim Logical. Return simulated Q instead of parameters of UH?
#'                DEFAULT: FALSE
#' @param type    Vector with two characters: type as in \code{\link{plot}},
#'                repeated if only one is given. 1st for obs, 2nd for sim.
#'                DEFAULT: c("o","l")
#' @param legx    legend position. DEFAULT: "center"
#' @param legy    legend position. DEFAULT: NULL
#' @param \dots   Arguments passed to optim
#'
lsc <- function(P,
Q,
area=50,
Qbase=Q[1],
n=2,
k=3,
x=1:length(P),
fit=1:length(Q),
plot=TRUE,
main="Precipitation and discharge",
plotsim=TRUE,
returnsim=FALSE,
type=c("o", "l"),
legx="center",
legy=NULL,
...)
{
# checking for wrong input:
if(length(P) != length(Q)) stop("Vectors P and Q are not of equal length.")
if(any(fit>length(Q) | fit<0)) stop("'fit' has to contain positive integers smaller than length(Q).")
# initial parameters (few storages, quickly drained catchment):
param <- c(n=n, k=k)
stop("stupid error you can easily remove.")
# root mean square error (RMSE) of Qsim to Q ist to be minimized:
minfun <- function(param)
   { # discrete values of UH:
   UH_val <- unitHydrograph(n=param["n"], k=param["k"], t=1:length(P))
   qsim <- superPos(P/10^3, UH_val) * area*10^6 /3600 + Qbase
   rmse( Q[fit], qsim[fit])
   }
# do the hard work:
optimized <- optim(par=param, fn=minfun, ...)$par
# calculate optimized UH:
finalUH <- unitHydrograph(optimized["n"], optimized["k"], t=1:length(P))
if(round(sum(finalUH), 1) !=1) warning("sum of UH is not 1, probably the time should be longer")
# simulate runoff:
Qsim <- superPos(P/10^3,  finalUH) * area*10^6 /3600 + Qbase
Qsim <- Qsim[1:length(Q)]
#
# runoff coefficient Psi:
# psi*P * A = Q * t
# psi = Qt / PA  # UNITS:  m^3/s * h * 3600s/h  / (mm=1/1000 m^3/m^2 * km^2)  /  1e6 m^2/km^2
psi <- mean(Q-Qbase, na.rm=TRUE) * length(Q) * 3600  /  ( sum(P)/1000 * area) / 1e6
if(psi>1) warning("Psi is larger than 1. The area given is not able to yield so much discharge. Consider the units in  ?lsc")
#
# graphic:
if(plot)
  {
  if(length(type)==1) type <- rep(type,2)
  op <- par(mar=rep(3,4))
  plot(x, P, type="h", col=4, yaxt="n", ylim=c(max(P)*2, 0), lwd=2, ann=FALSE)
  title(main=main)
  axis(2, pretty(P), col=4, las=1, col.axis=4)
  #
  par(new=TRUE); plot(x, Q, type=type[1], col=2, las=1, ylim=range(Q)*c(1,2), ann=FALSE, axes=FALSE)
  axis(4, pretty(Q), col=2, las=1, col.axis=2)
  #
  mtext("P [mm]", line=-2, col=4, adj=0.02, outer=TRUE)
  mtext("Q [m^3/s]", line=-2, col=2, adj=0.98, outer=TRUE)
  #
  if(plotsim)
     {
     lines(x, Qsim, type=type[2], col=8, lwd=2)
     lines(x[fit], Qsim[fit], col=1, lwd=2)
     legend(legx, legy, legend=c("Observed","Simulated"), lwd=c(1,2), pch=c(1,NA), col=c(2,1) )
     }
  par(op)
  } # end if plot
#
if(returnsim) return(Qsim)
else return(c(n=as.vector(optimized["n"]), k=as.vector(optimized["k"]), NSE=nse(Q, Qsim), psi=psi))
} # end of funtion
