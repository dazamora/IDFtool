#' InfoCRIT
#' 
#' This function determine values of the criteria Akaike's Information Criterion
#' (AIC) and Bayesian Information Criterion to fitted probability density 
#' functions with lmomco package
#'
#' @param Intensity: numeric vector with intensity values for a duration 
#' specific in different return periods. 
#' @param Parameters: list with three elements: i) type of distribution 
#' function ii) parameters fitted, and iii) source to call specfic function
#' in the package \pkg{lmomco}.
#'
#' @return Provided a numeric vector with the corresponding AIC and BIC values.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#' 
#' @export
#'
#' @examples
#' # Meteorology station in the Airport Farfan in Tulua, Colombia. 
#' data(Intgum5min)
#' data(Pargumbel)
#' InfoCRIT(Intensity=Intgum5min,Parameters=Pargumbel)
#' # Result: AIC= 86.73 & BIC= 86.62
#' 
InfoCRIT<- function (Intensity=...,Parameters=...){
  n <- length(Intensity)
  p <- length(Parameters$para)
  f <- lmomco::dlmomco(Intensity, Parameters)
  l <- log(f)
  LogLik <- sum(l)
  AIC <- -2 * LogLik + 2 * p
  BIC <- -2 * LogLik + p * log(n)
  out <- c(AIC, BIC)
  names(out) <- c("AIC", "BIC")
  return(out)
}
