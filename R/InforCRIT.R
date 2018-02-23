#' InfoCRIT
#' 
#' This function determine values of the criteria Akaike's Information Criterion
#' (AIC), Bayesian Information Criterion (BIC), AICc with bias correction and  
#' to fitted probability density 
#' functions with lmomco package
#'
#' @param Intensity: numeric vector with intensity values for a duration 
#' specific in different return periods. 
#' @param Parameters: list with three elements: i) type of distribution 
#' function ii) parameters fitted, and iii) source to call specfic function
#' in the package \pkg{lmomco}.
#'
#' @return Provided a numeric vector with the corresponding AIC,  BIC, AICc and KIC values.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#' 
#' @references 
#' 
#' 
#' @export
#'
#' @examples
#' # Meteorology station in the Airport Farfan in Tulua, Colombia. 
#' data(Intgum5min)
#' data(Pargumbel)
#' InfoCRIT(Intensity=Intgum5min,Parameters=Pargumbel)
#' # Result: AIC= 86.73, BIC= 86.62, AICc= 89.73, KIC= 
#' 
#' 
#' 
#' 
InfoCRIT<- function (Intensity=...,Parameters=...){
  n <- length(Intensity)
  p <- length(Parameters$para)
  f <- lmomco::dlmomco(Intensity, Parameters)
  l <- log(f)
  Q <- lmomco::mle2par(Intensity, type = Parameters$type, hessian = TRUE)
  LogLik <- sum(l)
  AIC <- -2 * LogLik + 2 * p
  BIC <- -2 * LogLik + p * log(n)
  AICc = AIC + (2 * p * ( p + 1))/(n - p - 1);
  KIC <- LogLik - p * log(2*pi) - log(det(Q$optim$hessian))
  
  out <- c(AIC, BIC, AICc)
  names(out) <- c("AIC", "BIC", "AICc", "KIC")
  return(out)
}
