#' Title
#'
#' @param Intensity 
#' @param Parameters 
#'
#' @return
#' @export
#'
#' @examples
#' 
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
