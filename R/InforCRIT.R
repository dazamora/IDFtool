#' @name
#' InfoCRIT
#' 
#' @title 
#' Generic function calculating information criterias for fitted PDF
#' 
#' @description 
#' This function determine values of the criteria Akaike's Information Criterion
#' (AIC), Bayesian Information Criterion (BIC), AICc with bias correction and KIC 
#' to fitted probability density functions with \pkg{lmomco} package
#'
#' @param Intensity numeric vector with intensity values for a specific duration 
#' in different return periods. 
#' @param Parameters list with three elements: i) type of distribution 
#' function ii) fitted parameters, and iii) source to call specfic function
#' in the package \pkg{lmomco}.
#' @param ... Additional arguments for the \code{optim()} function and other uses
#' 
#' @return Provided a numeric vector with the corresponding AIC,  BIC, AICc and KIC values.
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
#' crit  <- InfoCRIT(Intensity = Intgum5min, Parameters = Pargumbel)
#' ## Result: AIC= 86.73, BIC= 86.62, AICc= 89.73, KIC= -35.23 
InfoCRIT <- function (Intensity, Parameters, ...){
  n <- length(Intensity)
  p <- length(Parameters$para)
  f <- lmomco::dlmomco(Intensity, Parameters)
  l <- log(f)
  
  mledz.optim <- function (x, type, para.int = NULL, silent = TRUE, null.on.not.converge = TRUE, 
                           ptransf = function(t) return(t), pretransf = function(t) return(t), 
                           ...) 
  {
    if (is.null(para.int)) {
      lmr <- lmomco::lmoms(x)
      para.int <- lmomco::lmom2par(lmr, type = type, ...)
    }
    if (is.null(para.int)) {
      warning("could not estimate initial parameters via L-moments")
      return(NULL)
    }
    if (length(para.int$para) == 1) {
      warning("function is not yet built for single parameter optimization")
      return(NULL)
    }
    "afunc" <- function(para, x = NULL, ...) {
      lmomco.para <- lmomco::vec2par(pretransf(para), type = type, 
                                     paracheck = TRUE)
      if (is.null(lmomco.para)) 
        return(Inf)
      pdf <- lmomco::par2pdf(x, lmomco.para)
      L <- -sum(log(pdf), na.rm = TRUE)
      return(L)
    }
    rt <- NULL
    try(rt <- stats::optim(par = ptransf(para.int$para), fn = afunc, 
                           x = x,...), silent = silent)
    if (is.null(rt)) {
      warning("optim() attempt is NULL")
      return(NULL)
    }
    else {
      if (null.on.not.converge & rt$convergence != 0) {
        warning("optim() reports convergence error")
        return(NULL)
      }
      lmomco.para <- lmomco::vec2par(pretransf(rt$par), type = type)
      lmomco.para$AIC <- 2 * length(rt$par) - 2 * (-1 * rt$value)
      lmomco.para$optim <- rt
      return(lmomco.para)
    }
  }
  Q <- mledz.optim(x=Intensity, type = Parameters$type, hessian = TRUE)
  LogLik <- sum(l)
  AIC <- -2 * LogLik + 2 * p
  BIC <- -2 * LogLik + p * log(n)
  AICc = AIC + (2 * p * ( p + 1))/(n - p - 1);
  KIC <- LogLik - p * log(2*pi) - log(det(Q$optim$hessian))
  
  out <- c(AIC, BIC, AICc, KIC)
  names(out) <- c("AIC", "BIC", "AICc", "KIC")
  return(out)
}
