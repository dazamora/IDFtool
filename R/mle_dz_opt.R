mledz.optim <-function (x, type, para.int = NULL, silent = TRUE, null.on.not.converge = TRUE, 
          ptransf = function(t) return(t), pretransf = function(t) return(t), 
          ...) 
{
  if (is.null(para.int)) {
    lmr <- lmoms(x)
    para.int <- lmom2par(lmr, type = type, ...)
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
    lmomco.para <- vec2par(pretransf(para), type = type, 
                           paracheck = TRUE)
    if (is.null(lmomco.para)) 
      return(Inf)
    pdf <- par2pdf(x, lmomco.para)
    L <- -sum(log(pdf), na.rm = TRUE)
    return(L)
  }
  rt <- NULL
  try(rt <- optim(par = ptransf(para.int$para), fn = afunc, 
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
    lmomco.para <- vec2par(pretransf(rt$par), type = type)
    lmomco.para$AIC <- 2 * length(rt$par) - 2 * (-1 * rt$value)
    lmomco.para$optim <- rt
    return(lmomco.para)
  }
}