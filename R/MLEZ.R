#' Title
#'
#' @param x 
#' @param type 
#' @param para.int 
#' @param silent 
#' @param null.on.not.converge 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
MLEZ <- function (x, type, para.int = NULL, silent = TRUE, null.on.not.converge = TRUE, 
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
  
  total.par <- abs(para.int$para)
  lower <- para.int$para*0.7*(-1)
  upper <- para.int$para*1.3
  try(rt <- DEoptim(afunc, lower, upper, x = x, 
                    control = DEoptim.control(NP = 50, itermax = 100, trace = FALSE, packages = "lmomco")),silent = silent)
  
  if (is.null(rt)) {
    warning("optim() attempt is NULL")
    return(NULL)
  }
  else {
    if (null.on.not.converge & rt$convergence != 0) {
      warning("optim() reports convergence error")
      return(NULL)
    }
    lmomco.para <- lmomco::vec2par(pretransf(rt$optim$bestmem), type = type)
    lmomco.para$AIC <- 2 * length(rt$optim$bestmem) - 2 * (-1 * rt$value)
    lmomco.para$optim <- rt
    return(lmomco.para)
  }
}
