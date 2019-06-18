#' @name 
#' MLEZ
#' 
#' @title 
#' Fitting parameters of PDF by means of maximum likelihood
#' 
#' @description 
#' This function fits parameters of probability distribution functions (see \code{\link{selecDIST}})
#' by means of maximum likelihood and performs evolutionary global optimization via the 
#' Differential Evolution algorithm (see \pkg{DEoptim} package).
#'
#' @param Intensity a numeric vector with intensity [mm/h] values of different
#' years for a specific time duration (\emph{e.g.} 5, 15, 120 minutes, \emph{etc}.).  
#' @param type a character specifying the name of the distribution function that will 
#' be employed: exponencial, gamma, gev, gumbel, log.normal3, normal, pearson, log.pearson3 and 
#' wakeby (see \code{\link{selecDIST}}).
#' @param para.int Initial parameters as a vector \eqn{\Theta}.
#' @param silent a logical to silence the \code{\link{try}} function wrapping the \code{\link{DEoptim}} function.
#' @param null.on.not.converge a logical to trigger simple return of NULL if the \code{\link{DEoptim}} function 
#' returns a nonzero convergence status.#' 
#' @param pretransf
#'
#' @return A list of:
#' 
#'  \itemize{
#'    \item \code{Parameters} a list with type of distribution fitted and values of its parameters
#'  }
#' 
#' @export
#'
#' @examples
#' 
#' data(inten) 
#' TEST.MLE <- MLEZ(Intensity = inten[,4], type = "Gumbel", 
#'                  para.int = NULL, silent = TRUE, 
#'                  null.on.not.converge = TRUE, 
#'                  pretransf = function(t) return(t))
#' ## Results: xi = 71.178 ; alpha = 15.204 
#' 
MLEZ <- function (Intensity, type, para.int = NULL, silent = TRUE, null.on.not.converge = TRUE, 
                 pretransf = function(t) return(t)) {
  
  x <- Intensity
  type <- tolower(type)
  type <- selecDIST(Type = type)
  
  if (is.null(para.int)) {
    lmr <- lmomco::lmoms(x)
    if(lmr$ratios[3] < 0 & type == "ln3"){
      lmr$ratios[3] <- -lmr$ratios[3]
    }
    para.int <- lmomco::lmom2par(lmr, type = type)
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
  lower <- total.par*0.7*(-1)
  upper <- total.par*1.3
  try(rt <- DEoptim::DEoptim(afunc, lower, upper, x = x, 
                    control = DEoptim::DEoptim.control(trace = FALSE, packages = "lmomco")),silent = silent)
  
  if (is.null(rt)) {
    warning("optim() attempt is NULL")
    return(NULL)
  } else {
    lmomco.para <- lmomco::vec2par(pretransf(rt$optim$bestmem), type = type)
    lmomco.para$AIC <- 2 * length(rt$optim$bestmem) - 2 * (-1 * rt$value)
    #lmomco.para$optim <- rt
    return(lmomco.para)
  }
}
