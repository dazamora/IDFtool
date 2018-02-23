#' pW Criteria
#' 
#' This function determine what probability distribution function has the best
#' goodness-of-fit to observations, number of parameters and the quality of 
#' model parameter estimates, taking in account different criterias information in a 
#' \emph(Siena et al., 2017). 
#' 
#'
#' @param metrics 
#' @param numod 
#' @param pdfnames 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
p.criteria <- function(metrics = .., numod = .., pdfnames = ..){
  
  if(is.matrix(metrics)){
    IC <- metrics 
  } else {
    stop("metrics is not a matrix")
  }
  
  pMk <- 1/dim(metrics)[1] 
  MinIC <- apply(IC,2,min)
  
  DIC <- IC - matrix(MinIC, nrow = pMk, ncol = dim(IC)[2])
  aux <- exp(-0.5 * DIC * pMk)
  SumDICi <- sum(aux,2)
  pW <- aux/rowsum(SumDICi)
  
  return(pW)
}