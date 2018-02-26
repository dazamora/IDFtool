#' pW Criteria
#' 
#' This function determine what probability distribution function has the best
#' goodness-of-fit to observations, number of parameters and the quality of 
#' taking in account different information criterias that are evaluated in 
#' the criteria proposed by (\emph{Siena et al., 2017}).
#' 
#'
#' @param metrics: a numeric matrix with the values of criterias AIC, BIC, AICc and KIC (rows)
#' evaluated by each PDFs (columns). 
#' @param critnames: a character vector with the names of information criteria evaluated.
#' @param pdfnames: a character vector with the names of PDFs evaluated. 
#'
#' @return Provided a matrix with values of pW criteria by each PDFs and information criterias.
#' 
#' @author Adriana Pina <appinaf@unal.edu.co> and David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#' 
#' @references 
#' 
#' 
#' @export
#'
#' @examples
#' # Example of five PDFs and their respectives values of four information criterias.
#' data(fractures.crit)
#' 
#' pW.1 <- p.criteria(metrics = fractures.crit, critnames = rownames(fractures.crit), pdfnames = colnames(fractures.crit))
#' pW.1
#' 
#' head(round(pW.1, 2))
#' # GEV: AIC= 1, BIC= 1, AICc= 1, KIC= 1 
p.criteria <- function(metrics = .., critnames = .., pdfnames = ..){
  
  if(is.matrix(metrics)){
    IC <- metrics 
  } else {
    stop("metrics is not a matrix")
  }
  
  pMk <- 1/dim(metrics)[2] 
  MinIC <- apply(IC,1,min)
  
  DIC <- IC - matrix(MinIC, nrow = dim(IC)[1], ncol = dim(IC)[2])
  aux <- exp(-0.5 * DIC * pMk)
  SumDICi <- rowSums(aux)
  pW <- aux/matrix(SumDICi, dim(IC)[1], dim(IC)[2])
  
  colnames(pW) <- pdfnames
  rownames(pW) <- critnames
  
  return(pW)
}

