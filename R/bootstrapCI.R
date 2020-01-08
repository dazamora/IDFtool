#' @name
#' bootstrapCI
#' 
#' @title 
#' Determine confidence intervals to rainfall intensities estimated by a PDF function.
#' 
#' @description 
#' Conducts bootstrap to randomly sample of intensity values 'n' times for a 
#' specified distribution to estimate the confidence interval for each given 
#' non-exceedance probability.
#'
#' @param Intensity a numeric vector with intensity [mm/h] values of different 
#' years for a specific time duration (\emph{e.g.} 5, 15, 120 minutes, \emph{etc.})
#' @param Parameters list with three elements: i) type of distribution function 
#' ii) fitted parameters, and iii) source to call specfic function in the \pkg{lmomco} package.
#' @param Type a character specifying a name of the probability distribution function fitted 
#' (see \code{\link{selecDIST}}) by \code{\link{fitDISTRI}} function.
#' @param Rsample An integer representing number of resamples to conduct when 
#' confidence interval will be computed.
#' @param Return.P a numeric vector with return periods like non-exceedance probabilities.
#' @param Conf.Inter level of the confidence interval.
#'
#' @return A list of:
#' 
#'  \itemize{
#'   \item \code{nonexceed.prob} a numeric vector with non-exceedance probabilities.
#'   \item \code{lower.lim} a numeric vector confidence bound lower for quantile estimates.
#'   \item \code{upper.lim} a numeric vector confidence bound upper for quantile estimates.
#'   \item \code{Para.set} a matrix containing estimated distribution parameters for each resample.
#'   \item \code{quantiles} a matrix of quantile estimates for each resample.
#'  }
#'  
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#' 
#' @export
#'
#' @examples 
#' 
#' data(inten)
#' data(Pargumbel)
#' Tp <- c(2, 3, 5, 10, 25, 50, 100)
#' FR <- 1 - 1/Tp
#' CI.test <- bootstrapCI(Intensity = inten, Parameters = Pargumbel, Type = "gumbel",
#' Rsample = 50, Return.P =FR, Conf.Inter = 0.90)
#' (CI.test$Conf.Inter[,2])
#' # result = 125.9501, 139.6132, 154.6232, 173.6736, 197.5024, 215.0634, 232.4939
#'  
bootstrapCI <- function(Intensity, Parameters ,Type = "Gumbel",
                        Rsample = 1000, Return.P , Conf.Inter = 0.95) {
  
  # extract fitted model parameters
  base.params <- Parameters
  
  # create output matrices to store parameter sets and quantile estimates
  param.sets <- matrix(NA, nrow = Rsample, ncol = length(base.params$para))
  quantile.estimates <- matrix(NA, nrow = Rsample, ncol = length(Return.P),
                               dimnames = list(NULL, Return.P) )
  
  # begin bootstrapping procedure
  for(i in 1:Rsample) {
    
    valid.moments <- FALSE
    j <- 0
    
    # allow up to 20 re-tries to re-sample
    while(!valid.moments & j < 20) {
      
      # sample 'n' random variates from base distribution
      data <- lmomco::rlmomco(n=length(Intensity), base.params)
      
      # compute sample l-moments
      sample.moms <- lmomco::lmoms(data)
      
      valid.moments <- lmomco::are.lmom.valid(sample.moms)
      j <- j + 1
    }
    
    # error handling
    if(!valid.moments) {
      stop("Bootstrapping failed to sample valid l-moments")
    } else {
      # estimate distribution parameters
      dist.par <- lmomco::lmom2par(sample.moms, base.params$type)
      
      # store the distribution parameters
      param.sets[i,] <- dist.par$para
      
      # estimate quantiles at NEP
      estimated <- lmomco::qlmomco(Return.P, dist.par)
      
      # convert quantile estimates to real values if
      # distribution was transformed
      if(Type=="log.pearson3") estimated <- 10^estimated
      
      # store the quantiles at the desired AEP values
      quantile.estimates[i,] <- estimated
    }
    
  }
  
  # now calculate confidence limits for quantiles
  p <- c((1-Conf.Inter)/2, (1+Conf.Inter)/2)
  Conf.Inter <- sapply(colnames(quantile.estimates),
                       FUN=function(x){
                         stats::quantile(quantile.estimates[,x], probs=p, na.rm=TRUE)})
  
  # now return list object containing output
  return(list(Conf.Inter = data.frame(
    nonexceed.prob = Return.P,
    lower.lim = as.vector(Conf.Inter[1,]),
    upper.lim = as.vector(Conf.Inter[2,])),
    Para.set = param.sets,
    quantiles = quantile.estimates)
  )

}
