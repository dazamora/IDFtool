

#' Title
#'
#' @param Intensity 
#' @param Parameters 
#' @param Type 
#' @param Rsample 
#' @param Return.P 
#' @param Conf.Inter 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
bootstrapCI <- function(Intensity =.., Parameters =..,Type =..,
                        Rsample = 1000, Return.P =..., Conf.Inter = 0.90) {

  # extract fitted model parameters and flag as to whether the
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
                         quantile(quantile.estimates[,x], probs=p, na.rm=TRUE)})

  # now return list object containing output
  return(list(Conf.Inter = data.frame(
    nonexceed.prob = Return.P,
    lower.lim = as.vector(Conf.Inter[1,]),
    upper.lim = as.vector(Conf.Inter[2,])),
    Para.set = param.sets,
    quantiles = quantile.estimates)
  )

}
