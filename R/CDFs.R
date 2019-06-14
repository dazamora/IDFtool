#' @name 
#' pexp1
#' 
#' @title
#' Cumulative probability of Exponential distribution
#' 
#' @description     
#' This function computes the cumulative probability or nonexceedance probability
#' of the Exponential distribution given parameters \eqn{\xi} and \eqn{\alpha}. Based in function 
#' of \pkg{lmomco} package.
#'
#' @param q: A real value vector.
#' @param xi: \eqn{\xi} is a location parameter (see \code{\link{parexp}} function).
#' @param alpha: \eqn{\alpha} s a scale parameter (see \code{\link{parexp}} function).
#'
#' @return Nonexceedance probability (\emph{F}) for q.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#'  
#' @export
#'
#' @examples
#' 
pexp1 <- function(q, xi, alpha) {
  # ----Exponencial----
  lmomco::cdfexp(q, list(type = "exp", para = c(xi, alpha), source = "parexp"))
}

#' @name 
#' pgam1
#' 
#' @title 
#' Cumulative probability of Gamma distribution
#' 
#' @description 
#' This function computes the cumulative probability or nonexceedance probability 
#' of the Gamma distribution given parameters (\eqn{\alpha} and \eqn{\beta}) 
#' computed. Based in function of \pkg{lmomco} package.
#'
#' @param q: a real value vector.
#' @param alpha: \eqn{\alpha} is a scale parameter (see \code{\link{pargam}} function).
#' @param beta: \eqn{\beta} is a location parameter (see \code{\link{pargam}} function).
#'
#' @return Nonexceedance probability (\emph{F}) for q.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIRE
#' 
#' @export
#'
#' @examples
#' 
pgam1 <- function(q, alpha, beta) {
  # ----Gamma----
  lmomco::cdfgam(q, list(type = "gam", para = c(alpha, beta), source = "pargam"))
}

#' @name
#' pgev1
#' 
#' @title 
#' Cumulative probability of GEV distribution
#' 
#' @description 
#' This function computes the cumulative probability or nonexceedance probability 
#' of the Generalized Extreme Value distribution given parameters (\eqn{\xi}, \eqn{\alpha} and \eqn{\kappa}). 
#' Based in function of \pkg{lmomco} package.
#'
#' @param q: a real value vector.
#' @param xi: \eqn{\xi} is a location parameter (see \code{\link{pargev}} function).
#' @param alpha: \eqn{\alpha} is a scale parameter (see \code{\link{pargev}} function).
#' @param kappa: \eqn{\kappa} is a shape parameter (see \code{\link{pargev}} function).
#'
#' @return Nonexceedance probability (\emph{F}) for q.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIRE
#' 
#' @export
#'
#' @examples
#' 
pgev1 <- function(q, xi, alpha, kappa) {
  # ----GEV----
  lmomco::cdfgev(q, list(type = "gev", para = c(xi, alpha, kappa), source = "pargev"))
}

#' @name 
#' pgum1
#' 
#' @title 
#' Cumulative probability of Gumbel distribution
#' 
#' @description 
#' This function computes the cumulative probability or nonexceedance probability 
#' of the Gumbel distribution given parameters (\eqn{\xi} and \eqn{\alpha}). Based in function of 
#' lmomco package.
#'
#' @param q: a real value vector.
#' @param xi: \eqn{\xi} is a location parameter (see \code{\link{pargum}} function).
#' @param alpha: \eqn{\alpha} is a scale paramete (see \code{\link{pargum}} function).
#'
#' @return Nonexceedance probability (\emph{F}) for q.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIRE
#'  
#' @export
#'
#' @examples
#' 
pgum1 <- function(q, xi, alpha) {
  # ----Gumbel----
  lmomco::cdfgum(q, list(type = "gum", para = c(xi, alpha), source = "pargum"))
}


#' @name 
#' pln31
#' 
#' @title 
#' Cumulative probability of LogNormal distribution
#' 
#' @description 
#' This function computes the cumulative probability or nonexceedance probability 
#' of the Log-Normal3 distribution given parameters (\eqn{\zeta}, \eqn{\mu} and \eqn{\sigma}). Based in function 
#' of \pkg{lmomco} package.
#'
#' @param q: a real value vector.
#' @param zeta: \eqn{\zeta} is a lower bounds (see \code{\link{parln3}} function).
#' @param mulog: \eqn{\mu} is a location parameter (see \code{\link{parln3}} function).
#' @param sigmalog: \eqn{\sigma} is a scale parameter (see \code{\link{parln3}} function).
#' 
#' @return Nonexceedance probability (\emph{F}) for q.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIRE
#'
#' @export
#'
#' @examples
#' 
pln31 <- function(q, zeta, mulog, sigmalog) {
  # ----Log Normal 3----
  lmomco::cdfln3(q, list(type = "ln3", para = c(zeta, mulog, sigmalog), source = "parln3"))
}

#' @name 
#' pnor1
#' 
#' @title 
#' Cumulative probability of Normal distribution
#' 
#' @description 
#' This function computes the cumulative probability or nonexceedance 
#' probability of the Normal distribution given parameters (\eqn{\mu} and \eqn{\sigma}). 
#' Based in function of \pkg{lmomco} package.
#'
#' @param q: a real value vector.
#' @param mu: \eqn{\mu} is the arithmetic mean (see \code{\link{parnor}} function).
#' @param sigma: \eqn{\sigma} is the standard deviation (see \code{\link{parnor}} function).
#'
#' @return Nonexceedance probability (\emph{F}) for q.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIRE
#'
#' @export
#'
#' @examples
#' 
pnor1 <- function(q, mu, sigma) {
  # ----Normal----
  lmomco::cdfnor(q, list(type = "nor", para = c(mu, sigma), source = "parnor"))
}

#' @name 
#' ppe31
#' 
#' @title 
#' Cumulative probability of Pearson Type III distribution
#' 
#' @description 
#' This function computes the cumulative probability or nonexceedance 
#' probability of the Pearson Type III distribution given parameters 
#' (\eqn{\mu}, \eqn{\sigma}, and \eqn{\gamma}). Based in function of \pkg{lmomco} package.
#'
#' @param q: a real value vector.
#' @param mu: \eqn{\mu} is the mean (see \code{\link{parexp}} function).
#' @param sigma: \eqn{\sigma} is the standard deviation (see \code{\link{parexp}} function).
#' @param gamma: \eqn{\gamma} is the skew (see \code{\link{parexp}} function).
#'
#' @return Nonexceedance probability (\emph{F}) for q.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIRE
#'
#' @export
#'
#' @examples
#' 
ppe31 <- function(q, mu, sigma, gamma) {
  # ----Person III----
  lmomco::cdfpe3(q, list(type = "pe3", para = c(mu, sigma, gamma), source = "parpe3"))
}

#' @name 
#' pwak1
#' 
#' @title 
#' Cumulative probability of Wakeby distribution
#' 
#' @description 
#' This function computes the cumulative probability or nonexceedance 
#' probability of the Wakeby distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma} and \eqn{\delta}). 
#' Based in function of \pkg{lmomco} package.
#'
#' @param q: a real value vector.
#' @param xi: \eqn{\xi} is a location parameter (see \code{\link{parwak}} function).
#' @param alpha: \eqn{\alpha} is a location parameter (see \code{\link{parwak}} function).
#' @param beta: \eqn{\gamma} is a shape parameter (see \code{\link{parwak}} function).
#' @param gamma: \eqn{\gamma} is a shape parameter (see \code{\link{parwak}} function).
#' @param delta: \eqn{\delta} is a shape parameter (see \code{\link{parwak}} function).
#'
#' @return Nonexceedance probability (\emph{F}) for q.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIRE
#'
#' @export
#'
#' @examples
#' 
pwak1 <- function(q, xi, alpha, beta, gamma, delta) {
  # ----Wakeby----
  lmomco::cdfwak(q, list(type = "wak", para = c(xi, alpha, beta, gamma, delta), source = "parwak"))
}