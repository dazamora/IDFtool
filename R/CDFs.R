

#' Title
#'
#' @param q 
#' @param xi 
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
pexp1 <- function(q, xi, alpha) {
  # ----Exponencial----
  lmomco::cdfexp(q, list(type = "exp", para = c(xi, alpha), source = "parexp"))
}

#' Title
#'
#' @param q 
#' @param alpha 
#' @param beta 
#'
#' @return
#' @export
#'
#' @examples
pagm1 <- function(q, alpha, beta) {
  # ----Gamma----
  lmomco::cdfgam(q, list(type = "gam", para = c(alpha, beta), source = "pargam"))
}

#' Title
#'
#' @param q 
#' @param xi 
#' @param alpha 
#' @param kappa 
#'
#' @return
#' @export
#'
#' @examples
pgev1 <- function(q, xi, alpha, kappa) {
  # ----GEV----
  lmomco::cdfgev(q, list(type = "gev", para = c(xi, alpha, kappa), source = "pargev"))
}

#' Title
#'
#' @param q 
#' @param xi 
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
pgum1 <- function(q, xi, alpha) {
  # ----Gumbel----
  lmomco::cdfgum(q, list(type = "gum", para = c(xi, alpha), source = "pargum"))
}


#' Title
#'
#' @param q 
#' @param zeta 
#' @param mulog 
#' @param sigmalog 
#'
#' @return
#' @export
#'
#' @examples
pln31 <- function(q, zeta, mulog, sigmalog) {
  # ----Log Normal 3----
  lmomco::cdfln3(q, list(type = "ln3", para = c(zeta, mulog, sigmalog), source = "parln3"))
}

#' Title
#'
#' @param q 
#' @param mu 
#' @param sigma 
#'
#' @return
#' @export
#'
#' @examples
pnor1 <- function(q, mu, sigma) {
  # ----Normal----
  lmomco::cdfnor(q, list(type = "nor", para = c(mu, sigma), source = "parnor"))
}

#' Title
#'
#' @param q 
#' @param mu 
#' @param sigma 
#' @param gamma 
#'
#' @return
#' @export
#'
#' @examples
ppe31 <- function(q, mu, sigma, gamma) {
  # ----Person III----
  lmomco::cdfpe3(q, list(type = "pe3", para = c(mu, sigma, gamma), source = "parpe3"))
}

#' Title
#'
#' @param q 
#' @param xi 
#' @param alpha 
#' @param beta 
#' @param gamma 
#' @param delta 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
pwak1 <- function(q, xi, alpha, beta, gamma, delta) {
  # ----Wakeby----
  lmomco::cdfwak(q, list(type = "wak", para = c(xi, alpha, beta, gamma, delta), source = "parwak"))
}







