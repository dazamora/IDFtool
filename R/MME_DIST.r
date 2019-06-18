#' @name 
#' MME_DIST
#' 
#' @title 
#' Fitting parameters of PDF by means of moments method
#' 
#' @description 
#' This function allows to fit several distribution functions to observed data
#' by means of the moments method.
#'
#' @param Intensity: a numeric vector with intensity [mm/h] values of different
#' years for a specific time duration (\emph{e.g.} 5, 15, 120 minutes, \emph{etc}.).
#' @param Type: a character specifying the name of distribution function that it will 
#' be employed: exponencial, gamma, gev, gumbel, log.normal3, normal, pearson, log.pearson3 and 
#' wakeby (see \code{\link{selecDIST}}). 
#'
#' @return A list of:
#' 
#'  \itemize{
#'    \item \code{Parameters} a list with type of distribution fitted and values of its parameters
#'  }
#'  
#' @author Albeiro Figueroa <cafigueroao@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#'   
#' @export
#' 
#' @examples
#' 
#' # Meteorology station in the Farfan Airport in Tulua, Colombia.
#' data(inten)
#' Test.moments <- MME_DIST(Intensity = inten[,4], Type = "gumbel") 
#' ## Results: xi = 71.50 ; alpha = 14.21
#' 
MME_DIST <- function (Intensity, Type) {
  
  data <- Intensity
  distr <- selecDIST(Type)
  
  if (!is.character(distr)) {
    stop("distr must be a character string naming a distribution")
  }else{ distname <- distr }
  
  if (is.element(distname, c("nor", "ln3", "gam", "exp", "pe3", "lpe3","gum", "gev"))){ 
    meth <- "closed formula"
  }else{ stop("distribution not allowed") }
  
  if (!(is.numeric(data) & length(data) > 1)) {
    stop("data must be a numeric vector of length greater than 1")
  }
  n <- length(data)
  m <- mean(data)
  v <- var(data)
  ##################################################################
  ##   Distribución Normal (2P - norm),
  if (distname == "nor") {
    estimate <- c(mean = m, sd = sqrt(v))
    order <- 1:2
    PAR <- lmomco::vec2par(c(location = m, scale = sqrt(v)), type="nor")
    return(PAR)
  }
  #################################################################################
  ##   Distribución Log Normal (3P - ln3), from EnvStats with Method of Moments Estimation
  if (distname == "ln3") {
    if (any(data <= 0)){
      stop("values must be positive to fit a lognormal distribution")
    }
    
    m2 <- ((n - 1)/n) * v
    sqrt.b1 <- EnvStats::skewness(data, method = "moment")
    if (sqrt.b1 <= 0){
      print(paste("The sample skew is not positive. ",
                  "Admissible moment estimates do not exist"))
      PAR <- NULL
    } else {
      b1 <- sqrt.b1^2
      t1 <- 1 + b1/2
      t2 <- sqrt(t1^2 - 1)
      omega <- (t1 + t2)^(1/3) + (t1 - t2)^(1/3) - 1
      varlog <- log(omega)
      sdlog <- sqrt(varlog)
      meanlog <- 0.5 * log(m2/(omega * (omega - 1)))
      threshold <- m - exp(meanlog + varlog/2)
      
      PAR <- lmomco::vec2par(c(zeta = threshold,  mu = meanlog, sigma = sdlog), type = "ln3")
      # is.ln3(PAR)
    }
    return(PAR)
  }
  ##################################################################
  ##   Distribución Gamma (2P - gamma),
  if (distname == "gam") {
    shape <- m^2/v
    scale <- v/m
    PAR <- lmomco::vec2par(c(alpha=shape, beta=scale), type="gam")
    return(PAR)
  }
  ##################################################################
  ##   Distribución Exponencial (2P- exp), 
  if (distname == "exp") {
    location <- m - sqrt(v)
    scale <- sqrt(v)
    PAR <- lmomco::vec2par(c(xi=location, alpha=scale), type="exp")
    return(PAR)
  }
  ##################################################################
  ##   Pearson Tipo III (3P - pe3), 
  if (distname == "pe3") { 
    skew <- mean((data - mean(data))^3) / (sqrt(v)^3)
    shape  <- (4 / skew^2)
    scale <- sqrt(v) / sqrt(shape)
    location <- m - scale * shape
    PAR <- lmomco::vec2par(c(mu=location, sigma=scale, gamma=shape), type="pe3")
    return(PAR)
  }
  ##################################################################
  ##   Log Pearson Tipo III (3P - lpe3), 
  if (distname == "lpe3") { 
    skew <- mean((log(data) - mean(log(data)))^3) / (sqrt(var(log(data)))^3)
    shape  <- (4 / skew^2)
    scale <- sqrt(var(log(data))) / sqrt(shape)
    location <- mean(log(data)) - scale * shape
    PAR <- lmomco::vec2par(c(mu=location, sigma=scale, gamma=shape), type="pe3")
    return(PAR)
  }
  ##################################################################
  ##   Gumbel o Extremos Tipo I(2P - gum),
  if (distname == "gum") {
    scale <- (6^0.5) / pi *  sqrt(v)
    location <- m - (0.5772 * scale)
    PAR <- lmomco::vec2par(c(xi = location, alpha = scale), type="gum")
    return(PAR)
  }	
  ##################################################################
  ##   Distribución Valores Extremos Generalizados (3P - gev),
  if (distname == "gev"){
    scale <- sqrt(6 * var(data))/pi
    location <- mean(data) - 0.58 * scale
    k <- round(0.5 * n)
    mind <- min(data)
    if (mind < 0) 
      data <- data - mind
    data <- sort(data)
    delta <- log(data[n - 0:(k - 1)]/data[n - k])
    M1 <- sum(delta)/k
    M2 <- sum(delta^2)/k
    shape <- M1 + 1 - 0.5 * (M2/(M2 - M1^2))
    vec <- c(xi = location, alpha = scale, kappa = shape)
    z <- list(type = "gev", para = vec, source = "vec2par")
    names(z$para) <- c("xi", "alpha", "kappa")
    # PAR <- lmomco::vec2par(c(xi = location, alpha = scale, kappa = shape), nowarn = TRUE,type="gev")
    PAR <- z
    return(PAR)
  }
}
