#' @name 
#' selecDIST
#' 
#' @title 
#' Convert large names of distribution functions in a short one
#' 
#' @description 
#' This function converts large names of distribution functions in a short one 
#' to used it in functions of \pkg{IDFtool} package. 
#'
#' @param Type a character specifying the name of the distribution function that 
#' will be employed: exponential, gamma, gev, gumbel, log.normal3, normal, 
#' log.pearson3 and wakeby. 
#'
#' @return a character with an abbreviation of distribution function defined in \code{Type}.
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#' 
#' @export
#'
#' @examples
#' (selecDIST(Type = "log.normal3"))
#' # "ln3
#' 
selecDIST <- function(Type = "Gumbel"){
  if (Type == "exponential") {
    distribution <- "exp"
  } else if (Type == "gamma") {
    distribution <- "gam"
  } else if (Type == "gev") {
    distribution <- "gev"
  } else if (Type == "gumbel") {
    distribution <- "gum"
  } else if (Type == "log.normal3") {
    distribution <- "ln3"
  }else if (Type == "normal") {
    distribution <- "nor"
  }else if (Type == "pearson3") {
    distribution <- "pe3"
  }else if (Type == "log.pearson3") {
    distribution <- "pe3"
  }else if (Type == "wakeby") {
    distribution <- "wak"
  }else{ #log normal 2P
    print("En desarrollo")
  }
  
  return(distribution)
}
