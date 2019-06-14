#' @name 
#' resoPLOT
#' 
#' @title
#' Selection of the graphs resolution 
#' 
#' @param reso: a number to determine the resolution that the plot function will used to save 
#' graphs. It has two options: 300 and 600 ppi. 
#'
#' @return a vector with four values. They are used 
#' to save graphs by means of \code{\link{png}} function.
#' 
#' \itemize{
#'  \item \code{reso}
#'  \item \code{wth}
#'  \item \code{hth}
#'  \item \code{hth1}
#' }
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#' 
#' @export
#'
resoPLOT <- function(reso = 300){

  if (reso == 300) {
    reso <- 300
    wth <- 3*580
    hth <- wth
    hth1 <- wth/2^(.5)
  } else {
    reso <- 600
    wth <- 6*580
    hth <- wth
    hth1 <- wth/2^(.5)
  }

  out1 <- c(reso,wth,hth,hth1)
  names(out1) <- c("reso","wth","hth","hth1")
  return(out1)

}
