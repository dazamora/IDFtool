#' @name 
#' goodFIT
#'
#' @title
#' Computing of goodness-of-fit metrics and information criterias
#' 
#' @description    
#' This function computing of goodness-of-fit for continuous univariate 
#' distributions using tests: Kolmogorov-Smirnov, Anderson-Darling and 
#' Cramer-von Mises. It is based on \pkg{goftest} package. Moreover information criteria 
#' are evaluated: Akaike's Information Criterion and Bayesian Information 
#' Criterion, Akaike's Information Criterion with bias correction and Kashyap bayesian Information Criterio
#' by means of \code{\link{InfoCRIT}} function.
#'
#' @param Station: a character specifying a name or number of pluviographic 
#' station where data were measurement, and it use to save results in *.xls 
#' format. 
#' @param Type: a character specifying a name of probability distribution 
#' function fitted (see \code{\link{selecDIST}}) by \code{\link{fitDISTRI}} function. 
#' @param Intensity: a numeric vector with intensity values for a specific time 
#' duration in different return periods. 
#' @param Parameters: a list with three elements: i) Type of distribution function ii)
#' fitted parameters, and iii) source to call specfic function in the \pkg{lmomco} package.
#' @param M.fit: a character specifying a name of fit method employed on pdf, just three 
#' options are available: L-moments (\emph{Lmoments}), Probability-Weighted Moments (\emph{PWD}), 
#' and Maximum Likelihood (\emph{MLE}). 
#' @param Dura: a character specifying a time duration of the \code{Intensity}, (e.g. 30 min). 
#' This parameter is used to save results. 
#' @param Plot: a number (1) to determine if it will be plotted density curves both empirical 
#' as modeled (\emph{pdf}). If any other number is used graphs will not appear. 
#' @param Resolution: a number to determine the resolution that the plot function will used to save graphs. 
#' It has two options: 300 and 600 ppi. See \code{\link{resoPLOT}}. 
#' @param SAVE: a logical value. TRUE will save \code{Plot}, FALSE will just show it.  
#'
#' @return A data frame with statistics values of goodness of fit tests and its respective p-value, 
#' moreover information criteria are evaluated:
#' 
#' \itemize{
#'  \item \emph{Kolmogorov-Smirnov}: statistic= KS and p-value1
#'  \item \emph{Anderson-Darling}: statistic= AD and p-value2
#'  \item \emph{Cramer-von Mises}: statistic= Omega2 and p-value3
#'  \item \emph{Akaike's Information Criterion}: AIC
#'  \item \emph{Bayesian Information Criterion}: BIC
#'  \item \emph{Akaike's Information Criterion with bais correction}: AICc
#'  \item \emph{Kashyap bayesian Information Criterion}: KIC
#' }
#' 
#' @author David Zamora <dazamoraa@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#' 
#' @references 
#' Hurvich, C. M., & Tsai, C. L. (1989). Regression and time series model 
#' selection in small samples. Biometrika, 76(2), 297-307.
#' 
#' Kashyap, R. L. (1982). Optimal choice of AR and MA parts in autoregressive
#' moving average models. IEEE Transactions on Pattern Analysis and Machine
#' Intelligence, (2), 99-104.
#' 
#' @export
#'
#' @examples
#' 
#' # Meteorology station in the Airport Farfan in Tulua, Colombia.
#' data(Intgum5min)
#' data(Pargumbel)
#' # not plotted
#' test.fit <- goodFIT(Station = "2610516", Type = "Gumbel", Intensity = Intgum5min,
#'                     Parameters = Pargumbel,M.fit = "Lmoments", Dura ="5_min", Plot = 0)
#'  
goodFIT <- function(Station =..., Type =..., Intensity =..., Parameters =...,
                    M.fit =..., Dura =.., Plot =..., Resolution = 300, SAVE = FALSE){
  
  Type <- tolower(Type)
  Plot <- as.character(Plot)
  
  if(Type == "exponential"){
    # ----Exponencial----
    #----Kolmogorov-Smirnov----
    KST <- ks.test(Intensity, "pexp1", xi = Parameters$para[1], alpha = Parameters$para[2], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT <- goftest::ad.test(Intensity, "pexp1", xi = Parameters$para[1], alpha = Parameters$para[2])
    #--------Cramer-von Mises----
    CWT <- goftest::cvm.test(Intensity, "pexp1", xi = Parameters$para[1], alpha = Parameters$para[2])
    
  }else if(Type == "gamma"){
    # ----Gamma----
    #----Kolmogorov-Smirnov----
    KST <- ks.test(Intensity, "pgam1", alpha = Parameters$para[1], beta = Parameters$para[2], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT <- goftest::ad.test(Intensity, "pgam1", alpha = Parameters$para[1], beta = Parameters$para[2])
    #--------Cramer-von Mises----
    CWT <- goftest::cvm.test(Intensity, "pgam1", alpha = Parameters$para[1], beta = Parameters$para[2])
    
  }else if(Type == "gev"){
    # ----GEV----
    #----Kolmogorov-Smirnov----
    KST <- ks.test(Intensity, "pgev1", xi = Parameters$para[1], alpha = Parameters$para[2], kappa = Parameters$para[3], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT <- goftest::ad.test(Intensity, "pgev1", xi = Parameters$para[1], alpha = Parameters$para[2], kappa = Parameters$para[3])
    #--------Cramer-von Mises----
    CWT <- goftest::cvm.test(Intensity, "pgev1", xi = Parameters$para[1], alpha = Parameters$para[2], kappa = Parameters$para[3])
    
  }else if(Type == "gumbel"){
    # ----Gumbel----
    #----Kolmogorov-Smirnov----
    KST <- ks.test(Intensity, "pgum1", xi = Parameters$para[1], alpha = Parameters$para[2], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT <- goftest::ad.test(Intensity, "pgum1", xi = Parameters$para[1], alpha = Parameters$para[2])
    #----Cramer-von Mises----
    CWT <- goftest::cvm.test(Intensity, "pgum1", xi = Parameters$para[1], alpha = Parameters$para[2])
    
  }else if(Type == "log.normal3"){
    # ----Log Normal 3----
    #----Kolmogorov-Smirnov----
    KST <- ks.test(Intensity, "pln31", zeta = Parameters$para[1], mulog = Parameters$para[2], sigmalog = Parameters$para[3], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT <- goftest::ad.test(Intensity, "pln31", zeta = Parameters$para[1], mulog = Parameters$para[2], sigmalog = Parameters$para[3])
    #----Cramer-von Mises----
    CWT <- goftest::cvm.test(Intensity, "pln31", zeta = Parameters$para[1], mulog = Parameters$para[2], sigmalog = Parameters$para[3])
    
  }else if(Type == "normal"){
    # ----Normal----
    #----Kolmogorov-Smirnov----
    KST <- ks.test(Intensity, "pnor1", mu = Parameters$para[1], sigma = Parameters$para[2], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT <- goftest::ad.test(Intensity, "pnor1", mu = Parameters$para[1], sigma = Parameters$para[2])
    #----Cramer-von Mises----
    CWT <- goftest::cvm.test(Intensity, "pnor1", mu = Parameters$para[1], sigma = Parameters$para[2])

  }else if(Type == "pearson3" | Type == "log.pearson3"){
    # ----Person III----
    if(Type == "log.pearson3"){
      Intensity<-log10(Intensity)
    }
    #----Kolmogorov-Smirnov----
    KST <- ks.test(Intensity, "ppe31", mu = Parameters$para[1], sigma = Parameters$para[2], gamma = Parameters$para[3], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT <- goftest::ad.test(Intensity, "ppe31", mu = Parameters$para[1], sigma = Parameters$para[2], gamma = Parameters$para[3])
    #----Cramer-von Mises----
    CWT <- goftest::cvm.test(Intensity, "ppe31", mu = Parameters$para[1], sigma = Parameters$para[2], gamma = Parameters$para[3])

  }else if(Type == "wakeby"){
    # ----Wakeby----
    #----Kolmogorov-Smirnov----
    KST <- ks.test(Intensity, "pwak1", xi = Parameters$para[1], alpha = Parameters$para[2], beta = Parameters$para[3], gamma = Parameters$para[4], delta = Parameters$para[5],
                 alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT <- goftest::ad.test(Intensity, "pwak1", xi = Parameters$para[1], alpha = Parameters$para[2], beta = Parameters$para[3], gamma = Parameters$para[4], delta = Parameters$para[5])
    #----Cramer-von Mises----
    CWT <- goftest::cvm.test(Intensity, "pwak1", xi = Parameters$para[1], alpha = Parameters$para[2], beta = Parameters$para[3], gamma = Parameters$para[4], delta = Parameters$para[5])

  }else{ # log.normal2

  }

  # ----Evaua AIC y BIC----
  Test.Info <- InfoCRIT(Intensity, Parameters)

  # ----FIGURE----
  if(grepl("1", Plot)){
    # ---- Grafica el resultado de la PDF observada y la teorica----
    yd <- density(Intensity)
    xd <- seq(min(yd$x,  na.rm = TRUE),  max(yd$x, na.rm = TRUE),  , 100)
    yd2 <- lmomco::dlmomco(xd, Parameters)
    limy <- c(min(yd$y, yd2, na.rm = T), max(yd$y, yd2, na.rm = T))

    if (SAVE) {
      if(file.exists(paste(".", "FIGURES", Station, M.fit, Type, sep = "/"))){
        path.fig <- paste(".", "FIGURES", Station, M.fit, Type, sep = "/")
      }else{
        dir.create(paste(".", "FIGURES", Station, M.fit, Type, sep = "/"), recursive = TRUE)
        path.fig <- paste(".", "FIGURES", Station, M.fit, Type, sep = "/")
      }

      custom <- resoPLOT(reso = Resolution)

      png(filename = paste(path.fig, "/Den", "_", Type, "_", M.fit, "_", Dura, ".png", sep = ""),
          width = custom[2], height = custom[4], pointsize = 10,  res = custom[1], bg = "transparent")
    }
    par(mar = c(4.1, 3.5, 2.2, 2) + 0.1)
    par(mgp = c(2.5, 0.2, 0))
    plot(yd,  main = paste("Empirical density vs PDF ", Type, "-", M.fit, "\n", Station, sep = ""),
         ylim = limy, xaxt = "n", yaxt = "n", bty = "n", col = "red", lty = 1, lwd = 1,
         ylab = "Density", cex.lab = 1, cex.main = 0.9)
    lines(xd,  yd2,  lty = 2,  col = "blue",  lwd = 1.5)
    box(col = "gray56", lwd = 0.7)
    magicaxis::magaxis(1, tcl=-0.5, cex.axis=0.9, col = "gray56")
    magicaxis::magaxis(2, las=2, cex.axis=0.9, col = "gray56")
    mtext(side = 1, text = paste("Intensity [mm/h] - ", Dura, sep = ""), line = 1.6, cex = 1.1)
    legend("topright",  legend = c("Empirical",  "Modeling"),  col = c("red",  "blue"),
           lty = c(1, 2),  lwd = c(1,  1.5),  bty = "n", cex=0.8)
    if (SAVE) {
      dev.off()
    }
  }

  # -----RESULTS----
  return(data.frame(KS = KST$statistic, p.value1 = KST$p.value,
                    AD = ADT$statistic, p.value2 = ADT$p.value,
                    Omega2 = CWT$statistic, p.value3 = CWT$p.value,
                    AIC = Test.Info[1], BIC = Test.Info[2], AICc = Test.Info[3], KIC = Test.Info[4])
  )
}
