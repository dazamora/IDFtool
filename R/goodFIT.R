#' Title
#'
#' @param Station 
#' @param Type 
#' @param Intensity 
#' @param Parameters 
#' @param M.fit 
#' @param Dura 
#' @param Plot 
#' @param Resolution 
#' @param SAVE 
#'
#' @return
#' @export
#'
#' @examples
goodFIT <- function(Station =..., Type =..., Intensity =..., Parameters =...,
         M.fit =..., Dura =.., Plot =..., Resolution = 300, SAVE = FALSE){

  Type <- tolower(Type)
  Plot <- as.character(Plot)

  if(Type == "exponencial"){
    # ----Exponencial----
    #----Kolmogorov-Smirnov----
    KST<-ks.test(Intensity, "pexp1", xi = Parameters$para[1], alpha = Parameters$para[2], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT<-goftest::ad.test(Intensity, "pexp1", xi = Parameters$para[1], alpha = Parameters$para[2])
    #--------Cramer-von Mises----
    CWT<-goftest::cvm.test(Intensity, "pexp1", xi = Parameters$para[1], alpha = Parameters$para[2])

  }else if(Type == "gamma"){
    # ----Gamma----
    #----Kolmogorov-Smirnov----
    KST<-ks.test(Intensity, "pgam1", alpha = Parameters$para[1], beta = Parameters$para[2], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT<-goftest::ad.test(Intensity, "pgam1", alpha = Parameters$para[1], beta = Parameters$para[2])
    #--------Cramer-von Mises----
    CWT<-goftest::cvm.test(Intensity, "pgam1", alpha = Parameters$para[1], beta = Parameters$para[2])

  }else if(Type == "gev"){
    # ----GEV----
    #----Kolmogorov-Smirnov----
    KST<-ks.test(Intensity, "pgev1", xi = Parameters$para[1], alpha = Parameters$para[2], kappa = Parameters$para[3], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT<-goftest::ad.test(Intensity, "pgev1", xi = Parameters$para[1], alpha = Parameters$para[2], kappa = Parameters$para[3])
    #--------Cramer-von Mises----
    CWT<-goftest::cvm.test(Intensity, "pgev1", xi = Parameters$para[1], alpha = Parameters$para[2], kappa = Parameters$para[3])

  }else if(Type == "gumbel"){
    # ----Gumbel----
    #----Kolmogorov-Smirnov----
    KST<-ks.test(Intensity, "pgum1", xi = Parameters$para[1], alpha = Parameters$para[2], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT<-goftest::ad.test(Intensity, "pgum1", xi = Parameters$para[1], alpha = Parameters$para[2])
    #----Cramer-von Mises----
    CWT<-goftest::cvm.test(Intensity, "pgum1", xi = Parameters$para[1], alpha = Parameters$para[2])

  }else if(Type == "log.normal3"){
    # ----Log Normal 3----
    #----Kolmogorov-Smirnov----
    KST<-ks.test(Intensity, "pln31", zeta = Parameters$para[1], mulog = Parameters$para[2], sigmalog = Parameters$para[3], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT<-goftest::ad.test(Intensity, "pln31", zeta = Parameters$para[1], mulog = Parameters$para[2], sigmalog = Parameters$para[3])
    #----Cramer-von Mises----
    CWT<-goftest::cvm.test(Intensity, "pln31", zeta = Parameters$para[1], mulog = Parameters$para[2], sigmalog = Parameters$para[3])

  }else if(Type == "normal"){
    # ----Normal----
    #----Kolmogorov-Smirnov----
    KST<-ks.test(Intensity, "pnor1", mu = Parameters$para[1], sigma = Parameters$para[2], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT<-goftest::ad.test(Intensity, "pnor1", mu = Parameters$para[1], sigma = Parameters$para[2])
    #----Cramer-von Mises----
    CWT<-goftest::cvm.test(Intensity, "pnor1", mu = Parameters$para[1], sigma = Parameters$para[2])

  }else if(Type == "pearson3" | Type == "log.pearson3"){
    # ----Person III----
    if(Type == "log.pearson3"){
      Intensity<-log10(Intensity)
    }
    #----Kolmogorov-Smirnov----
    KST<-ks.test(Intensity, "ppe31", mu = Parameters$para[1], sigma = Parameters$para[2], gamma = Parameters$para[3], alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT<-goftest::ad.test(Intensity, "ppe31", mu = Parameters$para[1], sigma = Parameters$para[2], gamma = Parameters$para[3])
    #----Cramer-von Mises----
    CWT<-cgoftest::vm.test(Intensity, "ppe31", mu = Parameters$para[1], sigma = Parameters$para[2], gamma = Parameters$para[3])

  }else if(Type == "wakeby"){
    # ----Wakeby----
    #----Kolmogorov-Smirnov----
    KST<-ks.test(Intensity, "pwak1", xi = Parameters$para[1], alpha = Parameters$para[2], beta = Parameters$para[3], gamma = Parameters$para[4], delta = Parameters$para[5],
                 alternative  =  "two.sided")
    #----Anderson-Darling----
    ADT<-goftest::ad.test(Intensity, "pwak1", xi = Parameters$para[1], alpha = Parameters$para[2], beta = Parameters$para[3], gamma = Parameters$para[4], delta = Parameters$para[5])
    #----Cramer-von Mises----
    CWT<-goftest::cvm.test(Intensity, "pwak1", xi = Parameters$para[1], alpha = Parameters$para[2], beta = Parameters$para[3], gamma = Parameters$para[4], delta = Parameters$para[5])

  }else{ # log.normal2

  }

  # ----Evaua AIC y BIC----
  Test.Info<-InfoCRIT(Intensity, Parameters)

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
    plot(yd,  main = paste("Densidad Empirica vs FDP ", Type, "-", M.fit, "\n", Station, sep = ""),
         ylim = limy, xaxt = "n", yaxt = "n", bty = "n", col = "red", lty = 1, lwd = 1,
         ylab = "Densidad", cex.lab = 1, cex.main = 0.9)
    lines(xd,  yd2,  lty = 2,  col = "blue",  lwd = 1.5)
    box(col = "gray56", lwd = 0.7)
    magicaxis::magaxis(1, tcl=-0.5, cex.axis=0.9, col = "gray56")
    magicaxis::magaxis(2, las=2, cex.axis=0.9, col = "gray56")
    mtext(side = 1, text = paste("Intensity [mm/h] - ", Dura, sep = ""), line = 1.6, cex = 1.1)
    legend("topright",  legend = c("Empirico",  "Modelado"),  col = c("red",  "blue"),
           lty = c(1, 2),  lwd = c(1,  1.5),  bty = "n", cex=0.8)
    if (SAVE) {
      dev.off()
    }
  }

  # -----RESULTS----
  return(data.frame(KS = KST$statistic, p.value1 = KST$p.value,
                    AD = ADT$statistic, p.value2 = ADT$p.value,
                    Omega2 = CWT$statistic, p.value3 = CWT$p.value,
                    AIC = Test.Info[1], BIC = Test.Info[2])
  )
}
