#' Title
#'
#' @param Intensity 
#' @param Type 
#' @param Plot 
#' @param M.fit 
#' @param Periods 
#' @param Dura 
#' @param Station 
#' @param CI 
#' @param iter 
#' @param goodtest 
#' @param Resolution 
#' @param SAVE 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
fitDISTRI <- function(Intensity =..., Type ="Gumbel", Plot = 2, M.fit = "MLE",
         Periods =..., Dura =..., Station =..., CI = FALSE, iter = ..,
         goodtest = FALSE, Resolution = 300, SAVE = FALSE){

  # ----Fix variables----
  M.fit <- tolower(M.fit)
  Type <- tolower(Type)
  Tp <- Periods
  DR <- colnames(Intensity)
  Plot <- as.character(Plot)

  # ----Fit distribution----
  if(Type=="log.pearson3"){
    Intensity<-log10(Intensity)
  }

  distribution<-selecDIST(Type = Type)
  FR<-lmomco::T2prob(Tp)

  if(M.fit=="lmoments"){
    LMOM<-lmomco::lmoms(Intensity)
    Parameters<-lmomco::lmom2par(LMOM, type = distribution)
    INT<-lmomco::par2qua(FR,Parameters)
  }else if(M.fit=="pwd"){
    PMP<-lmomco::pwm(Intensity)
    Parameters<-lmomco::lmom2par(pwm2lmom(PMP),type = distribution)
    INT <-lmomco::par2qua(FR,Parameters)
  }else if(M.fit=="mle"){
    Parameters <- lmomco::mle2par(Intensity,type=distribution)
    INT <- lmomco::par2qua(FR,Parameters)
  }else{
    stop("Error on selecction fit model")
  }

  if(Type=="log.pearson3"){
    INT <- 10^INT
    Intensity <- 10^Intensity
  }
  names(INT) <- as.character(Tp)

  # ----Computed goodness-fit tests-----
  test.fit <- NULL
  if (goodtest) {
    test.fit <- goodFIT(Station = Station, Type = Type, Intensity = Intensity, Parameters = Parameters,
                        M.fit = M.fit,Dura = Dura, Plot = Plot,  Resolution = Resolution ,SAVE = SAVE)
  }
  # ---- Computed confidence interval----
  CI.result <- NULL
  if (CI) {
    Ttick <- c(1.001, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 150)
    FR.plot <- lmomco::T2prob(Ttick)
    CI.result <- bootstrapCI(Intensity = Intensity, Parameters = Parameters, Type = Type,
                             Rsample = iter, Return.P = FR.plot, Conf.Inter = 0.95)
  }
  # ----Plot Frecuency versus Intesity----
  if (grepl("2", Plot)) {

    if (SAVE) {
      if(file.exists(paste(".", "FIGURES", Station, M.fit, Type, sep = "/"))){
        path.fig <- paste(".", "FIGURES", Station, M.fit, Type, sep = "/")
      }else{
        dir.create(paste(".", "FIGURES", Station, M.fit, Type, sep = "/"), recursive = TRUE)
        path.fig <- paste(".", "FIGURES", Station, M.fit, Type, sep = "/")
      }

      custom <- resoPLOT(reso = Resolution)
      png(filename = paste(path.fig, "/IvsF", "_", Type, "_", M.fit, "_", Dura, ".png", sep = ""),
          width = custom[2], height = custom[4], pointsize = 10, res = custom[1], bg = "transparent")
    }

    lim.max <- max(CI.result$Conf.Inter[ ,3], na.rm = T)
    lim.min <- min(CI.result$Conf.Inter[ ,2], na.rm = T)
    lim.vert <- c(lim.min, lim.max)

    par(mar = c(4.1, 3.5, 2.2, 2) + 0.1)
    par(mgp = c(2.2, 0.2, 0))
    plot(Ttick,lmomco::par2qua(FR.plot,Parameters),type = "n",
         main = paste("Intensity vs Frequency ", Type, "-", M.fit, "\n", Station, sep = ""),
         ylim = lim.vert, xaxt = "n", yaxt = "n", bty = "n", xlab = "Retorn periods [year]",
         ylab = "Intensity [mm/h]", cex.lab = 1, cex.main = 0.9,log = "yx")

    abline(v = axTicks(1), h = axTicks(2), col = "gray80", lty = 3)
    polygon(c(Ttick, rev(Ttick)), c(CI.result$Conf.Inter[ ,2], rev(CI.result$Conf.Inter[ ,3])),
            col = scales::alpha("gray75", alpha = 0.2), border = scales::alpha("gray75", alpha = 0.2))
    lines(Ttick,CI.result$Conf.Inter[ ,2], lty = 2, lwd = 0.7, col = "gray68")
    lines(Ttick,CI.result$Conf.Inter[ ,3], lty = 2, lwd = 0.7, col = "gray68")
    lines(Ttick,lmomco::par2qua(FR.plot,Parameters), lty = 4, lwd = 1.1, col = "red")

    NEP.obs <- lmomco::par2cdf(Intensity,Parameters)
    Tp.obs <- lmomco::prob2T(NEP.obs)
    points(Tp.obs, Intensity, pch = 21, bg = "blue", col = "blue",cex = 0.8)
    legend("bottomright", c("CI 95 %", "Fit pdf", "Observed"), col = c("gray68", "red", "blue"),
           pt.bg = c(scales::alpha("gray75", alpha = 0.2), NA ,"blue"), pch = c(22,-1,21),
           lty = c(-1, 4, -1), lwd = 1, bty = "n", border = NULL, cex = 0.7, pt.cex = 1.1,
           title = paste("Duration\n", Dura, sep = ""), title.col = "magenta")

    magicaxis::magaxis(1, labels = FALSE, col = "gray56")
    axis(1,axTicks(1),axTicks(1), tick = FALSE, line = -0.05,cex.axis = 0.9)
    magicaxis::magaxis(2, las = 2, col = "gray56")
    #axis(2,axTicks(2),axTicks(2), lwd.ticks = 0.5, tck = 0.02, las = 2, cex.axis = 0.9)
    box(col = "gray56", lwd = 0.7)

    if (SAVE) {
      dev.off()
    }
  }

  return(list(Parameters = Parameters,
              Int.pdf = INT,
              Conf.Inter = CI.result,
              goodness.fit = test.fit,
              Info.PDF = c(Type,M.fit))
  )
}
