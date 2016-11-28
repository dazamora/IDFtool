
#' IDFCurve
#' 
#' An Intensity-Duration-Frequency curve (IDF Curve) is a graphical representation 
#' of the probability that a given average rainfall intensity will occur. This function allows fit different 
#' probability distribution functions (see \code{\link{selecDIST}}) by means of four fit methods (see \code{\link{fitDISTRI}})
#' to determined intensity [mm/h] for different return periods and per specific time durations. Finally, it compute parameters of the equations 
#' of the IDF curves (see \code{\link{regIDF}})
#' 
#' 
#' @param Data: a numeric matrix with years in the first column and other columns has intensity [mm/h] values by each a the \code{Duration}s time.
#' @param Station: a string with a name to identify source of \code{Data}.
#' @param Duration: a logical value or numeric vector. If it is TRUE will be use the durations (in minutes) 
#' by default: 5, 10 ,15, 20, 30, 60, 120 y 360. In case a numeric vector the durations must be in minutes.
#' @param Periods: a logical value or numeric vector.
#' @param Type: a character specifying the name of distribution function that it will 
#' be employed: exponencial, gamma, gev, gumbel, log.normal3, normal, pearson, log.pearson3 and 
#' wakeby (see \code{\link{selecDIST}}).
#' @param M.fit: a character specifying a name of fit method employed on pdf, just three 
#' options are available: L-moments (\emph{Lmoments}), Probability-Weighted Moments (\emph{PWD}), 
#' Maximum Likelihood (\code{\link{MLEZ}}) and Moments (\emph{MME}) (see \code{\link{MME_DIST}}).
#' @param Plot: it is a number of one to four digits. a number (1) to determine if it will be plotted density curves 
#' both empirical as modeled (\emph{pdf}). a number (2) to determine if it will be 
#' plotted curves between return \code{Periods} and intensity computed by \emph{pdf} fitted. 
#' Or use both numbers to get these graphs. If you use other number the graphs 
#' will not appear. he number three (3) determined if it will plot IDF curves (\code{Durations} versus \code{Intensity})
#' for all return \code{Periods}. The number four (4) determined if it will plot IDF curve each for return
#' period with its confidence and prediction intervals. Or use both numbers to get these graphs. If you use other number the graphs will not appear.
#' @param Strategy: a numeric vector used to identify Strategies to compute IDF curves with different datas sets: 1 just data from Ideam, 
#' 2 just data from HIDFUN tool and 3 used this data sets.
#' @param logaxe: a character to plot axis in log scale: x, y or both (xy). In other case used "".
#' @param CI: a logical value specifying whether confidence and prediction intervals will be computed.
#' @param iter: an integer representing number of resamples to conduct when 
#' confidence interval will be computed (see \code{\link{bootstrapCI}}). Use it only if 
#' CI is equal to TRUE.
#' @param goodtest: a logical value specifying whether goodness-fit tests should be 
#' cumputed to \emph{pdf} fitted by means of \code{\link{goodfit}} function.
#' @param Resolution: a number to determine resolution that the plot function used to save graphs. 
#' It can have two options: 300 and 600 ppi. See \code{\link{resoPLOT}}.
#' @param SAVE: a logical value. TRUE will save \code{Plot} but if is FALSE just show \code{Plot}.
#' @param name: a logical value. TRUE will use a default names to identify strategies: (1) "HIDFUN", (2) "IDEAM", (3) "AMBOS". In other case FALSE
#' allows: i) selected years of data sets, ii) insert durations to do IDF curves (in minutes)
#'
#' @return A list of:
#'
#'  \itemize{
#'    \item \code{Intesities} a numeric matrix of intensities values per each return \code{Periods} compute by \emph{pdf} fitted.
#'    \item \code{Models} a list with results of the function \code{\link{regIDF}}.
#'    \item \code{Test.fit} a list with results of the function \code{\link{goodFIT}}.
#'    \item \code{Distribution} a list with results of the function \code{\link{fitDISRTI}}.
#'  }
#' @author David Zamora <dazamoraa@unal.edu.co>
#' Albeiro Figueroa <cafigueroao@unal.edu.co> 
#' Water Resources Engineering Research Group - GIREH
#'
#' @export
#'
#' @examples
#' 
#' # Meteorology station in the Farfan Airport in Tulua, Colombia.
#' data(inten)
#' Test.idftool <- IDFCurve(Data = inten, Station='2610516', Duration = FALSE,
#' Periods = FALSE, Type = "gumbel", M.fit = "lmoments",
#' Plot = 1234, Strategy = 1, logaxe = "", CI = FALSE, iter = 100,
#' goodtest = FALSE, Resolution = 300, SAVE = FALSE, name = TRUE)
#' 
IDFCurve<-function(Data =..., Station='2610516', Duration = FALSE,
                   Periods = FALSE, Type = "gumbel", M.fit = "lmoments",
                   Plot = 1234, Strategy = 1:3, logaxe = "", CI = FALSE, iter = 500,
                   goodtest = FALSE, Resolution = 300, SAVE = FALSE, name = TRUE){
  
  # ----INPUTS----
  if (name) {
    name <- c("HIDFUN","IDEAM","AMBOS")
    option.alt <- 0
  } else {
    name <- "Select-data"
    option.alt <- 4
  }
  # ----OUTPUTS----
  Output<-list()
  
  # ----Load data----
  if (is.list(Data)) {
    input <- Data[[Station]]  
  } else {
    input <- Data
  }
  
  Int.total <- input[ ,2:dim(input)[2]]
  id.info <- which(is.na(Int.total[ ,1]) == TRUE)
  Plot <- as.character(Plot)
  
  # ----Definir Duraciones----
  if (length(Duration) == 1) {
    duration <- c(5/60, 10/60, 15/60, 20/60, 30/60, 1, 2, 6)   # durations in hr for idf analysis
  } else {
    duration <- Duration/60
  }
  # ----Definir Periods----
  if (length(Periods) == 1) {
    Tp<-c(2, 3, 5, 10, 25, 50, 100)
  } else {
    Tp<-Periods                            # return periods in years for plotting idf curves
  }
  
  # ----Corre los diferentes escenarios----
  for (estr in Strategy) {
    
    if (length(Strategy) < 2) {
      print("Just compute a strategy")
    } else {
      Mgoodness <- list()
      Midf <- list()
    }
    
    if (option.alt == 4) {
      print(paste("Min. Year: ", min(input[ ,1]), " - ", "Max. Year: ", max(input[ ,1]), sep = ""))
      cat("\n Select the range of years")
      years <- scan(what = numeric(), nmax=2, quiet = TRUE)
      print("Time duration in hours")
      print(duration)
      cat("\n Are you going to work with durations defined already? (Yes or No)")
      ans <- scan(what = character(), nmax = 1, quiet = TRUE)
      
      if (tolower(ans) == "yes"){
        print("Right")
      } else if (tolower(ans) == "no") {
        print("Select some of the next durations (minutes)")
        print(colnames(input[,-1])[c(3,5:8)])
        duration <- scan(what = numeric(), n = 5, quiet = TRUE)
      } else{
        stop("Incorrect selection")
      }
      idcol <- which(is.element(colnames(input[ ,-1]),as.character(duration)) == TRUE)
      idrow <- which(input[ ,1] >= years[1] & input[ ,1] <= years[2]) 
      intensities <- input[idrow,-1][ ,idcol]
      durations <- duration
      nom.dura <- paste(colnames(intensities), " min", sep = "")
    } else {
      if (estr == 1) { # Solo datos digitalizados
        if (length(Duration) == 1){
          intensities <- Int.total[-id.info, ]
        } else {
          in.dura <- is.element(colnames(Int.total),as.character(Duration))
          intensities <- Int.total[-id.info,in.dura]
        }
        durations <- duration
        nom.dura <- paste(duration*60, " min", sep = "")
      } else if (estr == 2) { # Solo datos Ideam
        intensities <- Int.total[id.info,c(3,5:8)]
        durations <- duration[c(3,5:8)]
        nom.dura <- paste(duration[c(3,5:8)]*60, " min", sep = "")
      } else { # Ambos conjunto de datos
        intensities <- Int.total[,c(3,5:8)]
        durations <- duration[c(3,5:8)]
        nom.dura <- paste(duration[c(3,5:8)]*60, " min", sep = "")
      }
    }
    # ----Ajusta la distribucion y calcula intesidades por duracion y periodo de retorno----
    nd <- length(durations)
    nTp <- length(Tp)
    
    #options(warn = 0) # Para mostrar wanings 0
    distri <- list() # Almacena todos los resultados de la funcion fitDISTRI
    idf <- matrix(nrow = nd, ncol = nTp)
    M.test.fit <- matrix(NA, nrow = nd, ncol = 8)
    CI.pdf.lower <- c()
    CI.pdf.upper <- c()
    
    for(i in 1:nd){
      distri[[nom.dura[i]]]<-fitDISTRI(Intensity = intensities[,i], Type = Type, Plot = Plot, M.fit = M.fit,
                                       Periods = Tp, Dura = nom.dura[i], Station = Station, CI = CI, iter = iter,
                                       goodtest = goodtest, Resolution = Resolution, SAVE = SAVE)
      idf[i, ] <- distri[[i]]$Int.pdf
      if(!is.null(distri[[i]]$goodness.fit)) {
        M.test.fit[i, ]<- unlist(distri[[i]]$goodness.fit)
      } else {
        M.test.fit <- NULL
      }
      if(!is.null(distri[[i]]$Conf.Inter)) {
        CI.pdf.lower <- cbind(CI.pdf.lower, distri[[i]]$Conf.Inter$Conf.Inter$lower.lim)
        CI.pdf.upper <- cbind(CI.pdf.upper, distri[[i]]$Conf.Inter$Conf.Inter$upper.lim)
      } else {
        CI.pdf.lower <- NULL
        CI.pdf.upper <- NULL
      }
    }
    
    colnames(idf) <- as.character(Tp)
    rownames(idf) <- nom.dura
    colnames(M.test.fit) <- names(distri[[i]]$goodness.fit)
    rownames(M.test.fit) <- nom.dura
    names.periods <- round(lmomco::prob2T(distri[[1]]$Conf.Inter$Conf.Inter$nonexceed.prob),0)
    
    if (CI) {
      colnames(CI.pdf.lower) <- nom.dura
      rownames(CI.pdf.lower) <- as.character(names.periods)
      colnames(CI.pdf.upper) <- nom.dura
      rownames(CI.pdf.upper) <- as.character(names.periods)
    }
    # ----Compute idf equations per each time durations----
    if(option.alt == 4){
      durations <- durations/60
    }
    
    Output[[name[estr]]] <- regIDF(Intensity = idf, Periods = Tp, Durations= durations*60, logaxe = logaxe,
                                   Plot = Plot, Resolution = Resolution, SAVE = SAVE, Strategy = estr, Intervals = CI,
                                   M.fit = M.fit, Type = Type, name = name, Station = Station)
    
    if (length(Strategy) > 2) {
      Mgoodness[[name[estr]]] <- M.test.fit
      Midf[[name[estr]]] <- idf
    }
    
    # ----Save results in worksheets of Excel----
    if (SAVE) {
      if (file.exists(paste(".", "RESULTS", Station, sep = "/"))) {
        path.result<-paste(".", "RESULTS", Station, sep = "/")
      } else {
        dir.create(paste(".", "RESULTS", Station, sep = "/"), recursive = TRUE)
        path.result <- paste(".", "RESULTS", Station, sep = "/")
      }
      
      xlsx::write.xlsx(idf, file = paste(path.result,"/", "IDF_", Station, "_", name[estr], ".xlsx",sep=""),
                       sheetName = "IDF.by.PDF", row.names = TRUE, col.names = TRUE)
      xlsx::write.xlsx(CI.pdf.lower, file = paste(path.result,"/", "IDF_", Station, "_", name[estr], ".xlsx",sep=""),
                       sheetName = "CIL-IDF.by.PDF", append = TRUE, row.names = TRUE, col.names = TRUE)
      xlsx::write.xlsx(CI.pdf.upper, file = paste(path.result,"/", "IDF_", Station, "_", name[estr], ".xlsx",sep=""),
                       sheetName = "CIU-IDF.by.PDF", append = TRUE, row.names = TRUE, col.names = TRUE)
      xlsx::write.xlsx(M.test.fit, file = paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
                       sheetName = "goodness.fit", append = TRUE, row.names = TRUE, col.names = TRUE)
      xlsx::write.xlsx(Output[[name[estr]]]$Coefficients, file = paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
                       sheetName = "Coefficients", append = TRUE, row.names = TRUE, col.names = TRUE)
      xlsx::write.xlsx(Output[[name[estr]]]$Predict,file = paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
                       sheetName = "Prediction.by.C-IDF", append = TRUE, row.names = TRUE, col.names = TRUE)
      xlsx::write.xlsx(Output[[name[estr]]]$test.fit.reg, file =  paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
                       sheetName="Performance-IDF.reg", append = TRUE, row.names = TRUE, col.names = TRUE)
      xlsx::write.xlsx(Output[[name[estr]]]$Confidence.Int, file = paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
                       sheetName="Conf.Int-IDF.reg", append = TRUE, row.names = TRUE, col.names = TRUE)
    }
  }
  
  # ----Salidas de la funcion----
  if (length(Strategy) < 2) {
    return(list(Intesities = idf, 
                Models = Output,
                Test.fit = M.test.fit, 
                Distribution = distri)
    )
  } else {
    return(list(Intesities = Midf, 
                Models = Output,
                Test.fit = Mgoodness, 
                Distribution = distri)
    )
  }
}
