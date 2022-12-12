#' @name
#' IDFCurve
#' 
#' @title 
#' Compute Intensity-Duration-Frequency curve
#' 
#' @description 
#' An Intensity-Duration-Frequency curve (IDF Curve) is a graphical representation 
#' of the probability that a given average rainfall intensity will occur. This function allows to fit different 
#' probability distribution functions (see \code{\link{selecDIST}}) by means of four fit methods (see \code{\link{fitDISTRI}})
#' to determine intensity [mm/h] for different return periods and per specific time durations. Finally, it computes equations parameters 
#' of the IDF curves (see \code{\link{regIDF}})
#' 
#' 
#' @param Data a numeric matrix with years in the first column and intensity values for each \code{Duration} time in the other columns.
#' @param Station a string with a name to identify source of \code{Data}.
#' @param Duration a logical value or numeric vector. If it is TRUE the durations (in minutes) 
#' by default will be used: 5, 10 ,15, 20, 30, 60, 120 y 360. In case of a numeric vector the durations must be in minutes.
#' @param Periods a logical value or numeric vector.
#' @param Type a character specifying the name of the distribution function that will 
#' be employed: exponencial, gamma, gev, gumbel, log.normal3, normal, pearson, log.pearson3 and 
#' wakeby (see \code{\link{selecDIST}}).
#' @param M.fit a character specifying a name of fit method employed on pdf, just three 
#' options are available: L-moments (\emph{Lmoments}), Probability-Weighted Moments (\emph{PWD}), 
#' Maximum Likelihood (\code{\link{MLEZ}}) and Moments (\emph{MME}) (see \code{\link{MME_DIST}}).
#' @param Plot it is a number of one to four digits. a number (1) to plot density curves 
#' both empirical as modeled (\emph{pdf}).(2) to plot if it will be 
#' plotted curves between return \code{Periods} and intensity computed by \emph{pdf} fitted. 
#' Or use (12) to get both graphs. (3) to plot IDF curves for all return periods: (\code{Durations} versus \code{Intensity})
#' for all return \code{Periods}. (4) to plot IDF curve each for return
#' \code{Periods} with its confidence and prediction intervals. Or use (34) to get both graphs. If other is used the graphs will not appear.
#' @param Strategy a numeric vector used to identify Strategies to compute IDF curves with different data sets: 1 just data from HIDFUN tool, 
#' 2 just data from IDEAM and 3 used this data sets. These strategies were created just for Colombian pluviograph data. For your data
#' use strategy number 1.
#' @param logaxe a character to plot axis in log scale: x, y or both (xy). In other case used "".
#' @param CI a logical value specifying whether confidence and prediction intervals will be computed in IDF curves.
#' @param CIpdf a logical value specifying whether confidence of pdf will be computed.
#' @param iter an integer representing number of resamples to conduct when 
#' confidence interval is computed (see \code{\link{bootstrapCI}}). Use it only if 
#' CI is equal to TRUE.
#' @param goodtest a logical value specifying whether goodness-fit tests should be 
#' cumputed to \emph{pdf} fitted by means of \code{\link{goodFIT}} function.
#' @param Resolution a number to determine the resolution that the plot function will used to save graphs. 
#' It has two options: 300 and 600 ppi. See \code{\link{resoPLOT}}.
#' @param SAVE a logical value. TRUE will save \code{Plot}, FALSE will just show \code{Plot}.
#' @param name a logical value. TRUE will use a default names to identify strategies: (1) "HIDFUN", (2) "IDEAM", (3) "AMBOS". In other case FALSE
#' allows: i) selected years of data sets, and ii) to insert durations to do IDF curves (in minutes)
#'
#' @return A list of:
#'
#'  \itemize{
#'    \item \code{Intesities} a numeric matrix of intensities values per each return \code{Periods} computed by \emph{pdf} fitted.
#'    \item \code{Models} a list with results of the function \code{\link{regIDF}}.
#'    \item \code{Test.fit} a list with results of the function \code{\link{goodFIT}}.
#'    \item \code{Distribution} a list with results of the function \code{\link{fitDISTRI}}.
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
#' Plot = 1234, Strategy = 1, logaxe = "", CI = FALSE, CIpdf = TRUE, iter = 50,
#' goodtest = FALSE, Resolution = 300, SAVE = FALSE, name = TRUE)
#' 
IDFCurve <- function(Data, Station = '2610516', Duration = FALSE,
                   Periods = FALSE, Type = "gumbel", M.fit = "lmoments",
                   Plot = 1234, Strategy = 1:3, logaxe = "", CI = FALSE, CIpdf = TRUE, iter = 500,
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
  if(dim(input)[2]==2){
    Int.total <- as.matrix(Int.total)
    colnames(Int.total) <- colnames(input)[2]
  }
  id.info <- which(is.na(Int.total[ ,1]) == TRUE)
  Plot <- as.character(Plot)
  
  # ----Definir Duraciones----
  if (is.logical(Duration)) {
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
          if(!(sum(id.info) == 0)){
            intensities <- Int.total[-id.info, ]
          } else {
            intensities <- Int.total
          }
        } else {
          
          if(length(Duration) == 1){
            in.dura <- is.element(colnames(Int.total),as.character(duration*60))
          } else {
            in.dura <- is.element(colnames(Int.total),as.character(Duration))
          }
          
          if(!(sum(id.info) == 0)){
            intensities <- Int.total[-id.info,in.dura]
          } else {
            intensities <- Int.total[ ,in.dura]
          }
        }
        durations <- duration
        nom.dura <- paste(duration*60, " min", sep = "")
      } else if (estr == 2) { # Solo datos Ideam
        intensities <- Int.total[id.info,c(3,5:8)]
        durations <- duration[c(3,5:8)]
        nom.dura <- paste(duration[c(3,5:8)]*60, " min", sep = "")
      } else { # Ambos conjunto de datos
        in.dura <- is.element(as.character(duration*60),colnames(Int.total))
        intensities <- as.matrix(Int.total[,in.dura])
        durations <- duration[in.dura]
        nom.dura <- paste(duration[in.dura]*60, " min", sep = "")
      }
    }
    # ----Ajusta la distribucion y calcula intesidades por duracion y periodo de retorno----
    nd <- length(durations)
    nTp <- length(Tp)
    
    #options(warn = 0) # Para mostrar wanings 0
    distri <- list() # Almacena todos los resultados de la funcion fitDISTRI
    idf <- matrix(nrow = nd, ncol = nTp)
    M.test.fit <- matrix(NA, nrow = nd, ncol = 10)
    CI.pdf.lower <- c()
    CI.pdf.upper <- c()
    
    for(i in 1:nd){
      distri[[nom.dura[i]]]<-fitDISTRI(Intensity = intensities[,i], Type = Type, Plot = Plot, M.fit = M.fit,
                                       Periods = Tp, Dura = nom.dura[i], Station = Station, CI = CIpdf, iter = iter,
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
    if(is.null(M.test.fit)){
      M.test.fit <- M.test.fit
    } else {
      colnames(M.test.fit) <- names(distri[[i]]$goodness.fit)
      rownames(M.test.fit) <- nom.dura
    }
    names.periods <- round(lmomco::prob2T(distri[[1]]$Conf.Inter$Conf.Inter$nonexceed.prob),0)
    
    if (CIpdf) {
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
      sheets.content <- list(idf, CI.pdf.lower, CI.pdf.upper, M.test.fit, Output[[name[estr]]]$Coefficients,
                             Output[[name[estr]]]$Prediction.Int, Output[[name[estr]]]$test.fit.reg,
                             Output[[name[estr]]]$Confidence.Int)
      sheets.names <- c("IDF.by.PDF", "CIL-IDF.by.PDF", "CIU-IDF.by.PDF", "goodness.fit",
                        "Coefficients", "Prediction.by.C-IDF", "Performance-IDF.reg", "Conf.Int-IDF.reg")
      openxlsx::write.xlsx(sheets.content, file = paste(path.result,"/", "IDF_", Station, "_", name[estr], ".xlsx",sep=""),
                       sheetName = sheets.names, row.names = TRUE, col.names = TRUE)
      
      if (CI) {
        wb <- openxlsx::loadWorkbook(paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""))  
        openxlsx::writeData(wb, 'Conf.Int-IDF.reg', Output[[name[estr]]]$Confidence.Int)
        openxlsx::saveWorkbook(wb, paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""), overwrite = T)
      }
    }
  }
  
  # ----Salidas de la funcion----
  if (length(Strategy) < 2) {
    return(list(Intensities = idf, 
                Models = Output,
                Test.fit = M.test.fit, 
                Distribution = distri))
  } else {
    return(list(Intensities = Midf, 
                Models = Output,
                Test.fit = Mgoodness, 
                Distribution = distri))
  }
}
