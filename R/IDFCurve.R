
#' IDFCurve
#'
#'
#'
#' @param Data: 1
#' @param Station: 11
#' @param Duration: 11
#' @param Periods: 11
#' @param Type: 11
#' @param M.fit: 11
#' @param Plot: 11
#' @param Strategy: 11
#' @param logaxe: 11
#' @param CI: 11
#' @param iter: 11
#' @param goodtest: 11
#' @param Resolution: 11
#' @param SAVE: 11
#' @param name: 12
#'
#' @return A list of:
#'
#'  \itemize{
#'    \item \code{Parameters} a list
#'    \item \code{Int.pdf}
#'    \item \code{Conf.Inter}
#'    \item \code{goodness.fit}
#'    \item \code{Info.PDF}
#'  }
#' @author David Zamora <dazamoraa@unal.edu.co>
#' Water Resources Engineering Research Group - GIREH
#'
#' @export
#'
#' @examples
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
    duration <- c(Duration[Duration < 1]/60, Duration[Duration >= 1])
  }
  # ----Definir Periods----
  if (length(Periods) == 1) {
    Tp<-c(2, 3, 5, 10, 25, 50, 100)
  }else{
    Tp<-Periods                            # return periods in years for plotting idf curves
  }
  
  # ----Corre los diferentes escenarios----
  for (estr in Strategy) {
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
      nom.dura <- colnames(intensities)
    } else {
      if (estr == 1) { # Solo datos digitalizados
        intensities <- Int.total[-id.info, ]
        durations <- duration
        nom.dura <- paste(c(5, 10, 15, 20, 30, 60, 120, 360), " min", sep = "")
      } else if (estr == 2) { # Solo datos Ideam
        intensities <- Int.total[id.info,c(3,5:8)]
        durations <- duration[c(3,5:8)]
        nom.dura <- paste(c(15, 30, 60, 120, 360), " min", sep = "")
      } else { # Ambos conjunto de datos
        intensities <- Int.total[,c(3,5:8)]
        durations <- duration[c(3,5:8)]
        nom.dura <- paste(c(15, 30, 60, 120, 360), " min", sep = "")
      }
    }
    # ----Ajusta la distribucion y calcula intesidades por duracion y periodo de retorno----
    nd <- length(durations)
    nTp <- length(Tp)
    
    options(warn = 0) # Para mostrar wanings 0
    distri <- list() # Almacena todos los resultados de la funcion fitDISTRI
    idf <- matrix(nrow = nd, ncol = nTp)
    M.test.fit <- matrix(NA, nrow = nd, ncol = 8)
    CI.pdf.lower <- c()
    CI.pdf.upper <- c()
    
    for(i in 1:nd){
      distri[[nom.dura[i]]]<-fitDISTRI(Intensity = intensities[,i], Type = Type, Plot = Plot, M.fit = M.fit,
                                       Periods = Tp, Dura = paste(as.character(durations[i]), " min.", sep = ""),
                                       Station = Station, CI = CI, iter = iter,
                                       goodtest = goodtest, Resolution = Resolution, SAVE = SAVE)
      idf[i, ] <- distri[[i]]$Int.pdf
      M.test.fit[i, ]<- distri[[i]]$goodness.fit
      CI.pdf.lower <- cbind(CI.pdf.lower, distri[[i]]$Conf.Inter$lower.lim)
      CI.pdf.upper <- cbind(CI.pdf.upper, distri[[i]]$Conf.Inter$upper.lim)
    }
    
    colnames(idf) <- as.character(Tp)
    rownames(idf) <- nom.dura
    colnames(M.test.fit) <- names(distri[[i]]$goodness.fit)
    rownames(M.test.fit) <- nom.dura
    names.periods <- round(lmomco::prob2T(distri[[1]]$Conf.Inter$nonexceed.prob),0)
    colnames(CI.pdf.lower) <- nom.dura
    rownames(CI.pdf.lower) <- as.character(names.periods)
    colnames(CI.pdf.upper) <- nom.dura
    rownames(CI.pdf.upper) <- as.character(names.periods)
    
    # ----Compute idf equations per each time durations----
    Output[[name[estr]]] <- regIDF(Intensity = idf, Periods = Tp, Durations= durations, logaxe = logaxe,
                                   Plot = Plot, Resolution = Resolution, SAVE = SAVE, Strategy = estr,
                                   M.fit = M.fit, Type = Type, name = name, Station = Station)
    
    # ----Save results in worksheets of Excel----
    if (SAVE) {
      if (file.exists(paste(".", "RESULTS", Station, sep = "/"))) {
        path.result<-paste(".", "RESULTS", Station, sep = "/")
      } else {
        dir.create(paste(".", "RESULTS", Station, sep = "/"), recursive = TRUE)
        path.result <- paste(".", "RESULTS", Station, sep = "/")
      }
      
      xlxs::write.xlsx(idf, file = paste(path.result,"/", "IDF_", Station, name[estr], ".xlsx",sep=""),
                       sheetName = "IDF.by.PDF", row.names = TRUE, col.names = TRUE)
      xlxs::write.xlsx(CI.pdf.lower, file = paste(path.result,"/", "IDF_", Station, name[estr], ".xlsx",sep=""),
                       sheetName = "CIL-IDF.by.PDF", append = TRUE, row.names = TRUE, col.names = TRUE)
      xlxs::write.xlsx(CI.pdf.upper, file = paste(path.result,"/", "IDF_", Station, name[estr], ".xlsx",sep=""),
                       sheetName = "CIU-IDF.by.PDF", append = TRUE, row.names = TRUE, col.names = TRUE)
      xlsx::write.xlsx(M.test.fit, file = paste(path.result, "/", "IDF_", Station, name[estr], ".xlsx", sep = ""),
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
  return(list(Intesidades = idf, 
              Modelos = Output,
              Ajuste = M.test.fit, 
              Distibucion = distri)
         )
  
}
