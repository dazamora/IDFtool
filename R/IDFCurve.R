

#' IDFCurve
#'
#'
#'
#' @param Data:
#' @param Station:
#' @param Duration:
#' @param Periods:
#' @param Type:
#' @param M.fit:
#' @param Plot:
#' @param Strategy:
#' @param logaxe:
#' @param CI:
#' @param iter:
#' @param goodtest:
#' @param Resolution:
#' @param SAVE:
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
#'
IDFCurve<-function(Data =..., Station='2610516', Duration = FALSE,
              Periods = FALSE, Type = "gumbel", M.fit = "lmoments",
              Plot = 1234, Strategy = 1:3, logaxe="", CI = FALSE, iter = 500,
              goodtest = FALSE, Resolution = 300, SAVE = FALSE){

  # ----INPUTS----
  name<-c("HIDFUN","IDEAM","AMBOS")

  # ----OUTPUTS----
  Output<-list()

  # ----Cargar datos----
  input <-LISTA.DURACION[[Station]]
  Int.total <-input[ ,2:dim(input)[2]]
  id.info<-which(is.na(Int.total[,1])==TRUE)
  Plot<-as.character(Plot)

  # ----Definir Durationes----
  if(length(Duration)==1){
    duration<-c(5/60, 10/60, 15/60, 20/60, 30/60, 1, 2, 6)   # durations in hr for idf analysis
  }else{
    duration<-c(Duration[Duration<1]/60,Duration[Duration>=1])
  }
  # ----Definir Periods----
  if(length(Periods)==1){
    Tp<-c(2,3,5,10,25,50,100)
  }else{
    Tp<-Periods                            # return periods in years for plotting idf curves
  }

# ----Corre los diferentes escenarios----
  for(estr in Strategy){

    if(estr==1){ # Solo datos digitalizados
      intensities<-Int.total[-id.info,]
      durations<-duration
      nom.dura<-paste(c(5,10,15,20,30,60,120,360)," min",sep="")
    }else if (estr==2){ # Solo datos Ideam
      intensities<-Int.total[id.info,c(3,5:8)]
      durations<-duration[c(3,5:8)]
      nom.dura<-paste(c(15,30,60,120,360)," min",sep="")
    }else{ # Ambos conjunto de datos
      intensities<-Int.total[,c(3,5:8)]
      durations<-duration[c(3,5:8)]
      nom.dura<-paste(c(15,30,60,120,360)," min",sep="")
    }

# ----Ajusta la distribucion y calcula intesidades por duracion y periodo de retorno----
    nd = length(durations)
    nTp = length(Tp)

    options(warn = 0) # Para mostrar wanings 0
    distri <- list() # Almacena todos los resultados de la funcion FIT_DISTRI
    idf <- matrix(nrow = nd, ncol = nTp)
    M.test.fit <- matrix(NA, nrow = nd,ncol = 8)

    for(i in 1:nd){
      distri[[nom.dura[i]]]<-fitDISTRI(Intensity = intensities[,i], Type ="Gumbel", Plot = Plot, M.fit = Type,
                                       Periods = Tp, Dura = paste(as.character(durations[i]), " min.", sep=""),
                                       Station = Station, CI = CI, iter = iter,
                                       goodtest = goodtest, Resolution = Resolution, SAVE = SAVE)
      idf[i,]<-distri[[i]]$Int.pdf
      M.test.fit[i,]<-distri[[i]]$goodness.fit

    }


  colnames(idf) <- as.character(Tp)
  rownames(idf) <- nom.dura
  colnames(M.test.fit) <- names(distri[[i]]$goodness.fit)
  rownames(M.test.fit) <- nom.dura

  # ----Compute idf equations per each time durations
  Output[[name[estr]]] <- regIDF(Intensity = idf, Periods = Tp, Durations= durations, logaxe = logaxe,
                                 Plot = Plot, Resolution = Resolution, SAVE = SAVE, Strategy = estr,
                                 M.fit = M.fit, Type = Type, name = name, Station = Station)


  if (SAVE) {

    if (file.exists(paste(".", "RESULTS", Station, sep = "/"))) {
      path.result<-paste(".", "RESULTS", Station, sep = "/")
    } else {
      dir.create(paste(".", "RESULTS", Station, sep = "/"), recursive = TRUE)
      path.result <- paste(".", "RESULTS", Station, sep = "/")
    }

    write.xlsx(idf, file = paste(path.result,"/", "IDF_", Station, name[estr], ".xlsx",sep=""),
               sheetName = "IDF.by.PDF")
    write.xlsx(idf, file = paste(path.result,"/", "IDF_", Station, name[estr], ".xlsx",sep=""),
               sheetName = "IDF.by.PDF")
    write.xlsx(M.test.fit, file = paste(path.result, "/", "IDF_", Station, name[estr], ".xlsx", sep = ""),
               sheetName = "goodness.fit", append = TRUE)
    write.xlsx(Output[[name[estr]]]$Coefficients, file = paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
               sheetName = "Coefficients", append = TRUE)
    write.xlsx(Output[[name[estr]]]$Predict,file = paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
               sheetName = "Prediction.by.C-IDF", append = TRUE)
    write.xlsx(Output[[name[estr]]]$test.fit.reg = paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
               sheetName="Performance-IDF.reg", append = TRUE)
    write.xlsx(Output[[name[estr]]]$Confidence.Int = paste(path.result, "/", "IDF_", Station, "_", name[estr], ".xlsx", sep = ""),
               sheetName="Conf.Int-IDF.reg", append = TRUE)

  }


}

  # Salidas de la funcion

    OUT.FINAL <- list(Intesidades = M.idf, Modelo = Output,
                    Ajuste = M.test.fit, Distibucion = distri)
    print("A")
  return(OUT.FINAL)

}
