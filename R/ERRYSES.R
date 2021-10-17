
#' Funcion para estimar el error y el sesgo de las estimaciones del TAI
#'
#'
#'
#' @param simData matriz o data frame con los valores reales de habilidad y los estimados
#' que seria la salida de la funcion TAIgeneric.
#' @param grilla grilla de puntos de habilidad para estimar el RMSE y el sesgo local
#' @return Returna una lista con:
#' * **RMGEN**: valores del sego y RMSE para la simulacion ingresada
#' * **RMBIATH**: estimaciones del RMSE y el sego para una grilla de valores de $\theta$
#' @export
#' @examples
#' # No correr
#'


ERRYSES = function(simData,grilla = seq(1:100)/100) {
  RMSE = sqrt(sum((simData[,1]-simData[,4])^2/nrow(simData)))
  BIAS = sum((simData[,1]-simData[,4])/nrow(simData))

  RMGEN = data.frame(RMSE,BIAS)

  #### Errores y sesgos locales
  ngrilla = length(grilla)
  RMBIATH = as.data.frame(matrix(NA,ncol=3,nrow=ngrilla))
  colnames(RMBIATH) = c('Grilla','RMSE','BIAS')

  for(i in 2:ngrilla){
    ii=grilla[i] - 1/ngrilla
    iii=grilla[i] + 1/ngrilla
    suaux= simData[,4] < iii & simData[,4] > ii
    RMBIATH[i,1]=grilla[i]
    RMBIATH[i,2]=sqrt(sum((simData[suaux,1]-simData[suaux,4])^2/length(sujtai[suaux])))
    RMBIATH[i,3]=sum((simData[suaux,1]-simData[suaux,4])/length(sujtai[suaux]))
  }

  return(list(raiz=RMGEN,rportheta=RMBIATH))
}
