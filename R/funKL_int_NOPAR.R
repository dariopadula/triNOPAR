
#' Estima la funcion paroxima la integral en un entorno de un punto para estimar la funcion KL
#'
#' Esta funcion define el integrando para estimar la funcion de KL  para una ICC estimada de forma NO parámetrica
#' o la curva isotonica en un punto fijo de la grilla
#'
#'
#' @param iccNP Estimacion de la ICC no parametrica en determinados puntos
#' @param sepGrilla Separecion entre los puntos de la grilla
#' @param entorno porcentaje que indica cuantos puntos tomar en el entorno del punto a estimar
#' @param jfijo posicion de la grilla donde se va a estimar la funcion

#' @return
#' Un vector con la funcion KL estimada en cada posicion *jfijo*.
#' @export
#' @examples
#' \dontrun{
#' # No correr
#' }



funKL_int_NOPAR <- function(iccNP,sepGrilla,entorno,jfijo) {

  totPuntos = length(iccNP)
  posInf = max(1,jfijo - floor(entorno*totPuntos))
  posMax = min(totPuntos,jfijo + floor(entorno*totPuntos))

  if(is.na(iccNP[1])) iccNP[1] = 0
  if(is.na(iccNP[totPuntos])) iccNP[totPuntos] = 1

  evalArg = iccNP[posInf:posMax]
  evalFija = iccNP[jfijo]

  if(evalFija == 1){
    ff<-sum(evalFija*log(evalFija/evalArg))*sepGrilla
  } else{
    ff<-sum(evalFija*log(evalFija/evalArg) + (1 - evalFija)*log((1 - evalFija)/(1 - evalArg)))*sepGrilla
  }
  return(ff)
}
