#' FUncion para integrar en el calculo de la funcion de KL para las curvas paramétrica en un punto dado
#'
#' Esta funcion estima es la que se integra en un entorno del punto para estimar KL
#'
#'
#' @param arg es la variable que se va a integrar
#' @param thj es el punto the la grilla de theta donde se va a calcular KL
#' @param a,b son los parámetros de discriminacion y dificultad del modelo parametrico 2PL

#' @return
#' Un vector con la funcion KL estimada en cada thj.
#' @export
#' @examples
#' # No correr


funKL_int_PAR <- function(arg,thj,a,b) {

  evalArg = exp(a*(arg-b))/(1+exp(a*(arg-b)))
  evalFija = exp(a*(thj-b))/(1+exp(a*(thj-b)))

  if(evalFija == 1){
    ff<-evalFija*log(evalFija/evalArg)
  } else{
    ff<-evalFija*log(evalFija/evalArg) + (1 - evalFija)*log((1 - evalFija)/(1 - evalArg))
  }
  return(ff)
}
