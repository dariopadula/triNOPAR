
#' Estima la funcion de KL para las curvas paramétrica en una grilla d epuntos
#'
#' Esta funcion estima la funcion de KL para para una ICC estimada de forma parámetrica en una grilla d epuntos
#'
#'
#' @param nomItem string con el nombre del item a estimar
#' @param paramsMIRT data frames con los parametros estimados para elmodelo parametrico
#' @param grillaEval grilla de puntos donde se va a estimar la funcion
#' @param distTrans distribucion para transformar la grilla en caso de que la grilla este definida en $[0,1]$.

#' @return
#' Un vector con la funcion KL estimada en cada punto de la grilla.
#' @export
#' @examples
#' \dontrun{
#' # No correr
#' }


KL_PAR = function(nomItem,paramsMIRT,grillaEval,distTrans = qnorm) {

  if(grillaEval[1] == 0) grillaEval[1] = (grillaEval[2]-grillaEval[1])*0.1
  if(grillaEval[length(grillaEval)] == 1) grillaEval[length(grillaEval)] = grillaEval[length(grillaEval)] - (grillaEval[length(grillaEval)]-grillaEval[length(grillaEval) - 1])*0.1

  if(!is.null(distTrans)) {
    grilla = distTrans(grillaEval)
  } else {
    grilla = grillaEval
  }

  #### Parametros
  a = paramsMIRT[nomItem,'a']
  b = paramsMIRT[nomItem,'b']

  KL = sapply(grilla,function(tt) {
    integrate(funKL_int_PAR,tt-0.1,tt+0.1,subdivisions=10,rel.tol = 0.03,
              abs.tol =0.05,stop.on.error = FALSE,thj = tt,a=a,b=b)$value
  })

  return(KL)
}


