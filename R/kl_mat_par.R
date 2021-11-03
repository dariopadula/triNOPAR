
#' Estima la funcion de KL para un conjunto de items estimados de forma parametrica
#'
#' Dado un na secuencia de items estimados de forma parametrica, estima la funcion KL para cada uno
#'
#'
#' @param paramsMIRT data frames con los parametros estimados para elmodelo parametrico. Los nombres
#' de las filas de este data rame tienen que ser los nomres de los items.
#' @param grillaEval grilla de puntos donde se va a estimar la funcion
#' @param distTrans distribucion para transformar la grilla en caso de que la grilla este definida en $[0,1]$.

#' @return
#' Un data frame con tantas columnas como items y filas como puntos de la grilla
#' con las estimaciones del KL
#' @export
#' @examples
#' \dontrun{
#' # No correr
#' }



kl_mat_par = function(paramsMIRT,grillaEval,distTrans = qnorm) {
  nomItems = rownames(paramsMIRT)
  kl_res = do.call(cbind,sapply(nomItems,
                                function(xx) {
                                  res = KL_PAR(nomItem = xx,paramsMIRT,grillaEval,distTrans)
                                  list(res)
                                }
  ))

  return(kl_res)
}
