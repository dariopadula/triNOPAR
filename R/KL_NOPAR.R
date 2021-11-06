
#' Estima la funcion de KL para las curvas NO paramétrica en una grilla de puntos
#'
#' Esta funcion estima la funcion de KL para para una ICC estimada de forma NO parámetrica en una grilla d epuntos
#'
#'
#' @param iccNP vecor con las estimaciones en una grilla
#' @param sepGrilla separacion entre puntos de la grilla
#' @param entorno porcentaje de puntos para tomar en el entorno


#' @return
#' Un vector con la funcion KL estimada en cada punto de la grilla.
#' @export
#' @examples
#' \dontrun{
#' # No correr
#' }


KL_NOPAR = function(iccNP,sepGrilla,entorno) {

  ###### Calcula el KL
  KL = sapply(1:length(iccNP),function(tt) {
    funKL_int_NOPAR(iccNP,sepGrilla,entorno,jfijo = tt)
  })

  return(KL)
}
