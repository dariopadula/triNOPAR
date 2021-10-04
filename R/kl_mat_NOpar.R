
#' Estima la funcion de KL para un conjunto de items estimados de forma NO parametrica
#'
#' Dado una secuencia de items estimados de forma NO parametrica (puede ser la isotonica),
#' estima la funcion KL para cada uno.
#'
#'
#' @param iccNP_mat data frames con las estimaciones de las ICC no parametricas. Los nombres de las columnas tienen que ser los nombres
#' de los items.
#' @param sepGrilla distancia entre los puntos consecutivos de la grilla. Mismos que la funcion **funKL_int_NOPAR**
#' @param entorno porcentaje que indica cuantos puntos tomar en el entorno del punto a estimar.
#' Mismos que la funcion **funKL_int_NOPAR**

#' @return
#' Un data frame con tantas columnas como items y filas como puntos de la grilla
#' con las estimaciones del KL
#' @export
#' @examples
#' # No correr
#' # KL_MAT_NP = kl_mat_NOpar(iccNP_mat,sepGrilla,entorno)



kl_mat_NOpar = function(iccNP_mat,sepGrilla,entorno) {
  nomItems = colnames(iccNP_mat)
  kl_res = do.call(cbind,sapply(nomItems,function(xx) {
    res = KL_NOPAR(iccNP_mat[,xx],sepGrilla,entorno)
    list(res)
  }
  ))

  return(kl_res)
}
