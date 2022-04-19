
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
#' @param puntuni grilla en el intervalo 0 1 donde estan estimados los puntos de la funcion icc no parametrica
#' @param shiftFun es la funcion que lleva a grilla 0 1 a la grilla de la distribucion que se indique, por defecto
#' se pone la normal

#' @return
#' Un data frame con tantas columnas como items y filas como puntos de la grilla
#' con las estimaciones del KL
#' @export
#' @examples
#' \dontrun{
#' # No correr
#' # KL_MAT_NP = kl_mat_NOpar(iccNP_mat,sepGrilla,entorno)
#' }



kl_mat_NOpar = function(iccNP_mat,sepGrilla,entorno,puntuni,shiftFun = qnorm) {

  ##### Calcula la ICCiso en una grilla de theta
  ## XXX Nuevo
  thet = shiftFun(puntuni)
  tmax = max(abs(thet))
  thetaG = seq(-tmax,tmax,length = length(puntuni))

  #################################################
  #################################################

  nomItems = colnames(iccNP_mat)
  kl_res = do.call(cbind,sapply(nomItems,function(xx) {


    ##### calculo la curva en la grilla de thetas
    ## XXX Nuevo
    iccis<-approx(thet,iccNP_mat[,xx],thetaG)$y
    res = KL_NOPAR(iccis,sepGrilla = abs(thetaG[1] - thetaG[2]),entorno)
    ##################################################

    ##### Reconstruyo los puntos anteriores que corresponden al int 01
    ## XXX Nuevo
    res<-approx(thetaG,res,thet)$y
    res[is.na(res)] = 0
    #################################################
    #################################################

    # res = KL_NOPAR(iccNP_mat[,xx],sepGrilla,entorno)
    list(res)
  }
  ))

  return(kl_res)
}

