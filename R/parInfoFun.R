#' Para una grilla de puntos en el intervalo (0,1) calcula la funcion de informacion del item
#'
#'
#'
#' @param grilla valores de theta donde se evalua la funcion
#' @param a parametro de discriminacion del la ICC
#' @param b parametro de dificultad del la ICC
#' @param c parametro de pseudo azar de la ICC
#' @return Retorna un vector con la funcion de informacion evaluada en los puntos de la grilla
#' @export
#' @examples
#' \dontrun{
#' a = 1; b = 0; c = 0.1
#' grillaUni = grilla =seq(0,1,0.01)
#' infor = parInfoFun(a,b,c,grillaUni,trnfFun = qnorm)
#' }

parInfoFun = function(a,b,c,grillaUni,trnfFun = qnorm) {
### Pasa de la grilla en 0-1 a la distribucion de los thetas
  grilla = eval(trnfFun(grillaUni))
### Capcula p y q
  p = iccFun(a,b,c,grilla)
  q = 1 - p
### Calcula la funcion de informacion
  InfFun = (a^2)*((p - c)^2/(1 - c)^2)*(q/p)
  return(InfFun)
}
