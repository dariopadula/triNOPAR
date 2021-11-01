#' Para una grilla de puntos calcula la curva caracteristica del item ICC
#'
#'
#'
#' @param grilla valores de theta donde se evalua la funcion
#' @param a parametro de discriminacion del la ICC
#' @param b parametro de dificultad del la ICC
#' @param c parametro de pseudo azar de la ICC
#' @return Retorna un vector con la funcion evaluada en los puntos de la grilla
#' @export
#' @examples
#' \dontrun{
#' a = 1; b = 0; c = 0.1
#' grilla = dnorm(seq(0,1,0.01))
#' p = iccFun(a,b,c,grilla)
#' }



iccFun = function(a,b,c,grilla) {
  p = c + (1-c)/(1 + exp(-a*(grilla - b)))
  return(p)
}
