#' Estima la derivada primera de la curva ISO y la funcion de informacion
#'
#' Esta funcion estima para las curvas isotonicas la derivada primera
#'
#'
#' @param icciso un vector con la estimacion isotonica
#' @param puntuni grilla donde esta evaluada la curva iso
#' @param nucleo alguna funcion de nucleo de dentro de las opciones en **nucleosAll**
#' @param hd ventan para la estimacion, por defecto es NULL, en ese caso se estima segun formula dada, si se engresa un valor,
#' se toma ese.
#' @return Una lista donde con:
#' * deriv: -Derivada de la curva isotonica
#' * Info: Funcion de informacion de la curva isotonica
#' @export
#' @examples
#' # No correr
#' # icciso = icciso1[['resfin']]
#' # puntuni = icciso1[['puntos']]
#' # resDer = derISO_Info(icciso,puntuni,nucleo = normal)



derISO_Info = function(icciso,puntuni,nucleo,hd = NULL) {

  ###### Calcula la derivada
  if(is.null(hd)) hd= 0.9*length(icciso)^(-1/5)*min(sd(icciso,na.rm = T),(summary(icciso)[5]-summary(icciso)[2])/1.364)

  iccisoFila = matrix(icciso,length(icciso),length(icciso),byrow = T)
  iccisoCol = t(iccisoFila)
  difer = iccisoFila - iccisoCol

  deriso =  (length(puntuni)*hd)/(apply(difer,2,function(xx) sum(nucleo(xx,h = hd,th = 0)$res,na.rm = T)))

  ###### Calcula la funcion de informacion
  inform = (deriso^2)/(icciso*(1 - icciso))

  return(list(deriv = deriso,Info = inform))
}
