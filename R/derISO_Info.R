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
#' \dontrun{
#' # No correr
#' # icciso = icciso1[['resfin']]
#' # puntuni = icciso1[['puntos']]
#' # resDer = derISO_Info(icciso,puntuni,nucleo = normal)
#' }



derISO_Info = function(iccis,puntuni,nucleo,hd = NULL) {

  ###### Calcula la derivada
  if(is.null(hd)) hd= 0.9*length(iccis)^(-1/5)*min(sd(iccis,na.rm = T),
                                                   (summary(iccis)[5]-summary(iccis)[2])/1.364)

  iccisoFila = matrix(iccis,length(iccis),length(iccis),byrow = T)
  iccisoCol = t(iccisoFila)
  difer = iccisoFila - iccisoCol

  deriso =  (length(puntuni)*hd)/(apply(difer,2,function(xx) sum(nucleo(xx,h = hd,th = 0)$res,na.rm = T)))

  ###### Calcula la funcion de informacion
  inform = (deriso^2)/(iccis*(1 - iccis))

  ###### Entorno
  refInf = max(min(iccis),0.4)
  refSup = min(max(iccis),max(0.6,refInf))

  pos40 = which.min(abs(iccis - refInf))
  pos60 = which.min(abs(iccis - refSup))

  if(pos40 == pos60) pos60 = pos60 + 20
  ##### Me quedo con esos valores y el resto lo mando a cero
  informEntorno = inform*0
  informEntorno[pos40:pos60] = inform[pos40:pos60]

  return(list(deriv = deriso,Info = inform,informEntorno = informEntorno))
}
