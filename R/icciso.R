#' Funcion para estimar las curvas isotonas en dimension 1
#'
#'
#'
#' @param icc1 estimados de la ICC no parametrica icc1 y puntosicc deben tener la misma longuitudla y curva estimada en el intervalo [0,1]
#' @param hd ventana para estimacion de la densidad
#' @param thetaiso vector donde se calculara la curva isotona
#' @param nt puntos en la grilla
#' @param puntosicc indica los valores donde esta calculada la icc1
#' @param nucleod nucle a usar para estimar la densidad
#' @return Returna un sub banco de items muestreados del banco original.
#' @export
#' @examples
#' \dontrun{
#' # No correr
#' #icciso1 = icciso(icc1,hd,theta,nt,puntosiso,nucleod)
#' }


icciso=function(icc1,hd,thetaiso,nt,puntosicc,nucleod){

  #Arma la grilla para aplicar la ec 2 del cap Mario
  t=c(1:nt/nt)
  t[length(t)]=0.9999
  fes=numeric(length(t))
  auxpu=c(0,puntosicc,1)
  auxic=c(icc1[1],icc1,1)
  #Se calcula la icc en los puntos i/nt
  ifi1<-approx(auxpu,auxic,t)$y


  fes = sapply(1:length(t),function(ii) {
    alu = sapply(ifi1, function(yy) {
      integrand=function(x){nucleod((yy-x)/hd,1,0)$res}
      resAlu = integrate(integrand,-Inf,t[ii],subdivisions=10,rel.tol = 0.03,abs.tol =0.05,stop.on.error = FALSE)[[1]]
    })
    resFes = (1/(length(ifi1)*hd))*sum(alu)
  })

  t=c(max(t[1]-((t[2]-t[1])/(fes[2]-fes[1]))*fes[1],0),t,1)
  fes[fes>=0.995]<-0.9995
  ### Arregla temas de borde
  fes=c(0.0005,fes,0.9995)

  ### Toma la reflexion respecto  a la bisectriz delcuadrado unidad
  resfin<-approx(fes,t,thetaiso)$y

  return(list(resfin=resfin,puntos=thetaiso))
}

