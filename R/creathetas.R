#' Simula las habilidades
#'
#' Dado un banco de items, esta funcion extrae items generados con distintos modelos. La cantidad de items de cada tipo de diseno se
#' determina en la matriz diseno.
#'
#' @param distrth Distribucion con la que se va a generar las habilidades
#' @param ... Parametros especificos de la funcion que se elija para simulra las habilidades
#' @return Un vector de habilidades de largo "n" donde el parametro n se define dentro de "..."
#' @export
#' @examples
#' \dontrun{
#' # thetas=creathetas("Normal",n=1000,mean=0,sd=1)
#' }


creathetas<-function(distrth,...){
  FUN = dplyr::case_when(distrth=="Normal"~parse(text="rnorm"),
                  distrth=="Uniforme"~parse(text="runif"),
                  distrth=="t"~parse(text="rt"),
                  distrth =="F"~parse(text="rf"),
                  distrth=="Chi2 "~parse(text="rchisq"),
                  distrth=="Beta"~parse(text="rbeta"),
                  distrth=="Gamma"~parse(text="rgamma"),
                  distrth=="Logistica"~parse(text="rlogis"),
                  distrth=="Lognormal"~parse(text="rlnorm"),
                  distrth=="Normtruncada"~parse(text="rtruncnorm"),
                  distrth=="Normsesgada"~parse(text="rsn"))
  thetas=eval(FUN)(...)
  return(thetas)
}
