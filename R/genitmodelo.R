#' Genera el banco de items para diferentes modelos con sus respectivos parametros
#'
#'
#'
#' @param caitems es la cantidad de items a ingresar en el banco.
#' @param modelo un modelo pertenciente a los siguientes por ahora 1PL = un parametro, 2PL = dos parametros, 3PL = tres parametros o CUB =
#' modelo cubico
#' @param parametros es una matriz que tiene dos columnas y tantas filas como parametros. Determina los limites inferiores y superiores para
#' simular los parametros.
#' @param banco es un data frame donde se van cargando los parametros y los nombres de los items. Este data frame es el generado por la funcion
#' genBancoDF
#' @return El mismo data frame con los nombres de los items y los parametros agregados debajo de los que ya existian.
#' @export
#' @examples
#' \dontrun{
#' # bancoIT = genBancoDF(nomFijas = c('NombreIt','Modelo'),nparam = 15)
#' # m2pl=matrix(c(0.2,1,-3,3),byrow = TRUE,ncol=2,nrow=2)
#' # bancoIT=genitmodelo(100,"2PL",m2pl,bancoIT)
#' }



genitmodelo=function(caitems,modelo,parametros,banco){
  #Primero comprueba que el modelo este dentro de los especificados

  control=dplyr::case_when(modelo=="1PL"~ TRUE, modelo=="2PL"~ TRUE,modelo=="3PL"~ TRUE,modelo=="CUB"~ TRUE)
  if(is.na(control)){stop("No existe este modelo")}
  if(ncol(parametros)==2){control=TRUE}else{stop("Numero incorrecto de columnas en la matriz de parametros")}
  pmod=dplyr::case_when(modelo=="1PL"~ 2, modelo=="2PL"~ 2,modelo=="3PL"~ 3, modelo=="CUB"~3)
  if(nrow(parametros)==pmod){control=TRUE}else{stop("Numero incorrecto de parametros para el modelo especificado")}
  filasbanco=dim(banco)[1]
  inicio=filasbanco+1
  fin=filasbanco+caitems
  for(j in inicio:fin){
    banco[filasbanco+1,1]=paste0("IT",j)
    banco[filasbanco+1,2]=modelo
    capara=dim(parametros)[1]
    colpar=3
    for(i in 1:capara){
      banco[filasbanco+1,colpar]=runif(1,parametros[i,1],parametros[i,2])
      colpar=colpar+1
    }
    filasbanco=filasbanco+1
  }
  return(banco)
}
