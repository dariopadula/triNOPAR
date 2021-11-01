#' Funcion para generar las funcions gs que se usan como
#' estadistico usando sumas simples o ponderadas sobre los items
#'
#' @param resp Matriz de ceros y unos que representan las respuestas a los items
#' @param grupos es un list que tiene tantas componentes como g's (estadisticos de ordenacion) y en las filas los items que se usan en la suma ponderada o simple
#' @param pesos ponderaciones en la suma
#' @return Un data frame con las respuestas y al final nuevas columnas con las estimaciones de las habilidades resultantes.
#' Una columna por cada dimension. Se identifican con la legra *g* seguido del numero de la dimension.
#' @export
#' @examples
#'\dontrun{
#' # Ejemplo 1 dimension solo suma
#' # itdim1=colnames(respu)[1:80]
#' # itemsgs=list(dim1=itdim1)
#' # pesos=list(dim1=rep(1,80))
#' # respG=estFunG_Simple(respu,itemsgs,pesos)
#'
#'
#' Ejemplo para dos dimensiones
#' itdim1=colnames(respu)[1:50]
#' itdim2=colnames(respu)[45:80]
#' itemsgs=list(dim1=itdim1,dim2=itdim2)
#' pesos=list(dim1=runif(50),dim2=runif(36))
#' respG=estFunG_Simple(respu,itemsgs,pesos)
#' }

estFunG_Simple = function(resp,grupos,pesos){
  ndims=length(grupos)
  scores=c(1:dim(resp)[1])
  for(j in 1:ndims){
    auxit=resp[,grupos[[j]]]
    auxscore=auxit%*%pesos[[j]]
    scores=cbind(scores,auxscore)
  }
  nomCols = c('ID',paste0('g',1:ndims))
  colnames(scores) = nomCols
  respG = as.data.frame(cbind(resp,scores))

  return(respG)
}
