
#' Genera un data frame vacio con un determinado numero de calumnas fijas de characters y un
#' conjunto de columnas numericas
#'
#'
#'
#' @param nomFijas es un vector con los nombres de las columnas que seran character. Los nombres de estas columnas
#' se llamaran igual lo que se ingrese en nomFijas
#' @param nparam entero que indica cuantas columnas numericas se van a cargar. Por defecto estas columnas se nombran con la letra
#' P y un numero secuencial.
#'
#' @return un data frame vacio.
#' @export
#' @examples
#' bancoIT = genBancoDF(nomFijas = c('NombreIt','Modelo'),nparam = 15)


genBancoDF = function(nomFijas = c('NombreIt','Modelo'),nparam = 15) {
  nomPar = paste0('P',1:nparam)
  colNames = c(nomFijas,nomPar)
  banco = stats::setNames(data.frame(matrix(0,nrow = 0,ncol = length(colNames))),colNames)
  for(ii in nomFijas) banco[,ii] = as.character(banco[,ii])
  for(ii in nomPar) banco[,ii] = as.double(banco[,ii])

  return(banco)
}
