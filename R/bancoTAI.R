#' Selecciona un banco de items al azar
#'
#' Dado un banco de items, esta funcion extrae items generados con distintos modelos. La cantidad de items de cada tipo de diseno se
#' determina en la matriz diseno.
#'
#' @param banco un banco de items.
#' @param diseno Una matriz de dos columnas y tantas filas como modelos distintas se quieran extraer del banco
#' @return Returna un sub banco de items muestreados del banco original.
#' @export
#' @examples
#' bancoIT = genBancoDF(nomFijas = c('NombreIt','Modelo'),nparam = 15)
#' m2pl=matrix(c(0.2,1,-3,3),byrow = TRUE,ncol=2,nrow=2)
#' bancoIT=genitmodelo(100,"2PL",m2pl,bancoIT)
#' diseno=matrix(c(25,"2PL"),ncol=2,byrow = FALSE)
#' bancTeor = bancoTAI(bancoIT,diseno)
#' este es porque Dario saco el otro


# Este cambio lo hizo dario

bancoTAI<-function(banco,diseno){

  nti = nrow(diseno)
  banTAI = NULL
  for(i in 1:nti){
    aux = banco %>% dplyr::filter(Modelo == diseno[i,2]) %>%
      dplyr::sample_n(size = as.numeric(diseno[i,1]), replace = FALSE)

    banTAI = rbind(banTAI,aux)
  }
  return(banTAI)
}
