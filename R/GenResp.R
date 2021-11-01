#' Genera respuestas dada la habilidad theta
#'
#' Para un banco de items parametrico simula una matriz de ceros (respuesta incorrecta) y unos (respuesta correcta) para cada uno de los
#' valores del vector theta.
#'
#' @param banco es un data frame con tantas filas como items tenga el Banco, en este dataframe se indica el modelo y los parámetros del ítem.
#' @param thetas es un vector con habilidades simuladas con la funcion *creathetas*.
#' @return Retorna una matriz con tantas filas como habilidades (individuos) y tantas columnas como items haya en el banco.
#' @export
#' @examples
#'\dontrun{
#' # Genero el banco vacio
#' bancoIT = genBancoDF(nomFijas = c('NombreIt','Modelo'),nparam = 15)
#' m2pl=matrix(c(0.2,1,-3,3),byrow = TRUE,ncol=2,nrow=2)
#' bancoIT = genitmodelo(100,"2PL",m2pl,bancoIT)
#' # Especifico que voy a extraer del banco completo, bancoTAI, 25 items de dos parametros
#' diseno=matrix(c(25,"2PL"),ncol=2,byrow = FALSE)
#' bancTeor = bancoTAI(bancoIT,diseno)
#'# Genero las habilidades
#' thetas=creathetas("Normal",n=1000,mean=0,sd=1)
#' respu=GenResp(banco = bancTeor,thetas)
#' }


GenResp<-function(banco,thetas){

  modelos=banco$Modelo
  D=1
  resp <- matrix(nrow = length(thetas), ncol = dim(banco)[1])
  colnames(resp)=paste0(banco$NombreIt,banco$Modelo)

  for(i in 1:length(modelos)){
    # Para cada modelo determina los parametros
    if (modelos[i]=="1PL"|modelos[i]=="2PL"|modelos[i]=="3PL") {
      a=banco$P1[i]
      b=banco$P2[i]
      c=banco$P3[i]
      if(is.na(c)) c=0
      Prob <- c + (1 - c) /(1 + exp(-D * a * (thetas - b)))
    }
    if(modelos[i]=="CUB"){
      a=banco$P1[i]
      b=banco$P2[i]
      c=banco$P3[i]
      Prob <- 1/ (1 + exp(-(a + b*thetas + c*thetas^3 )))
    }
    aux <- runif(length(thetas), 0, 1)
    resp[,i]=as.integer(Prob>aux)
  }
  return(resp)

}
