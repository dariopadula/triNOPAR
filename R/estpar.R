
#' Estima los parametros de los items usando MIRT
#'
#' Para una matriz de respuesta, en nuestro caso generada por la funcion GenResp, se estiman los parametros de los items asumiendo uno a mas
#' modelo para los distintos items. Se pueden usar los modelas verdadesros que son los que estan definidos en el banco o
#' se pueden modelar con un modelo distinto al que fueron generados
#'
#'
#' @param datos Matriz de ceros y unos que representan las respuestas a los items
#' @param modelo es lo mismo que se le asigna al parametro model de la funcion mirt
#' @param tipit un vector de caracteres que indica el modelo a aplicar a cada item. Si los items se estiman con distintos modelos, tipit es de
#' largo del numero de items, si se aplica el mismo modelo a todos los items, puede ser de largo uno. La forma de identificar los modelos de
#' los items tiene que ser de la forma que los usa la funcion mirt.
#' @return Retorna una lista con dos componentes. La primera es un data frame con los parametros estimados por la funcion mirt, la segunda
#' contiene informacion de la convergencia del modelo
#' @export
#' @examples
#'
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
#' ### Estima los parametros asumiendo unidimensionalidad y que los items son todos 2PL
#' modelo=1
#' tipit<-c(rep("2PL",25)) # En este caso podria se solo tipit = "2PL"
#' paramsEst = estpar(respu,modelo,tipit)

#Modificacion Mario 7 de octubre habilidades

estpar<-function(datos,modelo,tipit){
  estima = mirt::mirt(datos,model = modelo,itemtype=tipit,SE = TRUE)
  conv= c(mirt::extract.mirt(estima, what="converged") ,mirt::extract.mirt(estima, what="secondordertest") )
  para= mirt::coef(estima, simplify = TRUE, IRTpars = TRUE)
  hab=mirt::fscores(estima, method = 'ML', full.scores=TRUE,full.scores.SE = F)
  sal =list(parametros = para,habilidades=hab,convergencia = conv)
  return(sal)
}
