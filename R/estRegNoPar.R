#' Estima ICC usando regresion NO parametrica
#'
#' Esta funcion estima las curvas caracteristicas de los items usando regresion no parametrica tanto para items
#' unidimensionales como multidimensionales.
#'
#'
#' @param items Puede ser un vector de strings con los nombres de los items o las posiciones.
#' @param h Tamano de la ventana a usar dentro del nucleo
#' @param th_use Referencia a la estimacion de los thetas a usar del data frame proveniente de la funcion
#' *estthetaNP*, puede tomar los valores: pcg (usa los percentiles) o dtg (usa la tranformacion a alguna distribucion)
#' @param test Data frame o matriz con las respuestas de los individuos y las estimaciones de las habilidades. Salida
#' de la funcion *estthetaNP*.
#' @param puntos Es una grilla de puntos donde se va a hacer la estimacion no parametrica de las ICC, la grilla
#' se define dependiendo de las dimensiones de los items.
#' @param nucleo Es es la funcion Kernel utilizada para estimar la ICC
#' @param sigma Es un parametros para el nucleo, es un escalar positivo si la diemnsion es 1 y es una matriz en el caso
#' de mas de una dimension.
#' @return Una lista donde con:
#' * NPICC: matriz con la estimacion de las icc no parametrica
#' * pesos: los pesos usados para la regresion
#' * puntos: los puntos donde se realizo la estimacion
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
#' # Estimo las habilidades para una dimension usando AF
#' respG = estFunG(resp = respu,rotacion = c('oblimin'),corr = c('tet'),ndims = 1)
#' # Estimo los thetas de forma no parametrica y usando la transformacion a la distribucion normal
#' thetaest = estthetaNP(scores = respG,Dth = qnorm)
#'
#' # Estima una icc para un item especifico
#' NumItem = 6
#' # ICCNP = estRegNoPar(items = NumItem,h = 0.01,th_use = 'pcg',test = thetaest,puntos = grilla,nucleo=normal,sigma = 1)
#' # Estima para un conjunto de items
#' # ICCNP_todos = estRegNoPar(items = c(1:25),h = 0.01,th_use = 'pcg',test = thetaest,puntos = grilla,nucleo=normal,sigma = 1)


estRegNoPar = function(items,h = 0.2,th_use = 'pcg',test,puntos,nucleo=normal,sigma){

  varsTh = colnames(test)[grep(th_use,colnames(test))]
  dimension = length(varsTh)

  th=test[,varsTh]

  if(dimension==1){
    w = do.call(cbind,sapply(puntos,nucleo,h,th))
  }else{

    if(length(h) == 1) h = diag(dimension)*h
    if(length(sigma) == 1) sigma = diag(dimension)*sigma

    w = do.call(cbind,sapply(1:nrow(puntos),function(xx) {

      th = as.matrix(test[,varsTh])
      peval = matrix(c(puntos[xx,c(1)],puntos[xx,c(2)]),nrow(test),2,byrow = T)

      argumento<-(peval - th)%*%ginv(h)

      nunu<-(1/abs(det(h)))*nucleo(argumento, sigma)$res

      list(nunu)

    }))
  }

  pesos=prop.table(w,margin=2)
  CCInp = t(pesos)%*%as.matrix(test[,items])

  return(list(NPICC=CCInp, pesos=pesos,puntos = puntos))
}
