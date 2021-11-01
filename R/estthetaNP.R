#' Etimacion NO parametrica de las habilidades
#'
#' Esta funcion estima las habilidades de los individuos tomando como insumo algun estadistico (suma de scores, thetas estimados
#' con la funcion estFunG u otros).
#' ROmpe los scores con empates, calcula el percentil para cada individuo y realiza una transformacion llevando la uniforme al dominio
#' de alguna otra distribucion, por ejemplo la normal
#'
#'
#' @param scores Data frame resultante con las respuestas a los items (no importan) y con columnas que contengan la estimacion de
#' los scores de los individuos por agun metodo, ya sea suma, cargas factoriales etc.
#' @param Dth distribucion que se le va a aplicar a los percentiles para transformarlos. Hay que poner alguna funcion del estilo
#' qnorm por ejemplo.
#' @param ... Parametros adicianales que asociados a la funcion Dth. que se elija para hacer la transformacion.
#' @return Un data frame con las mismas columnas que tenia *scores*, pero se le agragan columnas identificadas como: seg (score g sin empates),
#' pcg (percentiles del score g sin empates) y dtg (es la transformacion aplicando la funcion Dth a pcg).
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
#' # Estimo las habilidades para una dimension usando AF
#' respG = estFunG(resp = respu,rotacion = c('oblimin'),corr = c('tet'),ndims = 1)
#' # Estimo los thetas de forma no parametrica y usando la transformacion a la distribucion normal
#' thetaest = estthetaNP(scores = respG,Dth = qnorm)
#' }


estthetaNP = function(scores,Dth = qnorm,...) {


  # encuentra los nombrs de las gs
  varsG = colnames(scores)[grep('^g',colnames(scores))]
  # recorre cada g para, desemparar, percentiles e inversa
  for(ii in varsG) {

    totRow = scores[,ii]
    orden = totRow[sort.int(totRow,index.return = T)$ix]
    dif = diff(orden)
    dmin = min(dif[dif > 0])
    scoreSE = totRow + runif(length(totRow),0,dmin)
    ## VARIABLE SIN EMPATES
    nomSE = gsub('^g','seg',ii)
    scores[,nomSE] = scoreSE
    ## PERCENTIL
    orden = order(scoreSE)
    nomPerc = gsub('^g','pcg',ii)
    scores[orden,nomPerc] = (1:length(scoreSE))/length(scoreSE)
    ## Arregla el ultimo para que no de infinito
    scores[dplyr::last(orden),nomPerc] = scores[dplyr::last(orden),nomPerc] - runif(1,0,1/length(scoreSE))
    ## Aplica la inversa
    nomDist = gsub('^g','dtg',ii)
    aux = scores[,nomPerc]
    scores[,nomDist] = Dth(aux,...)

  }
  return(scores)

}
