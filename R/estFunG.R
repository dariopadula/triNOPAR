#' AF para estimar las habilidades de los individuos
#'
#' Esta funcion utiliza analisis factorial (AF) para estimar las habilidades para cada fila de la matriz de respuesta
#' a un determinado conjunto de items.
#' Se pueden asumir una o mas dimensiones y utilizar distintos tipos de rotacion. Las habiliddes se estiman
#' segun las cargas factoriales de los individuos a los factores.
#'
#' @param resp Matriz de ceros y unos que representan las respuestas a los items
#' @param rotacion se indiqua que rotacion usar para el AF. Puede ser una de estas: varimax (defecto), oblimin o none (ninguna)
#' @param corr correlacion a utilizar: tet (tetracorica), poly (policorica), mixed, cor (Pearson)
#' @param ndims numero de dimensiones a estimar.
#' @return Un data frame con las respuestas y al final nuevas columnas con las estimaciones de las habilidades resultantes del AF.
#' Una columna por cada dimension. Se identifican con la legra *g* seguido del numero de la dimension.
#' @export
#' @examples
#'\dontrun{
#' # Genero el banco vacio
#' # bancoIT = genBancoDF(nomFijas = c('NombreIt','Modelo'),nparam = 15)
#' # m2pl=matrix(c(0.2,1,-3,3),byrow = TRUE,ncol=2,nrow=2)
#' # bancoIT = genitmodelo(100,"2PL",m2pl,bancoIT)
#' # Especifico que voy a extraer del banco completo, bancoTAI, 25 items de dos parametros
#' # diseno=matrix(c(25,"2PL"),ncol=2,byrow = FALSE)
#' # bancTeor = bancoTAI(bancoIT,diseno)
#'# Genero las habilidades
#' # thetas=creathetas("Normal",n=1000,mean=0,sd=1)
#' # respu=GenResp(banco = bancTeor,thetas)
#' # Estimo las habilidades para una dimension
#' # respG = estFunG(resp = respu,rotacion = c('oblimin'),corr = c('tet'),ndims = 1)
#' }



estFunG = function(resp,
                   rotacion = c('varimax','oblimin','none'),
                   corr = c('tet','poly','mixed','cor'),
                   ndims = 1) {


  ## Se queda con la primera por defecto
  if(length(rotacion) > 1) rotacion = dplyr::first(rotacion)
  if(length(corr) > 1) corr = dplyr::first(corr)
  # Genera las gs
  fit = psych::fa(resp,nfactors = ndims,rotate = rotacion,cor = corr)
  # Pone el ID
  scores = data.frame(ID = 1:nrow(resp),fit$scores)
  # Nombra el id y las gs
  nomCols = c('ID',paste0('g',1:ndims))
  colnames(scores) = nomCols

  respG = cbind(resp,scores)

  return(respG)
}
