#' Estima las ventanas para items unidimensionales
#'
#' Esta funcion utiliza la libreria np
#'
#'
#' @param items Puede ser un vector de strings con los nombres de los items o las posiciones.
#' @param th_use Referencia a la estimacion de los thetas a usar del data frame proveniente de la funcion
#' *estthetaNP*, puede tomar los valores: pcg (usa los percentiles) o dtg (usa la tranformacion a alguna distribucion)
#' @param test Data frame o matriz con las respuestas de los individuos y las estimaciones de las habilidades. Salida
#' de la funcion *estthetaNP*.
#' @param nucleodes Es es la descripción del Kernel utilizado.
#' Los valores permitidos son  gaussian, epanechnikov, o uniform
#' @param muestra cantidad de casos tomados para la estimacion, por defecto son todos.
#' @return Una vector con los anchos de ventana vantanas
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
#'
#' # Estima la ventana para un item especifico
#' h=ventana1D(items = 1,th_use = 'dtg',test = thetaest,nucleodes="gaussian",1000)
#' }

ventana1D = function(items,th_use = 'pcg',test=NULL,nucleodes="gaussian",muestra="TODO"){

varsTh = colnames(test)[grep(th_use,colnames(test))]
dimension = length(varsTh)

if(muestra!="TODO"){
  test=test[sort(sample(dim(test)[1],muestra)),]
  th=test[,varsTh]
  }


if(dimension==1){
    h=rep(NA,length(items))
    for(j in items){
      haux=npregbw(formula=test[,j]~th,ckertype=nucleodes)
      h[j]=haux$bandwidth$x }
}else {stop("Dimension mayor a 1")}


return(h)
}
