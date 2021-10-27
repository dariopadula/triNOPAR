


#' Funcion para simular la aplicacion de un TAI para items NO parametricos
#'
#'
#'
#' @param sujtai Habilidades verdaderas de los estudiantes con distribucion normal
#' @param epsilon error maximo en la estimacion con respecto a la anterior
#' @param minit minima cantidad de items a aplicar en el TAI
#' @param maxit maxima contidad de items a aplicar
#' @param curvaNOPAR curvas estimadas de forma no parametrica (isotonicas o por regresion)
#' @param parametros parametros verdaderos de los items
#' @param itemsSelec Metodo utilizado para seleccionar los items:
#'
#' * *InfoIso*: pseudo informacion
#' * *KL*: Kullback Leiber
#' * *ESH*: Esperanza de la entropia de shannon
#' * *random*: Seleccion aleatoria
#' @param matrizSelect en los casos de usar el metodo de seleccion InfoIso o KL
#' aqui se pone la matriz de informacion o de KL para todos los items en cada valor de la grilla
#' @param seqTheta Es la grilla de theta donde estan calculadas las matrices de info y KL
#' @return Returna una matriz con el verdadero valor de theta, el estimado,
#' el error cometido y el total de items aplicado.
#' @export
#' @examples
#' # No correr
#'


######## FUNCION GENERICA PARA APLICAR TAI
TAIgeneric <- function(sujtai,
                       epsilon, # Error máximo
                       minit, # Minimo numero de items a aplicar
                       maxit, # Maximo numero de items a aplicar
                       curvaNOPAR, # CUrvas caracterisiticas de los items NO PAR ISOTONICA
                       parametros, # Data frame con los parámetros de los items verdaderos
                       itemsSelec = c('InfoIso','KL','ESH','Random'), # Criterio de seleccion de items
                       matrizSelect = MINOPARAUX, # Matriz para elegir los items (solo se pone en el caso de InfoIso y KL)
                       # El nombre de las columnas tiene que ser el nimbre de los items
                       seqTheta  # Secuencia de thetas en formato uniforme (tiene que tener una correspondencia con
                       # los valores que se calculan las matrices de informacion y KL)
)
{
  estima=matrix(NA,nrow=length(sujtai),ncol=4)
  for(s in 1:length(sujtai)){
    if((s %% 1000) == 0) print(s)

    metodoSel = itemsSelec[1] # TOma el primero por defecto

    ## XXX CAMBIO
    rownames(parametros) = colnames(curvaNOPAR) ## Le pongo el nombre de las filas (las columnas de la info)
    resp=numeric(0)
    vero=1
    th=sujtai[s]
    th_est=runif(1,0.3,0.7)
    delta=1
    control=0
    seleit=NULL

    if(metodoSel %in% c('InfoIso','KL')) matSelAux = matrizSelect
    #seleccion del item
    while(delta > epsilon  & length(seleit)<maxit){
      th_vie=th_est

      #### Metodo de informacion de item isotonico  o KL
      if(metodoSel %in% c('InfoIso','KL')) {
        auxiliar = which.min(abs(seqTheta - th_est)) # Selecciona fila (theta para buscar en la matriz)
        itsel = colnames(matSelAux)[which.max(matSelAux[auxiliar,])] # Selecciona el nombre del items con max Info o KL
        matSelAux = matSelAux[,colnames(matSelAux) != itsel]
      }

      #### Metodo Shannon Entrophy
      if(metodoSel %in% c('ESH')) {
        veroPrev = vero
        ## Transforma verosimilitud en una matriz con la dimension de curvaNOPAR
        #veroPrev = matrix(veroPrev,ncol = ncol(curvaNOPAR),nrow = nrow(curvaNOPAR),byrow = F)
        veroPrev = matrix(veroPrev,ncol = 1,nrow = nrow(curvaNOPAR),byrow = F)
        ## Normaliza para obtener el pi_n (Ecuacion 1)
        #pi_n = veroPrev/matrix(colSums(veroPrev),ncol = ncol(veroPrev),nrow = nrow(veroPrev),byrow = T)
        pi_n = veroPrev/sum(veroPrev)

        ## Numerador de la ecuacion (2)
        pi_yj1_num = curvaNOPAR*pi_n
        ## Denominador de la ecuacion (2) (Es tambien lo que aparece en la ecuacion (3)) se devería de cancelar
        Pj1 = colSums(pi_yj1_num)
        ## Ecuacion (2)
        # pi_yj1 = pi_yj1_num/matrix(Pj1,ncol = ncol(pi_yj1_num),nrow = nrow(pi_yj1_num),byrow = T)
        pi_yj1 = pi_yj1_num/Pj1
        # Lo mismo para armar las ecuaciones 4 y 5
        pi_yj0_num = (1 - curvaNOPAR)*pi_n
        Pj0 = colSums(pi_yj0_num)
        # pi_yj0 = pi_yj0_num/matrix(Pj0,ncol = ncol(pi_yj0_num),nrow = nrow(pi_yj0_num),byrow = T)
        pi_yj0 = pi_yj0_num/Pj0

        # Terminos de la esperanza de shanon
        SHE0 = -colSums(log(pi_yj0^pi_yj0))*Pj0
        SHE1 = -colSums(log(pi_yj1^pi_yj1))*Pj1

        # Esperanza de la entropia de Shannon
        SHE = data.frame(matrix(SHE0 + SHE1,nrow = 1,ncol = ncol(curvaNOPAR)))
        colnames(SHE) = colnames(curvaNOPAR)

        if(length(seleit) !=0){
          SHE=SHE[,!colnames(SHE) %in% seleit]
        }
        itsel = colnames(SHE)[which.min(SHE)]
      }

      ### Seleccion Aleatoria
      if(metodoSel %in% c('Random')) {
        if(length(seleit) !=0){
          itEliegibles = colnames(curvaNOPAR)[!colnames(curvaNOPAR) %in% seleit]
        }
        itsel = sample(colnames(curvaNOPAR),1)
      }

      seleit = c(seleit,itsel) # Guarda el item
      ### Metodo: KL

      #simula respuesta
      a=parametros[itsel,1]
      b=parametros[itsel,2]
      Psuj=exp(a*(th-b))/(1+exp(a*(th-b)))
      raux=runif(1,0,1)
      if(raux<Psuj){respu=1}else{respu=0}
      resp=c(resp,respu)
      #calcula verosimilitud
      vero=vero*(curvaNOPAR[,itsel]^respu)*(1-curvaNOPAR[,itsel])^(1-respu)
      #plot(vero)
      #estimacion de theta
      if( length(unique(unlist(resp))) == 1 ){## solo aciertos o solo errores
        delta=1
        if(resp[1]==1){ th_est = 0.75+length(resp)*0.05
        if(th_est>1){th_est=0.99}
        } else { th_est = 0.25-length(resp)*0.05
        if(th_est<0){th_est=0.01}
        }
      } else { ## si hay aciertos y errores
        th_est=puntuni[which.is.max(vero)]
        if(length(resp)<minit){delta=1} else {delta=abs(th_est-th_vie)}
      }
    }
    estima[s,1]=th_est
    estima[s,2]=delta
    estima[s,3]=length(resp)
    estima[s,4]=pnorm(th)
  }
  return(estima)
}