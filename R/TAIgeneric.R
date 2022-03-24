


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
#' @param parEst Data frame con los parametros estimados de los items. Las filas tienes que ser nombradas
#' con los nombres de los items.
#' @param itemsSelec Metodo utilizado para seleccionar los items:
#'
#' * *InfoFun*: pseudo informacion
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
#' \dontrun{
#' # No correr
#' }
#'


######## FUNCION GENERICA PARA APLICAR TAI
TAIgeneric <- function(sujtai,
                       epsilon,
                       minit,
                       maxit,
                       curvaNOPAR,
                       parametros,
                       parEst,
                       itemsSelec = c('InfoFun','KL','ESH','Random'),
                       matrizSelect = MINOPARAUX,
                       seqTheta) {

  estima=matrix(NA,nrow=length(sujtai),ncol=4)
  for(s in 1:length(sujtai)){
    if((s %% 1000) == 0) print(s)

    metodoSel = itemsSelec[1] # TOma el primero por defecto


    ### Nombre de los items
    itemsNoms = rownames(parametros)

    resp=numeric(0)
    vero=1
    th=sujtai[s]
    th_est=runif(1,0.3,0.7)
    delta=1
    control=0
    seleit=NULL

    if(metodoSel %in% c('InfoFun','KL')) matSelAux = matrizSelect
    #seleccion del item
    while(delta > epsilon  & length(seleit)<maxit){
      th_vie=th_est

      #### Metodo de informacion de item isotonico  o KL
      if(metodoSel %in% c('InfoFun','KL')) {
        auxiliar = which.min(abs(seqTheta - th_est)) # Selecciona fila (theta para buscar en la matriz)
        itsel = colnames(matSelAux)[which.max(matSelAux[auxiliar,])] # Selecciona el nombre del items con max Info o KL
        matSelAux = matSelAux[,colnames(matSelAux) != itsel]
      }

      #### Metodo Shannon Entrophy
      if(metodoSel %in% c('ESH')) {
        veroPrev = vero
        ## Transforma verosimilitud en una matriz con la dimension de curvaNOPAR
        #veroPrev = matrix(veroPrev,ncol = 1,nrow = nrow(curvaNOPAR),byrow = F)
        veroPrev = matrix(veroPrev,ncol = ncol(curvaNOPAR),nrow = nrow(curvaNOPAR),byrow = F)
        ## Normaliza para obtener el pi_n (Ecuacion 1)

        # pi_n = matrix(veroPrev/sum(veroPrev),nrow = nrow(curvaNOPAR),ncol = 1,byrow = F)
        pi_n = matrix(veroPrev/sum(veroPrev),ncol = ncol(curvaNOPAR),nrow = nrow(curvaNOPAR),byrow = F)

        ## Numerador de la ecuacion (2)
        pi_yj1_num = curvaNOPAR*pi_n
        ## Denominador de la ecuacion (2) (Es tambien lo que aparece en la ecuacion (3)) se devería de cancelar
        Pj1 = colSums(pi_yj1_num)
        ## Ecuacion (2)
        pi_yj1 = pi_yj1_num/Pj1
        # Lo mismo para armar las ecuaciones 4 y 5
        pi_yj0_num = (1 - curvaNOPAR)*pi_n
        Pj0 = colSums(pi_yj0_num)
        pi_yj0 = pi_yj0_num/Pj0

        # Terminos de la esperanza de shanon
        SHE0 = -colSums(log(pi_yj0^pi_yj0))*Pj0
        SHE1 = -colSums(log(pi_yj1^pi_yj1))*Pj1

        # Esperanza de la entropia de Shannon
        SHE = data.frame(matrix(SHE0 + SHE1,nrow = 1,ncol = ncol(curvaNOPAR)))
        colnames(SHE) = itemsNoms

        if(length(seleit) !=0){
          SHE=SHE[,!colnames(SHE) %in% seleit]
        }
        itsel = colnames(SHE)[which.min(SHE)]
      }

      ### Seleccion Aleatoria
      if(metodoSel %in% c('Random')) {
        if(length(seleit) !=0){
          itEliegibles = itemsNoms[!itemsNoms %in% seleit]
          itsel = sample(itEliegibles,1)
        }
        itsel = sample(itemsNoms,1)
      }

      seleit = c(seleit,itsel) # Guarda el item


      #simula respuesta
      a=parametros[itsel,'P1']
      b=parametros[itsel,'P2']
      c=parametros[itsel,'P3']

      Psuj = iccFun(a,b,c,th)
      raux=runif(1,0,1)
      if(raux<Psuj){respu=1}else{respu=0}
      resp=c(resp,respu)

      # plot(seqTheta,matrizSelect[,itsel])
      # abline(v = th_est)

      # Calcula la curva dependiendo si es NP o Parametrica
      if(!is.null(curvaNOPAR)) {
        curva = curvaNOPAR[,itsel]
      } else {
        aEst = parEst[itsel,'a']
        bEst = parEst[itsel,'b']
        cEst = parEst[itsel,'c']
        curva = iccFun(aEst,bEst,cEst,grilla = qnorm(seqTheta))
      }
      #calcula verosimilitud
      vero=vero*(curva^respu)*(1-curva)^(1-respu)
      # plot(qnorm(seqTheta),curva)
      # print(itsel)
      #estimacion de theta
      if( length(unique(unlist(resp))) == 1 ){## solo aciertos o solo errores
        delta=1
        if(resp[1]==1){
          th_est = min(0.75+length(resp)*0.05,0.99)
        } else {
          th_est = max(0.25-length(resp)*0.05,0.01)
        }
      } else { ## si hay aciertos y errores
        th_est=seqTheta[which.max(vero)]
        if(length(resp)<minit){delta=1} else {delta=abs(th_est-th_vie)}
      }
      # print(th_est)
    }
    estima[s,1]=th_est
    estima[s,2]=delta
    estima[s,3]=length(resp)
    estima[s,4]=pnorm(th)
  }
  return(estima)
}
