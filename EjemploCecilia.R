#' Script de ejemplo de uso de las funciones
#'
#'
library(devtools)
# library(np)
#'libary(here)
load_all()
library(dplyr)
#'library(mirt)
#'library(psych)
#'library(MASS)

#******************************************************************
#Actualizando 2022_sincro
#Primero generamos data frame vacio donde iran los diferentes
#items que componen nuestro banco
#Esto se debe inicializar una sola vez a menos que se quiera cambiar
#el banco completo

bancoITEMS<- genBancoDF(nparam=15)

#Generamos un banco de 300 items
set.seed(12345)
m1pl=matrix(c(1,1,-2.5,1.5),byrow=T,ncol=2,nrow=2)#Baje la dificultad
bancoITEMS=genitmodelo(1000,"1PL",m1pl,bancoITEMS)
m2pl=matrix(c(0.5,2,-2.5,1.5),byrow=T,ncol=2,nrow=2)
bancoITEMS=genitmodelo(1000,"2PL",m2pl,bancoITEMS)
m3pl=matrix(c(0.5,2,-2.5,1.5,0.2,0.4),byrow=T,ncol=2,nrow=3)
bancoITEMS=genitmodelo(1000,"3PL",m3pl,bancoITEMS)

#diseno=matrix(c(100,"2PL",100,"1PL",100,"3PL"),ncol=2,byrow = T)
diseno=matrix(c(200,"1PL"),ncol=2,byrow = T)
bancTeor = bancoTAI(bancoITEMS,diseno)
thetas=creathetas("Normal",n=10000,mean=0,sd=1)
respu=GenResp(banco = bancTeor,thetas)

# Estimacion paramétrica de los parametros en 1 dimension de bancTeor
modelo=1
#tipit<-c(rep("2PL",100),rep("Rasch",100),rep("3PL",100))
tipit<-c(rep("Rasch",200))
paramsEst = estpar(respu,modelo,tipit)

#Estimados
aest=paramsEst$parametros$items[,"a"]
best=paramsEst$parametros$items[,"b"]
cest=paramsEst$parametros$items[,"g"]

#' Calcula la correlacion y diferencia entre estimado y real en el parametro b
cor(best,bancTeor[,4])
difepar=best-bancTeor[,4]



# Habilidades
thest=paramsEst$habilidades
mean(thest)
sd(thest)

difeth=thest-thetas

# Estimacion no parametrica Caso 1 dimension

itdim1=colnames(respu)[1:200]
itemsgs=list(dim1=itdim1)
pesos=list(dim1=rep(1,200))
respG=estFunG_Simple(respu,itemsgs,pesos)

thetaest = estthetaNP(scores = respG,Dth = qnorm)

cor(thetaest[,"dtg1"],thetas)
mean(thetaest[,"dtg1"])
sd(thetaest[,"dtg1"])

difeth2=thetaest[,"dtg1"]-thetas
unido=cbind(thetas,thetaest[,"dtg1"],difeth2)


#Estimacion NP de las ICC con la ventana elegida de forma optima
## Lo hago de forma paralela
library(parallel)
library(data.table)


puntosNP = seq(0,1,0.001)
#puntosNP = seq(-3,3,0.01)
th_use = 'pcg'
#th_use = 'dtg'
#nucelo=normal #epa #aca la dif con el ej mario


funEstNoPar = function(ii) {
  hT=ventana1D(items = ii,th_use = th_use,test = thetaest,nucleodes="gaussian",muestra=2000)
  iccnp=estRegNoPar(items = ii,h=last(hT),th_use = th_use,
                    test=thetaest,puntos=puntosNP,
                    nucleo=epa,sigma=1)

  return(as.vector(iccnp[["NPICC"]]))
}



trials = 1:200

cl <- makeCluster(detectCores())
clusterExport(cl, c("thetaest","puntosNP","th_use","ventana1D","estRegNoPar","normal", "epa"))
clusterEvalQ(cl, {
  library(tidyverse)
  library(np)})

tini = Sys.time()
system.time({
  results <- parallel::parLapply(cl,trials,funEstNoPar)
})


stopCluster(cl)
tfin = Sys.time()
difftime(tfin,tini,units = 'secs')

#### data frame con curvas ICC estimadas con regresion no parametrica
iccnp_mat = do.call(cbind,results)
colnames(iccnp_mat) = names(thetaest)[trials]

#Graficos de ICCNP
thepl=qnorm(seq(0,1,0.001))
#thepl=ICCNPT$puntos
D=1
parametros = bancTeor[,c('NombreIt','Modelo',paste0('P',1:3))]
a<-filter(parametros, NombreIt=="IT488")[,"P1"]
b<-filter(parametros, NombreIt=="IT488")[,"P2"]
c<-filter(parametros, NombreIt=="IT488")[,"P3"]
c<-0
plot(thepl,iccnp_mat[,"IT4881PL"],col="blue",xlim=c(-4,4),ylim=c(0,1))
Prob1 <- c + (1 - c) /(1 + exp(-D * a * (thepl - b)))
points(thepl,Prob1,col="red")

###Caso de item con difficultad =0.5
a<-filter(parametros, NombreIt=="IT789")[,"P1"]
b<-filter(parametros, NombreIt=="IT789")[,"P2"]
c<-filter(parametros, NombreIt=="IT789")[,"P3"]
c<-0
plot(thepl,iccnp_mat[,"IT7891PL"],col="blue",xlim=c(-4,4),ylim=c(0,1))
Prob1 <- c + (1 - c) /(1 + exp(-D * a * (thepl - b)))
points(thepl,Prob1,col="red")

################################ ############
############################################
############################################
##### Estima curvas np isotonicas de forma paralela

puntosicc = puntosNP
if(sum(abs(puntosNP) > 1) > 0) puntosicc = pnorm(puntosNP)

thetaiso = puntosNP

funEstNoParIso = function(ii) {

  iccnp = iccnp_mat[,ii]

  hdi= 0.9*length(iccnp)^(-1/5)*min(sd(iccnp),(quantile(iccnp,prob=0.75)-quantile(iccnp,prob=0.25))/1.364)
  nopariso=icciso(icc1=iccnp,hd=hdi,
                  thetaiso = thetaiso,nt=200,
                  puntosicc = puntosicc,nucleod=epa)

  return(nopariso$resfin)
}



cl <- makeCluster(detectCores())
clusterExport(cl, c("iccnp_mat","puntosicc","icciso","normal","epa","thetaiso"))
clusterEvalQ(cl, {
  library(tidyverse)
  library(np)})

tini = Sys.time()
system.time({
  resultsISO <- parallel::parLapply(cl,trials,funEstNoParIso)
})


stopCluster(cl)

tfin = Sys.time()
difftime(tfin,tini,units = 'secs')


#### data frame con curvas ICC estimadas de forma isotonica
icciso_mat = do.call(cbind,resultsISO)
#colnames(icciso_mat) = names(iccnp_mat)
colnames(icciso_mat) =names(thetaest)[trials]
sum(is.na(icciso_mat[1,])) #Tenemos 52 NAs al inicio
sum(is.na(icciso_mat[1001,]))#Todos los itemes tienen NA en el ultimo punto de la grilla

#Graficos de ICCNPISO
D=1
parametros = bancTeor[,c('NombreIt','Modelo',paste0('P',1:3))]
a<-filter(parametros, NombreIt=="IT488")[,"P1"]
b<-filter(parametros, NombreIt=="IT488")[,"P2"]
c<-filter(parametros, NombreIt=="IT488")[,"P3"]
c<-0
plot(thepl,icciso_mat[,"IT4881PL"],col="blue",xlim=c(-4,4),ylim=c(0,1))
Prob1 <- c + (1 - c) /(1 + exp(-D * a * (thepl - b)))
points(thepl,Prob1,col="red")

###Caso de item con difficultad =0.5
a<-filter(parametros, NombreIt=="IT789")[,"P1"]
b<-filter(parametros, NombreIt=="IT789")[,"P2"]
c<-filter(parametros, NombreIt=="IT789")[,"P3"]
c<-0
plot(thepl,icciso_mat[,"IT7891PL"],col="blue",xlim=c(-4,4),ylim=c(0,1))
Prob1 <- c + (1 - c) /(1 + exp(-D * a * (thepl - b)))
points(thepl,Prob1,col="red")

############################################
############################################
############################################
####### Estima la funcion de informacion para las curvas parametricas y las isotonicas

#### Items parametricos
parEst = paramsEst[[1]]$items
colnames(parEst) = c('a','b','c','d')

infoFunPar = do.call(cbind,
                   sapply(trials, function(xx) {
                      a = parEst[xx,'a']
                      b = parEst[xx,'b']
                      c = parEst[xx,'c']

                      info = parInfoFun(a,b,c,
                                        grillaUni = puntosicc,
                                        trnfFun = qnorm)
              list(info)
}))

colnames(infoFunPar) = rownames(parEst)[trials]

sum(is.na(infoFunPar))



### CURVAS ISOTONICAS

infoFunIso = do.call(cbind,
                     sapply(trials, function(xx) {

                       icciso = icciso_mat[,xx]
                       info = derISO_Info(icciso,
                                          puntuni = thetaiso,
                                          nucleo = normal,
                                          hd = NULL)
                       list(info$Info)
                     }))

colnames(infoFunIso) = colnames(icciso_mat)[trials]


sum(is.na(infoFunIso))#arrastramos los 252 NA

#################################################
#################################################
#################################################
##### Calcula la funcion KL para los items parametris los NP y lo NPiso
## Parametrico
KLFunPar = kl_mat_par(paramsMIRT = parEst[trials,],
                      grillaEval = thetaiso,
                      distTrans = qnorm)

colnames(KLFunPar) = rownames(parEst)[trials]

sum(is.na(KLFunPar))#no hay NA
##############################
### ICC no par
KLFunNoPar = kl_mat_NOpar(iccNP_mat = iccnp_mat[,trials],sepGrilla = 0.001,entorno = 0.1)
sum(is.na(KLFunNoPar))#no hay NA
##############################
### ICC no par isotonica
KLFunNoParIso = kl_mat_NOpar(iccNP_mat = icciso_mat[,trials],sepGrilla = 0.001,entorno = 0.1)
sum(is.infinite(KLFunNoParIso))
sum(is.na(KLFunNoParIso))

#IT31PL caso de un item que da infinito
plot(thepl,icciso_mat[,"IT31PL"],col="blue",xlim=c(-4,4),ylim=c(0,1))


####################################################################
####### INSUMOS PARA EL TAI (prueba)

## Parametros verdaderos (adecuo el df para que se interprete bien en la funcion)
parametros = bancTeor[,c('NombreIt','Modelo',paste0('P',1:3))]
rownames(parametros) = apply(parametros[,c('NombreIt','Modelo')],1,function(xx) paste(xx,collapse = ''))
parametros$P3 = ifelse(is.na(parametros$P3),0,parametros$P3)



#sujtai = rnorm(100) (ya teniamos definidos los reales con creatheta)
set.seed(2021)
sujtai = sample(thetas, 100)
epsilon = 0.01
minit = 10
maxit = 20
curvaNOPAR = NULL#si es nulo juega lo parametrico

##Simulacion con ICC parametricas
res1 = TAIgeneric(sujtai,
                       epsilon,
                       minit,
                       maxit,
                       curvaNOPAR,
                       parametros,
                       parEst,
                       itemsSelec = c('InfoFun'),
                       matrizSelect = infoFunPar,
                       seqTheta = puntosNP)
res2=TAIgeneric(sujtai,
                epsilon,
                minit,
                maxit,
                curvaNOPAR,
                parametros,
                parEst,
                itemsSelec = c('KL'),
                matrizSelect = KLFunPar,
                seqTheta = puntosNP)
res3=TAIgeneric(sujtai,
                epsilon,
                minit,
                maxit,
                curvaNOPAR,
                parametros,
                parEst,
                itemsSelec = c('Random'),
                #matrizSelect = infoFunPar,
                seqTheta = puntosNP)

##Simulación con ICC-NP
res4=TAIgeneric(sujtai,
                epsilon,
                minit,
                maxit,
                curvaNOPAR = iccnp_mat,
                parametros,
                parEst,
                itemsSelec = c('KL'),
                matrizSelect = KLFunNoPar,
                seqTheta = puntosNP)

res5=TAIgeneric(sujtai,
                epsilon,
                minit,
                maxit,
                curvaNOPAR = iccnp_mat,
                parametros,
                parEst = NULL,
                itemsSelec = c('ESH'),
                #matrizSelect = infoFunPar,
                seqTheta = puntosNP)

res6=TAIgeneric(sujtai,
                epsilon,
                minit,
                maxit,
                curvaNOPAR = iccnp_mat,
                parametros,
                parEst,
                itemsSelec = c('Random'),
                #matrizSelect = KLFunNoPar,
                seqTheta = puntosNP)

##Simulación con ICCISO

res7=TAIgeneric(sujtai,
                epsilon,
                minit,
                maxit,
                curvaNOPAR = icciso_mat,
                parametros,
                parEst,
                itemsSelec = c('InfoFun'),
                matrizSelect = infoFunIso,
                seqTheta = puntosNP)

res8=TAIgeneric(sujtai,
                epsilon,
                minit,
                maxit,
                curvaNOPAR = icciso_mat,
                parametros,
                parEst,
                itemsSelec = c('KL'),
                matrizSelect = KLFunNoParIso,
                seqTheta = puntosNP)


res9=TAIgeneric(sujtai,
                epsilon,
                minit,
                maxit,
                curvaNOPAR = icciso_mat,
                parametros,
                parEst,
                itemsSelec = c('ESH'),
                #matrizSelect = KLFunNoParIso,
                seqTheta = puntosNP)

res10=TAIgeneric_Ceci(sujtai,
                 epsilon,
                 minit,
                 maxit,
                 curvaNOPAR = icciso_mat,
                 parametros,
                 parEst,
                 itemsSelec = c('Random'),
                 #matrizSelect = KLFunNoParIso,
                 seqTheta = puntosNP)



###################################################
#### Calculo de los errores
errYses1 = ERRYSES(simData = res1,grilla = seq(1:100)/100) #para-MI
errYses1[["raiz"]]
errYses2 = ERRYSES(simData = res2,grilla = seq(1:100)/100)#para-KL
errYses2[["raiz"]]
errYses3 = ERRYSES(simData = res3,grilla = seq(1:100)/100)#para-Random
errYses3[["raiz"]]
errYses4 = ERRYSES(simData = res4,grilla = seq(1:100)/100)#NP-KL
errYses4[["raiz"]]
errYses5 = ERRYSES(simData = res5,grilla = seq(1:100)/100)#NP-ESH
errYses5[["raiz"]]
errYses6 = ERRYSES(simData = res6,grilla = seq(1:100)/100)#NP-RANDOM
errYses6[["raiz"]]
errYses7 = ERRYSES(simData = res7,grilla = seq(1:100)/100)#ISO-MI***
errYses7[["raiz"]]
errYses8 = ERRYSES(simData = res8,grilla = seq(1:100)/100)#ISO-MI***
errYses8[["raiz"]]

#Actualizado al 28_01_22
