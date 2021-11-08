#' Script de ejemplo de uso de las funciones
#'
#'
library(devtools)
# library(np)
#'libary(here)
load_all()
#'library(dplyr)
#'library(mirt)
#'library(psych)
#'library(MASS)

#******************************************************************
#Primero generamos data frame vacio donde iran los diferentes
#items que componen nuestro banco
#Esto se debe inicializar una sola vez a menos que se quiera cambiar
#el banco completo

bancoITEMS<- genBancoDF(nparam=15)

#Generamos un banco de 300 items
set.seed(12345)
m1pl=matrix(c(1,1,-2.5,2.5),byrow=T,ncol=2,nrow=2)
bancoITEMS=genitmodelo(1000,"1PL",m1pl,bancoITEMS)
m2pl=matrix(c(0.5,2,-2.5,2.5),byrow=T,ncol=2,nrow=2)
bancoITEMS=genitmodelo(1000,"2PL",m2pl,bancoITEMS)
m3pl=matrix(c(0.5,2,-2.5,2.5,0.2,0.4),byrow=T,ncol=2,nrow=3)
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



#Estimacion NP
# Estimo el item 1 NP
# uso la distribucion teorica de theta y la uniforme


hT=ventana1D(items = 1,th_use = 'dtg',test = thetaest,nucleodes="gaussian",muestra=2000)
hU=ventana1D(items = 1,th_use = 'pcg',test = thetaest,nucleodes="epanechnikov",muestra=2000)
ICCNPT=estRegNoPar(items=1,h=hT,th_use = 'dtg',test=thetaest,puntos=seq(-3,3,0.01),nucleo=normal,sigma=1)
ICCNPU=estRegNoPar(items=1,h=hU,th_use = 'pcg',test=thetaest,puntos=seq(0,1,0.001),nucleo=epa,sigma=1)

#Grafica del item 1
item=1
D=1
a=bancTeor[item,"P1"]
b=bancTeor[item,"P2"]
c=bancTeor[item,"P3"]
if(is.na(c)) c=0
thepl=ICCNPT$puntos
Prob1 <- c + (1 - c) /(1 + exp(-D * a * (thepl - b)))

estimados=paramsEst$parametros$items
ae=estimados[item,"a"]
be=estimados[item,"b"]
ce=estimados[item,"g"]
Prob2 <- ce + (1 - ce) /(1 + exp(-D * ae * (thepl - be)))

plot(thepl,Prob1,col="red",xlim=c(-4,4),ylim=c(0,1))
points(thepl,Prob2,col="cyan")
#points(thepl,icc)
points(thepl,ICCNPT$NPICC[,item],col="blue")

#Graficas 2
thepl=qnorm(seq(0,1,0.001))
plot(thepl,ICCNPU$NPICC,col="blue",xlim=c(-4,4),ylim=c(0,1))
Prob1 <- c + (1 - c) /(1 + exp(-D * a * (thepl - b)))
points(thepl,Prob1,col="red")


# Estimacion is'otona del item 1


hdi= 0.9*length(ICCNPT$NPICC)^(-1/5)*min(sd(ICCNPT$NPICC),(quantile(ICCNPT$NPICC,prob=0.75)-quantile(ICCNPT$NPICC,prob=0.25))/1.364)


ISINPT=icciso(icc1=ICCNPT$NPICC,hd=hdi,thetaiso=seq(0,1,0.001),nt=200,puntosicc=pnorm(ICCNPT$puntos),nucleod=epa)




thepl=qnorm(seq(0,1,0.001))
plot(thepl,ICCNPU$NPICC,col="blue",xlim=c(-4,4),ylim=c(0,1))
Prob1 <- c + (1 - c) /(1 + exp(-D * a * (thepl - b)))
points(thepl,Prob1,col="red")
points(thepl,ISINPT$resfin,col = 'orange')

#A mano

icc1=ICCNPT$NPICC
hd=0.4
thetaiso=seq(0,1,0.001)
nt=1000
puntosicc=pnorm(ICCNPT$puntos)
nucleod=epa

#######################
#### yo dario pongo esto

