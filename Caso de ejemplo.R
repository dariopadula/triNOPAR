library(dplyr)
library(mirt)
library(psych)
library(MASS)
library(np)

#1)Genera el dataframe de itemes
bancoIT = genBancoDF(nomFijas = c('NombreIt','Modelo'),nparam = 15)

#2) Genera el banco de items
set.seed(2021)
m2pl=matrix(c(0.2,1,-3,3),byrow = TRUE,ncol=2,nrow=2)
bancoIT=genitmodelo(100,"2PL",m2pl,bancoIT)


#3)Diseño del banco para TAI (aqui entrarian los experimentos)
#Experimento 1= genero un banco 25 itemes 2PL, estimo 2PL
diseno_exp1=matrix(c(25,"2PL"),ncol=2,byrow = FALSE)
bancTeor = bancoTAI(bancoIT,diseno_exp1)

#4) Generacion de habilidades reales teoricas-Normales
thetas=creathetas("Normal",n=1000,mean=0,sd=1)

#5)Se simula matriz de respuestas dado habilidades reales
respu=GenResp(bancTeor,thetas)

#6)Estimacion parametrica de los parametros de las ICC
# Estima los parametros asumiendo unidimensionalidad y que los items son todos 2PL
modelo=1
tipit<-c(rep("2PL",25))
paramsEst = estpar(respu,modelo,tipit)

#7)Generación de funciones gs
itdim1=colnames(respu)[1:25]
itemsgs=list(dim1=itdim1)
pesos=list(dim1=rep(1,25))
respG=estFunG_Simple(respu,itemsgs,pesos)
respG2 = estFunG(resp = respu,rotacion = c('oblimin'),corr = c('tet'),ndims = 1)# AF para habilidades

#8)Estimación NP de habilidades
thetaest = estthetaNP(scores = respG,Dth = qnorm)
thetaest2 = estthetaNP(scores = respG2,Dth = qnorm)

#9)Estimación ICC NP
#no estaba definida la grilla y la defini asi siguiendo lo mismo de icciso
nt=1000
grilla=c(1:nt/nt)
grilla[length(grilla)]
grilla[length(grilla)]=0.9999
ICCNP_todos = estRegNoPar(items = 1,h = 0.01,th_use = 'pcg',test = thetaest,puntos = grilla,nucleo=epa,sigma = 1)
ICCNP_todos2 = estRegNoPar(items = c(1:25),h = 0.01,th_use = 'pcg',test = thetaest2,puntos = grilla,nucleo=epa,sigma = 1)

#10)Estimación de ventanas
h=ventana1D(items = 1,th_use = 'dtg',test = thetaest,nucleodes="gaussian",1000)#ventana para los el primer item
ht=ventana1D(items = c(1:25),th_use = 'dtg',test = thetaest,nucleodes="gaussian",1000) #ventana para los 25 items

#11)ICC_param en grilla
a = 1
b = 0
c = 0.1
grilla = dnorm(seq(0,1,0.01))
p = iccFun(a,b,c,grilla)

#12)ICC_ISO

puntosiso=puntosicc=ICCNP_todos[[3]]
length(puntosiso)
puntosiso[length(puntosiso)]
plot(ICCNP_todos[[3]],icc1[,1],type='l',col='black',ylim=c(0,1),xlab='grilla',ylab='prob',main='ICCNP-Ramsay-IT1',lwd=2)
plot(ICCNP_todos[[3]],icc1[,2],type='l',col='black',ylim=c(0,1),xlab='grilla',ylab='prob',main='ICCNP-Ramsay-IT2',lwd=2)

icc1=ICCNP_todos[[1]]
icciso1 = icciso(icc1,hd,theta,nt,puntosiso,nucleod) #hd no esta definido

nt=1000
t=c(1:nt/nt)
t[length(t)]=0.9999
auxpu=c(0,puntosicc,1)
auxic=c(icc1[1],icc1[,1],1)
fes=numeric(length(t))#vector de 0
fi1<-approx(auxpu,auxic,t)$y

length(auxpu)
length(auxic)
length(t)

#13)
