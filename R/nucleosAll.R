
#' Condunto de nucleos unidimensionales y multidimensionelaes
#'
#' Los nucleos multidimensionales son: Normal multiplicativo (multinorm), Nucleo Epanechnikov esferico (epaesf),
#'  Normal multivariada (multinorm)
#'
#' @param x vector al que se le aplica el nucleo
#' @param h amplitud de la ventana.
#' @param th puntos donde se tiene valores de la funcion.
#' @param sigma matriz de varianza y covarianza para el caso de el nucleo *multinorm*
#' @return un vector de largo *x* con el resultado de la aplicacion del nucleo
#' @export
#' @examples
#' set.seed(1234)
#' xeval = seq(0,1,0.01)
#' h = 0.05
#' th = runif(1000)
#' res = normal(x = xeval,h = h,th = th)



normal=function(x,h,th){
  res=numeric(length(x))
  res=(1/sqrt(pi*2))*exp(-(((x-th)/h)^2)/2)
  return(list(res=res))
}

#Nucleo Epanechnikov

epa=function(x){
  res=numeric(length(x))
  res=(3/4)*(1-((x-th)/h)^2)*(abs(((x-th)/h))<=1)
  return(list(res=res))
}


acuepa=function(x){
  inte=(3/4)*(((x-th)/h)-(((x-th)/h)^3)/3+2/3)*(abs(((x-th)/h))<=1)
  return(list(inte))
}



#Nucleo uniforme

forme=function(x){
  res=numeric(length(x))
  res=(1/2)*(abs(((x-th)/h))<=1)
  return(list(res=res))
}



acunif=function(x){
  inte=(0.5+(((x-th)/h)/2))*(abs(((x-th)/h))<=1)
  return(list(inte))
}


#Nucleo Triangular

trian=function(x){
  res=numeric(length(x))
  res=(3/4)*(1-abs(((x-th)/h)))*(abs(((x-th)/h))<=1)
  return(list(res=res))
}



acutrian=function(x){
  if(((x-th)/h)>0){
    inte=((((x-th)/h)^2)/2+((x-th)/h)+0.5)*(abs(((x-th)/h))<=1)
  }else{
    inte=(0.5+((x-th)/h)-(((x-th)/h)^2)/2)*(abs(((x-th)/h))<=1)
  }
  return(list(inte))
}


#Nucleo Biweight

Biwei=function(x){
  res=numeric(length(((x-th)/h)))
  res=(15/16)*(1-((x-th)/h)^2)^2*(abs(((x-th)/h))<=1)
  return(list(res=res))
}


acubiw=function(x){
  inte=(15/16)*(((x-th)/h)-(2/3)*((x-th)/h)^3+(1/5)*((x-th)/h)^5+(8/15))*(abs(((x-th)/h))<=1)
  return(list(inte))
}



#Nucleo Triweight

Triwei=function(x){
  res=numeric(length(((x-th)/h)))
  res=(35/32)*(1-((x-th)/h)^2)^3*(abs(((x-th)/h))<=1)
  return(list(res=res))
}


acutri=function(x){
  inte=(35/32)*(((x-th)/h)-((x-th)/h)^3+(3/5)*((x-th)/h)^5-(((x-th)/h)^7)/7+(16/35))*(abs(((x-th)/h))<=1)
  return(list(inte))
}


#Nucleo Tricubo

Trcubo=function(x){
  res=numeric(length(((x-th)/h)))
  res=(70/81)*(1-abs(((x-th)/h))^3)^3*(abs(((x-th)/h))<=1)
  return(list(res=res))
}


acucub=function(x){
  if(((x-th)/h)>0){
    inte=(70/81)*(((x-th)/h)+3*(((x-th)/h)^4)/4+3*(((x-th)/h)^7)/7+(((x-th)/h)^10)/10+(81/140))*(abs(((x-th)/h))<=1)
  }else{
    inte=(70/81)*(((x-th)/h)-3*(((x-th)/h)^4)/4+3*(((x-th)/h)^7)/7-(((x-th)/h)^10)/10+(81/140))*(abs(((x-th)/h))<=1)
  }
  return(list(inte))
}


#Nucleo Coseno

Cos=function(x){
  res=numeric(length(((x-th)/h)))
  res=(pi/4)*cos((pi/2)*((x-th)/h))*(abs(((x-th)/h))<=1)
  return(list(res=res))
}


acucos=function(x){
  inte=(1/2)*(sin((pi/2)*((x-th)/h))+1)*(abs(((x-th)/h))<=1)
  return(list(inte))
}



#Nucleo Logistico

Logis=function(x){

  res=numeric(length(((x-th)/h)))
  res=1/(exp(((x-th)/h))+2+exp(-((x-th)/h)))
  return(list(res=res))
}


aculogis=function(x){
  inte=(1/2)*(sinh(((x-th)/h))/(1+cosh(((x-th)/h)))+1)
  return(list(inte))
}



#Nucleo Silverman

Silver=function(x){
  res=numeric(length(((x-th)/h)))
  res=(1/2)*exp(-abs(((x-th)/h))/sqrt(2))*sin((abs(((x-th)/h))/sqrt(2))+(pi/4))
  return(list(res=res))
}


#Nucleo Sigmoideo

Sigmo=function(x){
  res=numeric(length(((x-th)/h)))
  res=(2/pi)*(1/(exp(((x-th)/h))+exp(-((x-th)/h))))
  return(list(res=res))
}


acusig=function(x){
  inte=(2/pi)*(atan(exp(((x-th)/h))))
  return(list(inte))
}


#Nucleo normal multivariado

multinorm=function(x,sigma){

  #library(MASS)
  x=as.matrix(x)
  res=numeric(length=(dim(x)[1]))
  n=ncol(sigma)
  for(i in 1:length(res)){
    res[i]=exp((t(x[i,])%*%ginv(sigma)%*%x[i,])/-2)/(2*pi)^(n/2)*det(sigma)^(1/2)
  }
  #indicep=abs(res)<=3000
  #return(list(res=res,indicep=indicep))
  return(list(res=res))
}



#Nucleo Epanechnikov esferico

epaesf=function(x){

  x=as.matrix(x)
  res=numeric(length=(dim(x)[1]))
  for(i in 1:dim(x)[1]){
    res[i]=(1-t(x[i,]) %*% x[i,])*(t(x[i,]) %*% x[i,] <=1)
  }
  return(list(res=res))
}


#Normal multiplicativo
normmult=function(x){
  res=rep(1,dim(x)[2])

  for(l in 1:dim(x)[2]){
    res=res*norm(x[,l])
  }
  return(list(res=res))
}
