% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parInfoFun.R
\name{parInfoFun}
\alias{parInfoFun}
\title{Para una grilla de puntos en el intervalo (0,1) calcula la funcion de informacion del item}
\usage{
parInfoFun(a, b, c, grillaUni, trnfFun = qnorm)
}
\arguments{
\item{a}{parametro de discriminacion del la ICC}

\item{b}{parametro de dificultad del la ICC}

\item{c}{parametro de pseudo azar de la ICC}

\item{grillaUni}{valores de theta donde se evalua la funcion}
}
\value{
Retorna un vector con la funcion de informacion evaluada en los puntos de la grilla
}
\description{
Para una grilla de puntos en el intervalo (0,1) calcula la funcion de informacion del item
}
\examples{
\dontrun{
a = 1; b = 0; c = 0.1
grillaUni = grilla =seq(0,1,0.01)
infor = parInfoFun(a,b,c,grillaUni,trnfFun = qnorm)
}
}
