\name{logistic.io}
\alias{logistic.io}
\alias{io.p01}
\alias{io.p10}
\alias{io.p11}
\alias{io.pdot}
\alias{logistic.pt.io01}
\alias{logistic.pt.io10}
\alias{logistic.pt.io11}
\alias{logistic.pt.iodot}
\alias{logistic.lt.io01}
\alias{logistic.lt.io10}
\alias{logistic.lt.io11}
\alias{logistic.lt.iodot}
\title{Functions for Logistic Detection Probabilities and Integral Calculations for Removal Distance Sampling Data}
\description{Computes functional values for capture histories 01,11 .1 (pdot) and values of those probabilities and
pi(y) for point or line removal distance sampling.  The calculations use p(y) and delta(y) (if PI) at a vector of
distances y to be used for calculation or integration.}
\usage{
io.p01(xmat1, xmat2, dmat,beta,gamma,indep,PI,use.offset,posdep)
io.p10(xmat1, xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
io.p11(xmat1, xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
io.pdot(xmat1, xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
logistic.pt.io01(distance,xmat1, xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.pt.io10(distance,xmat1, xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.pt.io11(distance,xmat1, xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.pt.iodot(distance,xmat1, xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep)
logistic.lt.io01(distance,xmat1, xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.lt.io10(distance,xmat1, xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.lt.io11(distance,xmat1, xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.lt.iodot(distance,xmat1, xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep)  
}
\arguments{
  \item{distance}{vector of distance values}
  \item{xmat1}{design matrix for p(y) for observer 1}
  \item{xmat2}{design matrix for p(y) for observer 2}
  \item{dmat}{design matrix for delta(y)}
  \item{beta}{parameters for p(y)}
  \item{gamma}{parameters for delta(y)}
  \item{width}{radius of point transect or half-width of line transect}
  \item{indep}{if TRUE, full independence assumed}
  \item{PI}{if TRUE, point independence assumed; if !indep & !PI -- limiting independence }
  \item{use.offset}{should always be default of TRUE}
  \item{posdep}{if TRUE, enforce positive dependence}
}
\details{}
\value{}
\author{Jeff Laake}
