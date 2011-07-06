\name{logistic.rem}
\alias{logistic.rem}
\alias{rem.p01}
\alias{rem.p11}
\alias{rem.pdot}
\alias{rem.pc1}
\alias{logistic.pt.rem01}
\alias{logistic.pt.rem11}
\alias{logistic.pt.remdot}
\alias{logistic.lt.rem01}
\alias{logistic.lt.rem11}
\alias{logistic.lt.remdot}
\title{Mark-Recapture Distance Sampling: Logistic Detection Probabilities and Integral Calculations}
\description{Computes probability values for capture histories 01,11 .1 (pdot) for removal and 10,01,11 .. (pdot) for io,
 and integrals across pi(y) for point or line.  The calculations use p(y) and delta(y) (if PI=TRUE) at a vector of
distances y to be used for calculation or integration.}
\usage{
rem.p01(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
rem.p11(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset)
rem.pdot(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
rem.pc1(xmat1,xmat2,dmat=NULL,beta,gamma,indep,PI,use.offset,posdep)
logistic.pt.rem01(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.pt.rem11(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.pt.remdot(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep)
logistic.lt.rem01(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.lt.rem11(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
logistic.lt.remdot(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep)  
}
\arguments{
  \item{distance}{vector of distance values}
  \item{xmat1}{design matrix for p1(y)}
  \item{xmat2}{design matrix for p2(y)}
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
