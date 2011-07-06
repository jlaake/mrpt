\name{p.single}
\alias{p.single}
\alias{detfct.pt}
\alias{detfct.lt}
\title{Single Observer Detection Probabilities}
\description{Computes probabilities for single observer}
\usage{
p.single(par,x,models,width,cutp,point,type,p0)
detfct.pt(distance,xmat,beta,width,type="hn")
detfct.lt(distance,xmat,beta,width,type="hn")
}
\arguments{
  \item{par}{parameter values}
  \item{x}{data for observer 1}
  \item{models}{list of formulae for p(y) (p.formula) and delta(y) (delta.formula)}
  \item{width}{radius of point transect or line transect half-width}
  \item{cutp}{cut points for distances}
  \item{point}{if TRUE for point counts; otherwise, line transect}
  \item{type}{detection function string which is either "hn","hr", or "logistic" for half-normal, hazard-rate, or logistic respectively; for detfct.pt and detfct.lt, only "hn" or "hr"}
  \item{p0}{probablity of detection at 0 distance p(y=0); this is required if type="logistic" and is unused otherwise}  
  \item{distance}{vector of distance values}
  \item{xmat}{design matrix of values for scale function}
  \item{beta}{parameter vector for scale function}
}
\details{}
\value{}
\author{Jeff Laake}
