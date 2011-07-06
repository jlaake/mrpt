\name{lnl.single}
\alias{lnl.single}
\title{Single Observer log-likelihood}
\description{Computes negative log-likelihood for single observer point/line sampling}
\usage{
lnl.single(par,xdf,xmat,nfreq,models,width,cutp,debug,point,type,p0)
}
\arguments{
  \item{par}{values for parameters}
  \item{xdf}{dataframe for observer}
  \item{xmat}{design matrix for dataframe}
  \item{nfreq}{frequency of observations for each row in x}
  \item{models}{list of formulae for p(y) (p.formula) and delta(y) (delta.formula)}
  \item{width}{radius of point transect or half-width of line transect}
  \item{cutp}{cut points for distances}
  \item{debug}{if TRUE, output given during iteration}
  \item{point}{if TRUE for point counts; otherwise, line transect}
  \item{type}{detection function string which is either "hn","hr", or "logistic" for half-normal, hazard-rate, or logistic respectively}
  \item{p0}{probablity of detection at 0 distance p(y=0); this is required if type="logistic" and is unused otherwise}  
}
\details{}
\value{}
\author{Jeff Laake}
