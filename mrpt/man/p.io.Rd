\name{p.io}
\alias{p.io}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mark Recapture Distance Sampling: Independent Configuration Detection Probabilities}
\description{Computes probabilities for capture histories (10,01,11) in independent observer configuration}
\usage{
p.io(par,x1,x2,models,width,cutp,indep,PI,use.offset,posdep,point)
}
\arguments{
  \item{par}{parameter values}
  \item{x1}{data for observer 1}
  \item{x2}{data for observer 2}
  \item{models}{list of formulae for p(y) (p.formula) and delta(y) (delta.formula)}
  \item{width}{radius of point transect or line transect half-width}
  \item{cutp}{cut points for distances}
  \item{indep}{if TRUE, full independence assumed}
  \item{PI}{if TRUE, point independence assumed; if !indep & !PI -- limiting independence }
  \item{use.offset}{should always be default of TRUE}
  \item{posdep}{if TRUE enforces positive dependence}
  \item{point}{if TRUE for point counts; otherwise, line transect}
}
\details{}
\value{}
\author{Jeff Laake}
