\name{logistic.single}
\alias{logistic.single}
\alias{s.1}
\alias{logistic.pt.s1}
\alias{logistic.lt.s1}
\title{Functions for Logistic Detection Probabilities and Integral Calculations for Single Observer}
\description{Computes logistic detection function values for fitting of a logistic detection function to single observer data.}
\usage{
s.1(xmat,beta,p0)
logistic.pt.s1(distance,xmat,beta,width,p0) 
logistic.lt.s1(distance,xmat,beta,width,p0) 
}
\arguments{
  \item{distance}{vector of distance values}
  \item{xmat}{design matrix for p(y)}
  \item{beta}{parameters for p(y)}
  \item{p0}{detection probability at 0; p(0))}
  \item{width}{radius of point transect or half-width of line transect}
}
\details{}
\value{}
\author{Jeff Laake}
