\name{fyhist.pt.plot}
\alias{fyhist.pt.plot}
\alias{fy.pt.plot}
\alias{gyhist.pt.plot}
\alias{gy.pt.plot}
\alias{gcyhist.pt.plot}
\alias{gcy.pt.plot}
\alias{delta.pt.plot}
\title{MRDS plotting functions}
\description{Functions called from plot.mrpt that produce plots of fitted models}
\usage{
fyhist.pt.plot(model,linesize=1)
fy.pt.plot(model,distance=NULL,x=NULL,pdot=NULL)
gyhist.pt.plot(model,linesize=1,single=TRUE)
gy.pt.plot(model,distance=NULL,x=NULL,pdot=NULL,single=TRUE)
gcyhist.pt.plot(pc,model,linesize=1)
gcy.pt.plot(model,distance=NULL,x=NULL,pdot=NULL,plot=TRUE)
delta.pt.plot(model,distance=NULL,x=NULL,pdot=NULL)
}
\arguments{
  \item{model}{fitted model}
  \item{distance}{vector of distance values for plotting}
  \item{linesize}{size for fitted line}
  \item{x}{unique values of double observer dataframe for fitting model}
  \item{single}{if TRUE plots detection function for single observer and for combined observers}
  \item{pdot}{vector of detection probability that it was seen by at least one observer}
  \item{plot}{if TRUE, produces plot or if FALSE, it returns the dataframe of computed values}
  \item{pc}{vector of conditional detection probabilities for histogram bar heights} 
}
\details{
\code{fyhist.pt.plot} produces the scaled histogram of distancesfor the probability density plot. \code{fy.pt.plot} produces the plotted
fitted probability density function. \code{gyhist.pt.plot} produces scaled histogram of distances for detection probability plot. 
\code{gy.pt.plot} produces fitted detection function line. \code{gcyhist.pt.plot} produces histogram showing fitted
conditional detection function using step function for distances. \code{gcy.pt.plot} produces fitted conditional detection
function. \code{delta.pt.plot} produces fitted or assumed shape of delta(y).
}
\author{Jeff Laake}

