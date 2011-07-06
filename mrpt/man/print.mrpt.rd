\name{print.mrpt}
\alias{print.mrpt}
\alias{plot.mrpt}
\title{Print and plots for for mrpt models}
\description{Prints summaries and plots of fitted models}
\usage{
\method{print}{mrpt}(x,...)
\method{plot}{mrpt}(x,distance=NULL,titles="",single=TRUE,...)
}
\arguments{
  \item{x}{fitted model}
  \item{distance}{vector of distance values for plotting}
  \item{titles}{main titles for each plot}
  \item{single}{if TRUE plots detection function for single observer and for combined observers}
  \item{...}{additional arguments included for S3 compliance. Not used at present}
}
\details{
\code{print} prints summary of fitted model providing characteristics of model, AICc, Nhat and se, and parameter estimates and standard errors.
\code{plot} provides various plots for examining fitted model against data.  For mark-recapture (MR) only model
it only prints the conditional detection probability for the primary observer against a histogram which are the
fitted values if distance is treated as a step function. For MRDS models, several plots are given on a single page. Working from
top to bottom and left to right, the first plot is the fitted probability density function (pdf) of the observed distances
from both observers and a histogram scaled to be a pdf. The second is the fitted detection function for both observers 
combined (p.(y)) and a scaled histogram of the data.  If \code{single=TRUE}, then the third plot is the same thing but for the
primary observer only (p1(y)).  The third or fourth (if \code{single=TRUE}), is the conditional detection function which is
delta(y)*p1(y).  The histogram for the plot is based on a fitted model based on MR only. The final plot is the fitted
or assumed delta(y) function.
}
\author{Jeff Laake}

