\name{fit.single}
\alias{fit.single}
\title{Distance Sampling: Single Observer Model Fitting}
\description{Fits single observer model with binned or unbinned distance sampling data}
\usage{
fit.single(x, par = NULL, p.formula = ~distance, width = 1, cutp = NULL, debug = FALSE, 
              method = "Nelder-Mead", hessian = TRUE, point=TRUE, type="hn",p0=NULL)
}
\arguments{
  \item{x}{dataframe of double observer data with fields distance, Obs1, and Obs2 defined below}
  \item{par}{initial values for parameters}
  \item{p.formula}{formula for probability of detection p(y)}
  \item{width}{radius of point transect or half-width of line transect}
  \item{cutp}{cut points for binned distances}
  \item{debug}{if TRUE, output given during iteration}
  \item{method}{optimization method used in optim; default is Nelder-Mead}
  \item{hessian}{if TRUE, computes and returns hessian for fit}
  \item{point}{if TRUE for point counts; otherwise, line transect}
  \item{type}{if "hn" fits half-normal; if "hr", hazard rate; otherwise a logistic}
  \item{p0}{p(y=0); if type="logistic", this specifies the intercept which cannot be fitted with single observer data} 
}
\details{Uses \code{\link{lnl.single}} and \code{\link{p.single}} with optim to fit distance sampling model with a half-normal (type="hn"),
hazard rate (type="hr") or logistic.  Model is specified by \code{p.formula} for detection function p(y).  For \code{type="hn"} or \code{"hr"}, \code{p.formula} 
specifies the scale function sigma and does not include distance.  Whereas, for \code{type="logistic"}, \code{p.formula} should include distance, but 
you also need to specify the value for p0.  As an example for the logistic, the model would be ~-1+distance and with the specified value of p0,
the functional value would be plogis(log(p0/(1-p0))+beta[1]*distance) where beta[1] is the slope parameter for distance.  
}
\value{
  \item{fit}{optim results}
  \item{AICc}{model selection value}
  \item{nfreq}{frequency of objects for each value in data}
  \item{data}{unique set of distance data used for fitting}
  \item{models}{list containing p.formula used to fit model}
  \item{control}{list of arguments used to control fit of model like width, cutp, etc}
}        
\author{Jeff Laake}
\examples{
# Simulated data
beta=-.05
n=300
w=100
xx=gen.points(n,w,sigma1=abs(1/beta),detfct="hn",type="single",epsilon=.1)
# use BFGS instead of Nelder-Mead because there is only one parameter and the latter is unreliable in that case
mod=fit.single(xx,width=100,cutp=c(0,25,50,75,100),par=log(abs(1/beta)),type="hn",p0=1,p.formula=~1,method="BFGS")
mod
# following should be close to 20 (1/.05) in range of about 19-21; will vary because of random selection of points
# by gen.points
exp(mod$fit$par)
}

