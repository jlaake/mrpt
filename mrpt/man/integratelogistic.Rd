\name{integratelogistic}
\alias{integratelogistic}
\alias{integratelogistic.io}
\alias{integratelogistic.rem}
\alias{integratelogistic.single}
\alias{logisticbyx.io}
\alias{logisticbyx.single}
\alias{logisticbyx.rem}
\title{Logistic integration for double observer models}
\description{Computes integral (numerically) over x from lower to upper of a logistic detection function. logisticbyx.rem or .io 
    treats logistic as a function of covariate z; for a given z it computes function fx with those covariate values at a range of distances.}
\usage{
integratelogistic.rem(fx, x1, x2, models, beta, gamma, width, lower = 0, upper = width, indep, PI, use.offset, posdep)
integratelogistic.io(fx,x1, x2, models, beta, gamma, width, lower = 0, upper = width, indep, PI, use.offset, posdep)
integratelogistic.single(fx, x, models, beta, width, lower=0,upper=width,p0)
logisticbyx.rem(distance, x1, x2, models, beta, gamma, width, indep, PI, use.offset, posdep, fx)
logisticbyx.io(distance, x1, x2, models, beta, gamma, width, indep, PI, use.offset, posdep, fx)
logisticbyx.single(distance,x,models,beta,width,fx,p0)
}
\arguments{
  \item{x}{covariate data for single observer model with multiple distance values for integration}
  \item{x1}{covariate data for observer 1 for io/removal model with multiple distance values for integration}
  \item{x2}{covariate data for observer 2 for io/removal model with multiple distance values for integration}
  \item{distance}{vector of distances specified by integration}
  \item{models}{list containing p.formula (formula for p(y)) and delta.formula (formula for delta(y))}
  \item{beta}{coefficients for p(y)}
  \item{gamma}{coefficients for delta(y)}
  \item{width}{radius of point count or half-width of line transect}
  \item{lower}{lower bound for integration}
  \item{upper}{upper bound for integration}
  \item{indep}{if TRUE, full independence assumed}
  \item{PI}{if TRUE, point independence assumed; if !indep & !PI -- limiting independence }
  \item{use.offset}{should always be default of TRUE}
  \item{fx}{function to be integrated}
  \item{p0}{probablity of detection at 0 distance p(y=0)}
  \item{posdep}{if TRUE, enforce positive dependence}
}
\details{
This is a set of functions that accomodate the integration of a logistic function for single observers (.single) and
removal (.rem) and independent observer (.io) configurations with double observers.  They are used to 
integrate various functions depending on the capture history and that value is passed as the argument \code{fx}.
The functions \code{integratelogistic.x} where .x is .io,.rem or .single are simple functions that call the  R function \code{integrate} 
to integrate the function \code{logisticbyx.x} where .x is .io,.rem or .single.  These functions in turn
call the specified \code{fx}. The set of \code{logisticbyx.x} functions handle the vectorized integration in which
integrate passes multiple distance values to be computed. But the integration is for a single set of values for the other covariates.  Those
functions use \code{expand.grid} to make a complete design matrix with all the different distance values with the 
same values of the covariates. Once that is done then it calls the function specified in \code{fx} with the 
appropriate design matrices and they are passed back to integrate.  

The functions used for \code{fx} are defined in \code{\link{logistic.rem}}, \code{\link{logistic.io}}, and \code{\link{logistic.single}}.
}
\value{The functions \code{integratelogistic.x} return a single value for the integral. The functions \code{logisticbyx.x} return
a vector of values for the function at the specified distances.}
\author{Jeff Laake}
\seealso{\code{\link{logistic.rem}}, \code{\link{logistic.io}}, \code{\link{logistic.single}}}
