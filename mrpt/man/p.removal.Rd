\name{p.removal}
\alias{p.removal}
\alias{p.removal.mr}
\title{Mark-Recapture Distance Sampling: Removal Configuration Detection Probabilities}
\description{Computes probabilities for capture histories (01,11) in removal observer configuration}
\usage{
p.removal(par, x1, x2, models, width, cutp = NULL, indep = FALSE, PI = TRUE, use.offset = TRUE, posdep = FALSE, point=TRUE)
p.removal.mr(par, x1, x2, models)
}
\arguments{
  \item{par}{values for parameters}
  \item{x1}{dataframe of double observer data for observer 1}
  \item{x2}{dataframe of double observer data for observer 2}
  \item{models}{lsit of formulae for p(y) and delta(y)}
  \item{width}{radius of point transect}
  \item{cutp}{cut points for distances}
  \item{indep}{if TRUE, full independence assumed}
  \item{PI}{if TRUE, point independence assumed; if !indep & !PI -- limiting independence }
  \item{use.offset}{should always be default of TRUE}
  \item{posdep}{if TRUE forces positive dependence}
  \item{point}{if TRUE for point counts; otherwise, line transect}
}
\details{
\code{p.removal} computes probabilities for each capture history when \code{mronly=FALSE} in call to \code{\link{lnl.removal}} for a composite
mark-recapture distance sampling and \code{p.removal.mr} is used to compute them when \code{mronly=TRUE}.

\code{p.removal.mr} is quite simple as it calls \code{\link{rem.p01}}, \code{\link{rem.p11}}, \code{\link{rem.pdot}}
to compute the probabilities directly from the logistic detection function with the \code{beta} parameter
values (in this case all of \code{par}) and the design matrices specified by \code{models$p.formula} and applied to \code{x1} and \code{x2} for the
two observers.  Those functions allow specification of delta(y) but in this case delta(y)=1 for all y.

\code{p.removal} is much more complicated because it has to handle both binned and unbinned cases and 
it involves integration over pi(y) which differs for line and point sampling.  The integration over distance (y) 
is complicated by the fact that integration must be done for each unique set of covariate values used in 
\code{models$p.formula} for p(y) and that formula could specify an interaction with y and another covariate.
For binned and unbinned data, \code{pdot} is an integration (over 0 to width of p.(y)pi(y)). It is computed
for each unbinned observation or group of observations with similar covariate values and distance bin for binned data. This could be 
made more clever for unbinned data. The integration is accomplished with the function \code{\link{integratelogistic.rem}} 
which integrates the function \code{\link{logistic.pt.remdot}} for points or \code{\link{logistic.lt.remdot}} 
for lines. For unbinned data the only part of the log-likelihood is the computation of probability of each
capture history (01,11) times pi(y).  The functions to compute these values are
\code{\link{logistic.pt.rem01}} and \code{\link{logistic.pt.rem11}} for points and 
\code{\link{logistic.lt.rem01}} and \code{\link{logistic.lt.rem11}} for lines.  They only differ
based on calculation of pi(y) and could have been combined but I did not want a conditional statement
to have to be tested each time they are called when they are integrated.  For binned
data, each of the 2 functions (for 01 and 11) are integrated over the range of the distance bin
through a call to  \code{\link{integratelogistic.rem}}.  So for binned data there are 3 integrals for
each unique combination of distance bin and covariate values and for unbinned distances there is
only one integration but currently it is done for each observation separately.

From these calculations, a list of 3 vectors is constructed with values p01,p11 and pdot for the 3 quantities for
both functions. For \code{p.removal}, note that for unbinned data p01 and p11 are not strictly probabilities but are the probability of the
capture history times pi(y), so they are strictly probability densities since y is continuous.  For binned data
they are probabilities over the interval averaged based on the expected distribution of pi(y).
}
\value{list of probabilities
  \item{p01}{probability of 01 capture history}
  \item{p11}{probability of 11 capture history}
  \item{pdot}{probability of 01 or 11 capture history}
}
\author{Jeff Laake}
\seealso{\code{\link{lnl.removal}}}
