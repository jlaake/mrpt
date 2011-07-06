\name{lnl.removal}
\alias{lnl.removal}
\title{Mark-Recapture Distance Sampling: log-likelihood for removal configuration}
\description{Computes negative log-likelihood of removal observer configuration for double-observer point/line sampling}
\usage{
lnl.removal(par, x1, x2, nfreq, models, width, cutp, debug = FALSE, indep = FALSE, 
                PI = TRUE, use.offset = TRUE, posdep = FALSE, point=TRUE, mronly=FALSE)
}
\arguments{
  \item{par}{values for parameters}
  \item{x1}{dataframe of double observer data for observer 1}
  \item{x2}{dataframe of double observer data for observer 2}
  \item{nfreq}{frequency of observations for each row in x}
  \item{models}{list of formulae for p(y) (p.formula) and delta(y) (delta.formula)}
  \item{width}{radius of point transect or half-width of line transect}
  \item{cutp}{cut points for distances}
  \item{debug}{if TRUE, output given during iteration}
  \item{indep}{if TRUE, full independence assumed}
  \item{PI}{if TRUE, point independence assumed; if !indep & !PI -- limiting independence }
  \item{use.offset}{should always be default of TRUE}
  \item{posdep}{if TRUE forces positive dependence}
  \item{point}{if TRUE for point counts; otherwise, line transect}
  \item{mronly}{if TRUE only uses mark-recapture portion of the likelihood}
}
\details{
The likelihood that is computed depends on the value of \code{mronly} and \code{cutp}.
If \code{mronly=TRUE} then it returns the value of L_omega specified in the unbinned section of
Laake and Collier (in prep) for either binned (\code{cutp}) or unbinned data.  Binned data are handled 
with either distance being a factor variable (\code{p.formula=~factor(distance)})
or by using the \code{distance} value given (presumably the bin midpoint) as a numeric predictor(\code{p.formula=~distance}).
The L_omega likelihood is equivalent to the likelihood for the Huggins estimator in MARK for 2 occasions.  With
the removal model, it is equivalent to setting the recapture probability c=1 to achieve the removal configuration.                
The specification of pi(y) is irrelevant for this likelihood, so it is the same for points or lines.
               
If \code{mronly=FALSE}, then the likelihood is composed of L_y and L_omega and the
equations in Laake and Collier (in prep) in the unbinned and binned sections are used depending on whether the
data are binned or not. The value of \code{point} determines the value of pi(y) which is 1/width for
lines and 2y/w^2 for points.  Because the calculations involve integrals, the calculation of the log-likelihood 
takes longer than with \code{mronly=TRUE}. Note that L_y and L_omega are not calculated separately in this
case because the product is simpler because p.(y_i) in the denominator of L_y cancels with the same value
in the numerator of L_omega.

The code for this function is rather simple because most of the computation is in \code{\link{p.removal}} for \code{mronly=FALSE}
and in \code{\link{p.removal.mr}} for \code{mronly=TRUE}.
}
\value{negative log-likelihood for the defined model}
\author{Jeff Laake}
\seealso{ \code{\link{p.removal.mr}}, \code{\link{p.removal}}, \code{\link{fit.removal}} }



