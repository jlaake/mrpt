\name{fit.removal}
\alias{fit.removal}
\title{Mark-Recapture Distance Sampling: Removal Model Fitting}
\description{Fits removal model with binned or unbinned mark-recapture distance sampling data}
\usage{
fit.removal(x, par = NULL, p.formula = ~distance, delta.formula = ~-1+distance, 
  width = max(c(x$distance,cutp)), cutp = NULL, debug = FALSE, indep = NULL, PI = NULL,  
  method = "Nelder-Mead", posdep = NULL, hessian = TRUE, point=TRUE, mronly=FALSE)
}
\arguments{
  \item{x}{dataframe of double observer data with fields distance, Obs1, and Obs2 defined below}
  \item{par}{initial values for parameters}
  \item{p.formula}{formula for probability of detection p(y)}
  \item{delta.formula}{formula for dependence function delta(y)}
  \item{width}{radius of point transect or half-width of line transect}
  \item{cutp}{cut points for binned distances}
  \item{debug}{if TRUE, output given during iteration}
  \item{indep}{if TRUE, full independence assumed}
  \item{PI}{if TRUE, point independence assumed; if not indep and not PI -- limiting independence }
  \item{method}{optimization method used in optim; default is Nelder-Mead}
  \item{posdep}{if TRUE forces positive dependence}
  \item{hessian}{if TRUE, computes and returns hessian for fit}
  \item{point}{if TRUE for point counts; otherwise, line transect}
  \item{mronly}{if TRUE, uses only mark-recapture portion of the likelihood (ala Huggins model)}
  
}
\details{Uses \code{\link{lnl.removal}} and \code{\link{p.removal}} with optim to fit mark-recapture distance sampling model
with removal configuration of observers for point or line transect.  Model is specified by \code{p.formula} for detection function p(y) and
by \code{delta.formula} for the dependence function delta(y). The independence assumption is determined by arguments \code{indep} and \code{PI}.  If \code{indep==TRUE}
then full independence is assumed at all y which means that delta(y)=1 for all y.  If \code{PI==TRUE} then point 
independence is assumed which means that delta(0)=1 but delta(y) for all other y can be non-unity.  If \code{posdep==TRUE} (the default if \code{PI==TRUE})
then it forces positive dependence which means that delta(y)>=1 for all y.  If \code{mronly=TRUE} then \code{indep=TRUE} and it only uses the mark-recapture
portion of the likelihood. If \code{PI=FALSE} and \code{indep=FALSE} then it will use the limiting independence approach 
but this has not been tested with this code. 

The dataframe \code{x} contains a record for each sighting for each observer with fields distance (distance value), observer=1 (primary) or 2(secondary) and detected=0 or 1
to indicate if the observer missed(0) or detected (1) the animal(object). For all secondary observer records, detected=1. If differences due to personnel is going to be examined where 
persons rotated between primary and secondary, then you should have a field called pair which defines the ordering of the persons for each point.  For example,
if I included a field called person with values A and B, then pair should be defined as AB when A is primary and BA when B was primary.  
The pairing doesn't have to be of two persons and instead could be characteristics of the persons like experience. But a variable named pair
must be given with the pairing of the characteristics. See \code{\link{GCWA2008}} for examples. 

If the distances are in bins, distance should specify a value in the bin and then these are binned into
intervals with the \code{cutp} argument which is either a vector of breaks (c(0,50,100)) or specified as a matrix matrix(c(0,50,50,100),nrow=2,ncol=2,byrow=TRUE).  The
former approach is easiest. Rather than specify \code{cutp} it is also possible to define fields in \code{x} named \code{dist.begin} and 
\code{dist.end} that are the beginning and ending end point for the distance interval.  This allows the intervals to vary (eg. varying altitudes
with angle bins in aerial survey) and also allows inclusion of observations with exact distances \code{dist.begin==dist.end} (e.g., when photographic 
surveys are also part of an aerial survey).  Note that this latter aspect has not yet been tested.

The arguments \code{p.formula} and \code{delta.formula} are R formulas like p.formula=~distance.  When PI is chosen, \code{delta.formula} should
not contain an intercept because delta(0)=1.  To specify this use delta.formula=~-1+distance.  If you want to use another covariate to describe 
dependence then it will have to be as an interaction with distance (~-1+sex:distance); otherwise it will specify an non-zero intercepts.  Non-zero 
intercepts are only allowed with limiting independence which has not been tested with this code.  Also note that limiting independence does not 
always yield reasonable results.  See Buckland et al. (2010) for a better description of limiting independence.

The argument \code{par} is a vector of parameters that are ordered with parameters for \code{p.formula}
followed by those for \code{delta.formula}.  For example, for the default formulas of \code{p.formula=~distance} and \code{delta.formula=~-1+distance},
there are 2 parameters for p.formula (intercept and slope for distance) and 1 for delta.formula for the slope for distance.  These are ordered in 
par as the first 2 for p.formula and the third for delta.formula.   Initial values for \code{par} are computed if this argument is not specified.
Initial values can be specified either as a vector of values of the correct length or by specifying a model that was previously run.  It will
match up the parameters that are in common and will use 0 for any in the new model that were not in the model specified with \code{par}.
  
}
\value{
  \item{fit}{optim results}
  \item{beta}{vector of parameters for p(y)}
  \item{gamma}{vector of parameters for delta(y)}
  \item{AICc}{model selection value}
  \item{nfreq}{frequency of objects for each value in data}
  \item{alldata}{dataframe x specified in function call}
  \item{data}{unique set of distance data used for fitting}
  \item{dm}{design matrices specified by the models; xmat1 for observer 1, xmat2 for observer 2 and dmat for delta}
  \item{models}{list containing p.formula and delta.formula used to fit model}
  \item{control}{list of arguments used to control fit of model like width, cutp, etc}
  \item{p}{fitted probabilities for data}
  \item{Nhat}{computed abundance in the sampled areas}
}
\author{Jeff Laake}
\references{
Buckland, S. T., J. Laake and D. Borchers (2010). "Double observer line transect methods: levels of independence." Biometrics 66: 169-177.

Laake, J. L. and B. A. Collier (in prep). "Mark-Recapture Distance Sampling for Points"
}
\seealso{ \code{\link{p.removal.mr}}, \code{\link{p.removal}}, \code{\link{lnl.removal}} }
\examples{
# Simulated data
set.seed(42)
xx=gen.points(n=300,w=100,sigma1=-.03,detfct="logistic",p0_1=0.8,type="removal",gamma=log(0.03),posdep=TRUE)
mod=fit.removal(xx,width=100,cutp=c(0,25,50,75,100),PI=TRUE)
mod
plot(mod)
}

