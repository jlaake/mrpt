\name{gen.points}
\alias{gen.points}
\title{Generates random sample for point counts}
\description{Generate simulated point transect data with optional capture history for 2 observers. No observer effect allowed p1=p2.}
\usage{gen.points(n,w,sigma1,sigma2=NULL,detfct="hn",p0_1=0.9999999,p0_2=NULL,type="single",gamma=NULL,epsilon=NULL,all=FALSE,posdep=FALSE, rotate=TRUE)}
\arguments{
  \item{n}{index to sample size (also depends on sigma,p0 etc)}
  \item{w}{radius of point transect}
  \item{sigma1}{for observer 1; sigma for half-normal detection function or sigma = distance slope for logistic}
  \item{sigma2}{for observer2; sigma for half-normal detection function or sigma = distance slope for logistic}
  \item{detfct}{if "hn" uses half-normal; otherwise uses logistic}
  \item{p0_1}{for observer 1;probability of detection at distance 0}
  \item{p0_2}{for observer 2;probability of detection at distance 0}
  \item{type}{single observer "single"; dependent dual observers "removal"; if not those then can be for dual observers with trial generator or 
           independent dual observers}
  \item{gamma}{slope for distance to create dependence except at y=0 or intercept and slope in delta0 to create dependence at all y}
  \item{epsilon}{if is.null(gamma), an additional form of heterogeneity is N(0,epsilon^2) added to logit of detection probability if epsilon single value
                      or N(0,(epsilon[1]+epsilon[2]*distance)^2) if not}
  \item{all}{if TRUE, then all objects within the point count are specified including those not seen}
  \item{posdep}{if TRUE, enforces positive dependence}
  \item{rotate}{if TRUE and either sigma or p0 values differ, then persons are rotated at half}
}
\details{
The sample size index works fine if there is no heterogeneity but will yield resulting samples
with fewer than n observations if there is heterogeneity.  Simply increase the index to achieve the n that you want.

Heterogeneity can be added to single observer using epsilon and for dual observer with either gamma or epsilon. In the latter
case, if gamma is non-null then it will be used instead of epsilon.  If \code{p1(y)} is the baseline detection probability at distance \code{y} defined
by \code{detfct} and \code{sigma}, then for each animal  \code{p(y)= plogis(log(p1(y)/(1-p1(y))) + rnorm(1,0,epsilon))} if heterogeneity is constant with
\code{y} (i.e., only one value of epsilong specified); otherwise, \code{p(y)= plogis(log(p1(y)/(1-p1(y)))+ rnorm(1,0,epsilon[1]+epsilon[2]*y))}.

For dual observers, heterogeneity can also be specified using values for gamma in \code{delta0(y)} as specified in Buckland et al.(2010).
If gamma has only one value then \code{delta0(y)=plogis(gamma[1]*y)} and if it has 2 values, then it is \code{delta0(y)=plogis(gamma[1]+gamma[2]*y)} which
yields dependence at \code{y=0}.

Note that because p0 is used in the logistic as \code{p1(y)= plogis(log(p0/(1-p0))+ sigma*y}, values of \code{p0} close to 1 will dampen the effect
of the slope and to get a decrease in detection probability with distance, it will need to be a large negative value.  If \code{p0} is set to 1 with 
a logistic detection function, it is reset to p0=0.9999999 to avoid divide by 0.

Differences between persons can be achieved by setting different p0_1 and p0_2 values or different sigma1 and sigma2 values.  If
\code{rotate=TRUE} then the persons are rotated between primary and secondary positions for half of the possible detections.

}
\value{dataframe with 2 records for each observation, one for primary and one for secondary
 \item{distance}{distance to observed animal}
 \item{observer}{1 for priamry and 2 for secondary}
 \item{detected}{1 if seen and 0 if missed}
 \item{person}{person identifier for observer}
 \item{pair}{person pairing in primary-secondary order}
}
\author{Jeff Laake}
\references{
Buckland, S. T., J. Laake and D. Borchers (2010). "Double observer line transect methods: levels of independence." Biometrics 66: 169-177.

Laake, J. L. and B. A. Collier (in prep). "Mark-Recapture Distance Sampling for Points"
}
\examples{
xx=gen.points(300,100,sigma1=20,detfct="hn",type="single",epsilon=c(0,.1))
xx=gen.points(350,100,-.05,detfct="logistic",type="removal",gamma=.02,p0_1=0.9)
}