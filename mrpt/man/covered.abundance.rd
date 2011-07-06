\name{covered.abundance}
\alias{covered.abundance}
\alias{se.covered.abundance}
\title{Abundance in area covered by samples}
\description{Computes abundance in covered area in samples and its standard error}
\usage{
covered.abundance(par,model)
se.covered.abundance(model,delta=0.001)
}
\arguments{
  \item{par}{values for parameters}
  \item{model}{fitted model}
  \item{delta}{proportion of parameter value used for step size in computing numeric first derivative}
}
\details{
Computes abundance as sum of counts divided by the estimated detection probability for that count. Depending on the 
model there are differences in how those probabilities and counts are constructed.  For mark-recapture only models, the  
count and probability is specific to the distance (binned or unbinned) and covariate values; whereas, for MRDS models the count is specific
to sets of covariate values excluding distance and the probability is averaged over the distribution
of distances, pi(y) which is 1/W for lines and 2y/w^2 for points.
}
\value{Nhat from \code{covered.abundance} and se(Nhat) from \code{se.covered.abundance}}
\author{Jeff Laake}

