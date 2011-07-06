################################################################################  
#  Functions to compute values for integrals/likelihood for removal distance sampling data 
#
#  Arguments:
#
#   distance       - vector of distance values
#   xmat           - design matrix for p(y)
#   beta           - coefficients for p(y)
#   p0             - p0 for intercept of logistic
#   width          - radius of point count or half-width of line transect
################################################################################  
s.1=function(xmat,beta,p0)
{
  p1=plogis(xmat%*%beta+log(p0/(1-p0)))
  return(p1)
}
logistic.pt.s1 <- function(distance,xmat,beta,width,p0) 
   return(s.1(xmat,beta,p0)*2*distance/width^2)

logistic.lt.s1 <- function(distance,xmat,beta,width,p0) 
   return(s.1(xmat,beta,p0)/width)
