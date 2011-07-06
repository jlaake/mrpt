################################################################################  
#  Computes integral (numerically) over y from lower to upper of a
#     logistic type detection function defined by f.
#
#  Arguments:
#
#   fx            - function to be integrated
#   x             - covariate data for observer1 with multiple distance values for integration
#   models        - list containing p.formula (formula for p(y)) and delta.formula (formula for delta(y))
#   beta           - coefficients for p(y)
#   gamma          - coefficients for delta(y)
#   width          - radius of point count
#   lower,upper    - bounds for integration
#   indep          - if TRUE, full independence assumed
#   PI             - if TRUE, point independence assumed; if !indep & !PI -- limiting independence 
#   use.offset     - should always be default of TRUE
#
################################################################################  
integratelogistic.single <-
function (fx, x, models, beta, width, lower=0,upper=width,p0)
{
  integrate(logisticbyx.single,lower=lower,upper=upper, subdivisions=10, rel.tol=0.001,
          stop.on.error=FALSE,fx=fx,x=x, models=models, beta=beta, 
          width=width,p0=p0)$value
}
integratelogistic.rem <-
function (fx, x1, x2, models, beta, gamma, width, lower=0,upper=width,indep,PI,use.offset,posdep)
{
  integrate(logisticbyx.rem,lower=lower,upper=upper, subdivisions=10, rel.tol=0.001,
          stop.on.error=FALSE,fx=fx,x1=x1,x2=x2, models=models, beta=beta, 
          gamma=gamma,indep=indep,PI=PI,use.offset=use.offset,width=width,posdep=posdep)$value
}
integratelogistic.io <-
function (fx, x1, x2, models, beta, gamma, width, lower=0,upper=width,indep,PI,use.offset,posdep)
{
  integrate(logisticbyx.io,lower=lower,upper=upper, subdivisions=10, rel.tol=0.001,
          stop.on.error=FALSE,fx=fx,x1=x1,x2=x2, models=models, beta=beta, 
          gamma=gamma,indep=indep,PI=PI,use.offset=use.offset,width=width,posdep=posdep)$value
}
################################################################################  
#  logisticbyx - treats logistic as a function of covariate z; for a given z it computes
#                function with those covariate values at a range of distances
#
#  Arguments:                                 
#
#   distance       - vector of distance values
#   x              - covariate data for observer1 with multiple distance values for integration
#   models         - list containing p.formula (formula for p(y)) and delta.formula (formula for delta(y))
#   beta           - coefficients for p(y)
#   gamma          - coefficients for delta(y)
#   width          - radius of point count
#   indep          - if TRUE, full independence assumed
#   PI             - if TRUE, point independence assumed; if !indep & !PI -- limiting independence 
#   use.offset     - should always be default of TRUE
#   fx             - function to be integrated
#
#  value: vector of values for integral calculation
#
################################################################################  
logisticbyx.single <-                                                                               
function (distance,x,models,beta,width,fx,p0)
{
   xlist1 <- as.list(x)
   xlist1$distance <-distance
   x1 <- expand.grid(xlist1)
   xmat=model.matrix(models$p.formula,x1)
   xx=fx(distance=distance,xmat=xmat,beta=beta,width=width,p0=p0)
   return(xx)
}
logisticbyx.rem <-
function (distance,x1,x2,models,beta,gamma,width,indep,PI,use.offset,posdep,fx)
{
   xlist1 <- as.list(x1)
   xlist1$distance <-distance
   x1 <- expand.grid(xlist1)
   xlist2 <- as.list(x2)
   xlist2$distance <-distance
   x2 <- expand.grid(xlist2)
   xmat1=model.matrix(models$p.formula,x1)
   xmat2=model.matrix(models$p.formula,x2)
   dmat=model.matrix(models$delta.formula,x1)
   xx=fx(distance=distance,xmat1=xmat1,xmat2=xmat2,dmat=dmat,beta=beta,gamma=gamma,width=width,
                   indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)
   return(xx)
}
logisticbyx.io <-
function (distance,x1,x2,models,beta,gamma,width,indep,PI,use.offset,posdep,fx)
{
   xlist1 <- as.list(x1)
   xlist1$distance <-distance
   x1 <- expand.grid(xlist1)
   xlist2 <- as.list(x2)
   xlist2$distance <-distance
   x2 <- expand.grid(xlist2)
   xmat1=model.matrix(models$p.formula,x1)
   xmat2=model.matrix(models$p.formula,x2)
   dmat=model.matrix(models$delta.formula,x1)
   xx=fx(distance=distance,xmat1=xmat1,xmat2=xmat2,dmat=dmat,beta=beta,gamma=gamma,
               width=width,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)
   return(xx)
}


