################################################################################  
#  Functions to compute values for integrals/likelihood for removal distance sampling data 
#
#  Arguments:
#
#   distance       - vector of distance values
#   xmat1          - design matrix for p1(y)
#   xmat2          - design matrix for p2(y) - only observer differences if rotated
#   dmat           - design matrix for delta(y)
#   beta           - coefficients for p(y)
#   gamma          - coefficients for delta(y)
#   width          - radius of point count or half-width of line transect
#   indep          - if TRUE, full independence assumed
#   PI             - if TRUE, point independence assumed; 
#                         if !indep & !PI -- limiting independence 
#   use.offset     - should always be default of TRUE 
#   posdep         - if TRUE, restricts to positive dependence
################################################################################  
rem.p01=function(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
{
  p1=plogis(xmat1%*%beta)
  p2=plogis(xmat2%*%beta)
  if(indep)
    return(p2*(1-p1))
  else
  {
    xx=p2-p1*p2*delta(dmat,p1,p2,gamma,PI=PI,use.offset=use.offset,posdep=posdep)
    xx[xx<0]=0                           
    xx[xx>1]=1
    return(xx)
  }
}
rem.p11=function(xmat1,xmat2=NULL,dmat=NULL,beta,gamma,indep,PI,use.offset)
{
  return(plogis(xmat1%*%beta))
}
rem.pc1=function(xmat1,xmat2,dmat=NULL,beta,gamma,indep,PI,use.offset,posdep)
{
  p1=plogis(xmat1%*%beta)
  p2=plogis(xmat2%*%beta)
  return(p1*delta(dmat,p1,p2,gamma,PI=PI,use.offset=use.offset,posdep=posdep))
}
rem.pdot=function(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
{
  p1=plogis(xmat1%*%beta)
  p2=plogis(xmat2%*%beta)
  if(indep)
    return(p1+p2-p1*p2)
  else
  {
    xx=p1+p2-p1*p2*delta(dmat,p1,p2,gamma,PI=PI,use.offset=use.offset,posdep=posdep)
    xx[xx<0]=0                           
    xx[xx>1]=1
    return(xx)
  }
}
logistic.pt.rem01 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(rem.p01(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)*2*distance/width^2)

logistic.pt.rem11 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(rem.p11(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset)*2*distance/width^2)

logistic.pt.remdot <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(rem.pdot(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)*2*distance/width^2)

logistic.lt.rem01 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(rem.p01(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)/width)

logistic.lt.rem11 <- function(distance,xmat1, xmat2, dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(rem.p11(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset)/width)

logistic.lt.remdot <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(rem.pdot(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)/width)

