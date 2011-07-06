################################################################################  
#  Functions to compute values for integrals/likelihood for io distance sampling data 
#
#  Arguments:
#
#   distance       - vector of distance values
#   xmat1          - design matrix for p(y) for position 1
#   xmat2          - design matrix for p(y) for position 2
#   dmat           - design matrix for delta(y)
#   beta           - coefficients for p(y)
#   gamma          - coefficients for delta(y)
#   width          - radius of point count or half-width of line transect
#   indep          - if TRUE, full independence assumed
#   PI             - if TRUE, point independence assumed; 
#                         if !indep & !PI -- limiting independence 
#   use.offset     - should always be default of TRUE 
################################################################################  
io.p01=function(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
{
  p1=plogis(xmat1%*%beta)
  p2=plogis(xmat2%*%beta)
  if(indep)
    return(p2*(1-p1))
  else
  {
    xx=p2-p1*p2*delta(dmat,p1,p2,gamma,PI,use.offset,posdep)
    xx[xx<0]=0                           
    xx[xx>1]=1
    return(xx)
  }
}
io.p10=function(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
{
  p1=plogis(xmat1%*%beta)
  p2=plogis(xmat2%*%beta)
  if(indep)
    return(p1*(1-p2))
  else
  {
    xx=p1-p1*p2*delta(dmat,p1,p2,gamma,PI,use.offset,posdep)
    xx[xx<0]=0                           
    xx[xx>1]=1
    return(xx)
  }
}

io.p11=function(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
{
  p1=plogis(xmat1%*%beta)
  p2=plogis(xmat2%*%beta)
  if(indep)
    return(p1*p2)
  else
  {
    xx=p1*p2*delta(dmat,p1,p2,gamma,PI,use.offset,posdep)
    xx[xx<0]=0                           
    xx[xx>1]=1
    return(xx)
  }
}

io.pdot=function(xmat1,xmat2,dmat,beta,gamma,indep,PI,use.offset,posdep)
{
  p1=plogis(xmat1%*%beta)
  p2=plogis(xmat2%*%beta)
  if(indep)
    return(p1+p2-p1*p2)
  else
  {
    xx=p1+p2-p1*p2*delta(dmat,p1,p2,gamma,PI,use.offset,posdep)
    xx[xx<0]=0                           
    xx[xx>1]=1
    return(xx)
  }
}
logistic.pt.io01 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(io.p01(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)*2*distance/width^2)

logistic.pt.io10 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(io.p10(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)*2*distance/width^2)

logistic.pt.io11 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(io.p11(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)*2*distance/width^2)

logistic.pt.iodot <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(io.pdot(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)*2*distance/width^2)

logistic.lt.io01 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(io.p01(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)/width)

logistic.lt.io10 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(io.p10(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)/width)

logistic.lt.io11 <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(io.p11(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)/width)

logistic.lt.iodot <- function(distance,xmat1,xmat2,dmat,beta,gamma,width,indep,PI,use.offset,posdep) 
   return(io.pdot(xmat1,xmat2,dmat,beta,gamma,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)/width)

