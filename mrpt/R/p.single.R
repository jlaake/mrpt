################################################################################
p.single <- function(par,x,models,width,cutp,point,type,p0)    
{
################################################################################
#  Computes probabilities for capture histories in single configuration
# 
#  Arguments:
#   par            - initial values for parameters
#   x              - dataframe of double observer data
#   models         - list containing p.formula (formula for p(y)) 
#   width          - radius of point transect or half-width of line transect
#   cutp           - cut points for distances
#   point          - if TRUE for point counts; otherwise, line transect
################################################################################
xmat=model.matrix(models$p.formula,x)
beta=par
p1=rep(0,length=nrow(x))
pdot=rep(0,length=nrow(x))
# Loop over each value in x
for (i in 1:nrow(x))
{
#  Uses either a set of pre-defined cut points or dist.begin and dist.end fields in data
#  Work this out and set to lower and upper
   if(is.null(cutp))
   {
     if(is.null(x$dist.begin) || any(is.na(x$dist.begin)) || is.null(x$dist.end) || any(is.na(x$dist.end)))
     {
         lower=x[i,"distance"]
         upper=lower

     }
     else
     {
         lower=x[i,"dist.begin"]
         upper=x[i,"dist.end"]
     }
   }
   else
   {
      lower=cutp[x[i,"distance"],1]
      upper=cutp[x[i,"distance"],2]
   }
#  If lower!=upper then this needs to be an integral over the cutpoint range   
   if(lower!=upper)
   {
     if(point)
     {
       if(type=="hn" | type=="hr")
          p1[i]=integrate(detfct.pt,lower=lower,upper=upper,xmat=xmat[i,,drop=FALSE],beta=beta,width=width,type=type)$value
       else
          p1[i]=integratelogistic.single(logistic.pt.s1,x[i,],models,beta,width,lower=lower,upper=upper,p0=p0)
     }
     else
     {
       if(type=="hn" | type=="hr")    
          p1[i]=integrate(detfct.lt,lower=lower,upper=upper,xmat=xmat[i,,drop=FALSE],beta=beta,width=width,type=type)$value
       else
          p1[i]=integratelogistic.single(logistic.lt.s1,x[i,],models,beta,width,lower=lower,upper=upper,p0=p0)
     }       
   }
#  Otherwise it is a computation rather than an integral (unbinned distance)
   else
   {   
     if(point)
     {
       if(type=="hn" | type=="hr")    
          p1[i]=detfct.pt(lower,xmat[i,,drop=FALSE],beta, width,type=type) 
       else        
          p1[i]=logistic.pt.s1(lower,xmat[i,],beta,width,p0=p0)    
     }
     else
     {
       if(type=="hn" | type=="hr")    
          p1[i]=detfct.lt(lower,xmat[i,,drop=FALSE],beta, width,type=type) 
       else
          p1[i]=logistic.lt.s1(lower,xmat[i,],beta,width,p0=p0)           
     }
   }
#  For all observations compute the integral over 0-W
   if(point)
   {
      if(type=="hn" | type=="hr")    
         pdot[i]=integrate(detfct.pt,lower=0,upper=width,xmat=xmat[i,,drop=FALSE],beta=beta,width=width,type=type)$value
      else
        pdot[i]=integratelogistic.single(logistic.pt.s1,x[i,],models,beta,width,lower=0,upper=width,p0=p0)            
   }
   else
   {
      if(type=="hn" | type=="hr")    
         pdot[i]=integrate(detfct.lt,lower=0,upper=width,xmat=xmat[i,,drop=FALSE],beta=beta,width=width,type=type)$value
      else
        pdot[i]=integratelogistic.single(logistic.lt.s1,x[i,],models,beta,width,lower=0,upper=width,p0=p0)      
   }
}
return(list(p1=p1,pdot=pdot))
}

detfct.pt=function(distance,xmat,beta,width,type="hn")
{
  if(type=="hn")
    p1=exp(-(distance^2/(2*exp(xmat%*%beta)^2)))*2*distance/width^2
  else
  {
    key.shape=beta[1]
    key.scale=exp(xmat%*%beta[2:length(beta)])
    p1=(1 - exp( - (distance/key.scale)^( - key.shape)))*2*distance/width^2
  }    
  return(p1)
}

detfct.lt=function(distance,xmat,beta,width,type="hn")
{
  if(type=="hn")
    p1=exp(-(distance^2/(2*exp(xmat%*%beta)^2)))/width
  else
  {
    key.shape=beta[1]
    key.scale=exp(xmat%*%beta[2:length(beta)])
    p1=(1 - exp( - (distance/key.scale)^( - key.shape)))/width
  }    
  return(p1)
}
