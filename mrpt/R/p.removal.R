################################################################################
p.removal <- function(par,x1,x2,models,width,cutp=NULL,indep=FALSE,
                   PI=TRUE,use.offset=TRUE,posdep=FALSE,point=TRUE)    
{
################################################################################
#  Computes probabilities for capture histories in removal configuration
# 
#  Arguments:
#   par            - initial values for parameters
#   x1             - dataframe of double observer data for position 1
#   x2             - dataframe of double observer data for position 2
#   models         - list containing p.formula (formula for p(y)) and delta.formula (formula for delta(y))
#   width          - radius of point transect or lalf-width of line transect
#   cutp           - cut points for distances
#   indep          - if TRUE, full independence assumed
#   PI             - if TRUE, point independence assumed; if !indep & !PI -- limiting independence 
#   use.offset     - should always be default of TRUE
#   posdep         - if TRUE forces positive dependence
#   point          - if TRUE for point counts; otherwise, line transect
################################################################################
xmat1=model.matrix(models$p.formula,x1)
xmat2=model.matrix(models$p.formula,x2)
dmat=model.matrix(models$delta.formula,x1)
# Extract parameter vectors: beta for detection and gamma for delta
beta=par[1:ncol(xmat1)]
# Extract/compute value of gamma dependening on posdep
if(!indep)
    gamma=par[(1+ncol(xmat1)):length(par)]
else
    gamma=0
p01=rep(0,length=nrow(x1))
p11=rep(0,length=nrow(x1))
pdot=rep(0,length=nrow(x1))
# Loop over each value in x
for (i in 1:nrow(x1))
{
#  Uses either a set of pre-defined cut points or dist.begin and dist.end fields in data
#  Work this out and set to lower and upper
   if(is.null(cutp))
   {
     if(is.null(x1$dist.begin) || any(is.na(x1$dist.begin)) || is.null(x1$dist.end) || any(is.na(x1$dist.end)))
     {
         lower=x1[i,"distance"]
         upper=lower

     }
     else
     {
         lower=x1[i,"dist.begin"]
         upper=x1[i,"dist.end"]
     }
   }
   else
   {
      lower=cutp[x1[i,"distance"],1]
      upper=cutp[x1[i,"distance"],2]
   }
#  If lower!=upper then this needs to be an integral over the cutpoint range   
   if(lower!=upper)
   {
     if(point)
     {
       p01[i]=integratelogistic.rem(logistic.pt.rem01,x1[i,],x2[i,],models,beta,gamma,width,
         lower=lower,upper=upper, indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)
       p11[i]=integratelogistic.rem(logistic.pt.rem11,x1[i,],x2[i,],models,beta,gamma,width,
         lower=lower,upper=upper,indep=indep,PI=PI,use.offset=use.offset) 
     }
     else
     {
       p01[i]=integratelogistic.rem(logistic.lt.rem01,x1[i,],x2[i,],models,beta,gamma,width,
          lower=lower,upper=upper, indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)
       p11[i]=integratelogistic.rem(logistic.lt.rem11,x1[i,],x2[i,],models,beta,gamma,width,
          lower=lower,upper=upper,indep=indep,PI=PI,use.offset=use.offset) 
     }
   }
#  Otherwise it is a computation rather than an integral (unbinned distance)
   else
   {   
     if(point)
     {
       p01[i]=logistic.pt.rem01(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width,
           indep=indep,PI=PI,use.offset=use.offset,posdep=posdep) 
       p11[i]=logistic.pt.rem11(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width,
           indep=indep,PI=PI,use.offset=use.offset) 
     }
     else
     {
       p01[i]=logistic.lt.rem01(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width,
          indep=indep,PI=PI,use.offset=use.offset,posdep=posdep) 
       p11[i]=logistic.lt.rem11(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width,
          indep=indep,PI=PI,use.offset=use.offset) 
     }
   }
#  For all observations compute the integral over 0-W
   if(point)
      pdot[i]=integratelogistic.rem(logistic.pt.remdot,x1[i,],x2[i,],models,beta,gamma,width,
        lower=0,upper=width, indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)
   else
      pdot[i]=integratelogistic.rem(logistic.lt.remdot,x1[i,],x2[i,],models,beta,gamma,width,lower=0,upper=width, 
         indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)
}
return(list(p11=p11,p01=p01,pdot=pdot))
}

