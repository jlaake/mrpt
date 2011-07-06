p.io <- function(par,x1,x2,models,width,cutp,indep,PI,use.offset,posdep,point)
{
# create design matrix for p1,p2 and delta
xmat1=model.matrix(models$p.formula,x1)
xmat2=model.matrix(models$p.formula,x2)
dmat=model.matrix(models$delta.formula,x1)
if(PI & all(dmat[,1]==1)) stop("\nError: No intercept allowed for delta(y) with PI model\n")
# extract parameter vectors: beta for detection and gamma for delta
beta=par[1:ncol(xmat1)]
if(!indep)
    gamma=par[(1+ncol(xmat1)):length(par)]
else
    gamma=0
# compute 2 integrals over cutpoints: g(r)r dr and delta(r)g(r)r dr
p01=rep(0,length=nrow(x1))
p10=rep(0,length=nrow(x1))
p11=rep(0,length=nrow(x1))
pdot=rep(0,length=nrow(x1))
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
       p01[i]=integratelogistic.io(logistic.pt.io01,x1[i,],x2[i,],models,beta,gamma,width,lower=lower,upper=upper,indep,PI,use.offset)
       p10[i]=integratelogistic.io(logistic.pt.io10,x1[i,],x2[i,],models,beta,gamma,width,lower=lower,upper=upper,indep,PI,use.offset)
       p11[i]=integratelogistic.io(logistic.pt.io11,x1[i,],x2[i,],models,beta,gamma,width,lower=lower,upper=upper,indep,PI,use.offset)
     }
     else
     {
       p01[i]=integratelogistic.io(logistic.lt.io01,x1[i,],x2[i,],models,beta,gamma,width,lower=lower,upper=upper,indep,PI,use.offset)
       p10[i]=integratelogistic.io(logistic.lt.io10,x1[i,],x2[i,],models,beta,gamma,width,lower=lower,upper=upper,indep,PI,use.offset)
       p11[i]=integratelogistic.io(logistic.lt.io11,x1[i,],x2[i,],models,beta,gamma,width,lower=lower,upper=upper,indep,PI,use.offset)
     }
   }
#  Otherwise it is a computation rather than an integral (unbinned distance)
   else
   {   
     if(point)
     {
       p01[i]=logistic.pt.io01(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width, indep=indep,PI=PI,use.offset=use.offset) 
       p10[i]=logistic.pt.io10(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width, indep=indep,PI=PI,use.offset=use.offset) 
       p11[i]=logistic.pt.io11(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width, indep=indep,PI=PI,use.offset=use.offset) 
     }
     else
     {
       p01[i]=logistic.lt.io01(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width, indep=indep,PI=PI,use.offset=use.offset) 
       p10[i]=logistic.lt.io10(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width, indep=indep,PI=PI,use.offset=use.offset) 
       p11[i]=logistic.lt.io11(lower,xmat1[i,,drop=FALSE],xmat2[i,,drop=FALSE],dmat[i,,drop=FALSE], beta, gamma, width, indep=indep,PI=PI,use.offset=use.offset) 
     }
   }
#  For all observations compute the integral over 0-W
   if(point)
      pdot[i]=integratelogistic.io(logistic.pt.iodot,x1[i,],x2[i,],models,beta,gamma,width,lower=0,upper=width,indep,PI,use.offset)
   else
      pdot[i]=integratelogistic.io(logistic.lt.iodot,x1[i,],x2[i,],models,beta,gamma,width,lower=0,upper=width,indep,PI,use.offset)
}
return(list(p11=p11,p01=p01,p10=p10,pdot=pdot))
}
