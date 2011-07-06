fit.single <-function(x,par=NULL,p.formula=~distance,width=1,cutp=NULL,debug=FALSE,
                           method="Nelder-Mead",hessian=TRUE,point=TRUE,type="hn",p0=NULL)
{
# Fits point transect single observer data with binned or unbinned data
#
# Arguments:
#
#   x              - dataframe of double observer data
#   par            - initial values for parameters
#   p.formula      - formula for probability of detection
#   width          - radius of point transect or half-width of line transect
#   cutp           - cut points for distances
#   debug          - if TRUE, output given during iteration
#   method         - optimization method used in optim
#   hessian        - if TRUE, computes and returns hessian for fit
#   point          - if TRUE point count; otherwise line transect
#   type           - if "hn" fits half-normal; if "hr", hazard rate; otherwise a logistic
#       
# Value:
# List with elements 
#   model - optim results
#   AICc  - model selection value
##################################################################################
# Make sure all data are available
  if(is.null(x$distance))
    stop("\ndistance field missing in x\n")
  if(any(is.na(x$distance)))
    stop("\nsome distance values are missing (NA) in x\n")
  fulldata=x
# Create cutp matrix from cutp breaks
  if(!is.null(cutp))
  {
    if(!is.matrix(cutp))
    {
      x$distance=as.numeric(cut(x$distance,breaks=cutp,include=TRUE))
      if(length(cutp)>2)
         cutp=matrix(c(cutp[1],rep(cutp[2:(length(cutp)-1)],each=2),cutp[length(cutp)]),ncol=2,byrow=TRUE)
      else
         stop("\nMore than one interval must be specified by cutp.\n")
    } 
    else
    {
      if(dim(cutp)[1]<2)
         stop("\nMore than one interval must be specified by cutp.\n")
      else
        if(dim(cutp)[2]!=2)
          stop("\nWith cutp specified as matrix, it must have only 2 columns.\n")
      x$distance=as.numeric(cut(x$distance,breaks=c(cutp[,1],cutp[nrow(cutp),2])))
    }
  }      
# Construct design matrices for p(y) - xmat and delta(y) - dmat
  if(class(p.formula)=="formula")
    xmat=model.matrix(p.formula,x) 
  else
    stop("\np.formula must be a valid R formula.\n")
  if(nrow(xmat)!=length(x$distance))
    stop("\nSome data must be missing (NA) because number of rows in design matrices \nfor p.formula do not match number of rows in data.\n")
# Construct design matrices for p(y) - xmat and delta(y) - dmat
  xmat=model.matrix(p.formula,x)
# Check length of initial value vector to make sure it matches with formula
  if(length(par)!=ncol(xmat))
    stop("\nLength of par ",length(par)," does not match needed number of parameters ",ncol(xmat)," based on formula\n")
  if(type!="hn" &type!="hr" & all(xmat[,1]==1)) stop("\nError: No intercept allowed for p(y) with logistic model\n")
# Paste the 2 together and work out the number of unique values and collapse to 
# the set of unique values and counts (nfreq).  This will reduce the amount of numerical 
# integration in computation of the likelihood
  cmat=apply(cbind(xmat,x$distance),1,paste,collapse="")
  nfreq=table(cmat)
  xmat= xmat[!duplicated(cmat),,drop=FALSE][order(unique(cmat)),,drop=FALSE]
  vars=unique(c("distance",all.vars(p.formula)))
  x=x[!duplicated(cmat),vars,drop=FALSE][order(unique(cmat)),,drop=FALSE]
  if(is.null(par)) par=rep(0,ncol(xmat))
  if(type=="hr")par=c(2,par)
# Call optim to minimize the negative log likelihood for the removal configuration
  if(length(par)==1)
  {
    if(type=="logistic")
       fit=optimize(lnl.single,lower=par-.1,upper=par+.1,x=x,xmat=xmat,nfreq=nfreq,models=list(p.formula=p.formula),
            width=width,cutp=cutp,debug=debug,point=point,type=type,p0=p0)
    else    
       fit=optimize(lnl.single,lower=log(width/100),upper=log(width*4),x=x,xmat=xmat,nfreq=nfreq,models=list(p.formula=p.formula),
            width=width,cutp=cutp,debug=debug,point=point,type=type,p0=p0)
    fit$value=fit$objective
    fit$par=fit$minimum
    fit$minimum=NULL
    fit$objective=NULL
    fit$convergence=0
  }
  else
  fit=optim(par,lnl.single,x=x,xmat=xmat,nfreq=nfreq,models=list(p.formula=p.formula),
            width=width,cutp=cutp,debug=debug,method=method,
            hessian=FALSE,control=list(maxit=5000),point=point,type=type,p0=p0)
  if(hessian) 
     fit$hessian=hessian(lnl.single,x=fit$par,method="Richardson",xdf=x,xmat=xmat,
                    nfreq=nfreq,models=list(p.formula=p.formula),width=width,
                    cutp=cutp,point=point,type=type,p0=p0,debug=debug)
# Return resulting fitted model and AICc value
  mod=list(fit=fit,beta=fit$par,gamma=0,AICc=2*fit$value+2*length(fit$par)*(sum(nfreq)/(sum(nfreq)-length(fit$par)-1)),nfreq=nfreq,
          alldata=fulldata,data=x,dm=list(xmat=xmat),models=list(p.formula=p.formula),
          control=list(width=width,cutp=cutp,debug=debug,point=point,type=type,p0=p0))
  mod$p=p.single(mod$fit$par,mod$data,mod$models,mod$control$width,mod$control$cutp,point=mod$control$point,type=type,p0=p0)        
  mod$Nhat=sum(mod$nfreq/mod$p$pdot)
  class(mod)=c("mrpt","single")
  return(mod)
}

