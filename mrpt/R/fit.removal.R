##################################################################################
fit.removal <-function(x,par=NULL,p.formula=~distance,delta.formula=~-1+distance,width=max(c(x$distance,cutp)),cutp=NULL,debug=FALSE,
            indep=NULL,PI=NULL,method="Nelder-Mead",posdep=NULL,hessian=TRUE,point=TRUE,mronly=FALSE)
{
##################################################################################
# Fits point transect mrds removal model with binned or unbinned data
#
# Arguments:
#
#   x              - dataframe of double observer data
#   par            - initial values for parameters
#   p.formula      - formula for probability of detection
#   delta.formula  - formula for dependence function
#   width          - radius of point transect or half-width of line transect
#   cutp           - cut points for distances
#   debug          - if TRUE, output given during iteration
#   indep          - if TRUE, full independence assumed
#   PI             - if TRUE, point independence assumed; if !indep & !PI -- limiting independence 
#   method         - optimization method used in optim
#   posdep         - if TRUE forces positive dependence
#   hessian        - if TRUE, computes and returns hessian for fit
#   point          - if TRUE point count; otherwise line transect
#   mronly         - if TRUE, uses only mark-recapture likelihood
#       
# Value:
# List with elements 
#   model - optim results
#   AICc  - model selection value
##################################################################################
  if(mronly)
  {
   indep=TRUE
   PI=FALSE
  }
  if(is.null(indep))
  {
    if(is.null(PI))
    {
      PI=TRUE
      indep=FALSE
    } else
    {
      if(!PI)
        indep=TRUE
      else 
        indep=FALSE
    }
  } 
  if(indep)
  {
    PI=FALSE
    posdep=FALSE
  } else
    PI=TRUE
  if(PI &is.null(posdep)) posdep=TRUE  
  use.offset=TRUE
# Make sure all data are available
  if(is.null(x$distance))
    stop("\ndistance field missing in x\n")
  if(any(is.na(x$distance)))
    stop("\nsome distance values are missing (NA) in x\n")
  if(is.null(x$observer))
    stop("\nobserver field missing in x\n")
  if(any(is.na(x$observer)))
    stop("\nsome 0bserver values are missing (NA) in x\n")
  if(is.null(x$detected))
    stop("\ndetected field missing in x\n")
  if(any(is.na(x$detected)))
    stop("\nsome detected values are missing (NA) in x\n")
  fulldata=x
  x1=x[x$observer==1,]
  x2=x[x$observer==2,]
# Create cutp matrix from cutp breaks
  if(!is.null(cutp))
  {
    if(!is.matrix(cutp))
    {
      x1$distance=as.numeric(cut(x1$distance,breaks=cutp,include=TRUE))
      x2$distance=x1$distance
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
      x1$distance=as.numeric(cut(x1$distance,breaks=c(cutp[,1],cutp[nrow(cutp),2]),include=TRUE))
      x2$distance=x1$distance
    }
  }
  if(mronly & length(grep("factor(distance)",as.character(p.formula)))==0 &!is.null(cutp)) 
  {
     x1$distance=apply(cutp,1,mean)[x1$distance]    
     x2$distance=x1$distance
  }
# Construct design matrices for p(y) - xmat and delta(y) - dmat
  if(class(p.formula)=="formula")
  {
    xmat1=model.matrix(p.formula,x1) 
    xmat2=model.matrix(p.formula,x2) 
  }
  else
    stop("\np.formula must be a valid R formula.\n")
  if(class(delta.formula)=="formula")
  {
    dmat=model.matrix(delta.formula,x1)
    if(PI & all(dmat[,1]==1)) stop("\nError: No intercept allowed for delta(y) with PI model\n")
  }
  else
    stop("\ndelta.formula must be a valid R formula.\n")
  if(nrow(xmat1)!=nrow(xmat2) | nrow(xmat1)!=nrow(dmat) | nrow(xmat1)!=length(x1$distance))
    stop("\nSome data must be missing (NA) because number of rows in design matrices \nfor p.formula and delta.formula do not match or do not match number of rows in data.\n")
# If null par create initial values
  if(is.null(par)) 
  {
     par=rep(0,ncol(xmat1))
     if(is.null(cutp))  
     {
        cut1=max(width/4,3*width/sqrt(nrow(x)/2))
        cutn=width-cut1
     }
     else
     {
        cut1=cutp[2]
        cutn=cutp[length(cutp)-1]
     }
     seenby2=nrow(x[x$observer==2&x$distance<cut1,,drop=FALSE])
     seenby1=sum(x$detected[x$observer==1&x$distance<cut1])
     p0=(seenby1/seenby2)
     if(p0==1)p0=0.99
     par[1]=log(p0/(1-p0))
     seenby2=nrow(x[x$observer==2&x$distance>cutn,,drop=FALSE])
     seenby1=sum(x$detected[x$observer==1&x$distance>cutn])
     pc=seenby1/seenby2
     pdist=fit.single(x=x[x$observer==1&x$detected==1,], par = -1/width, p.formula = ~-1+distance, width = width, cutp = cutp, 
                     method = "Nelder-Mead", hessian = FALSE, point=point, type="logistic",p0=p0)$fit$par
     pw=plogis(log(p0/(1-p0))+width*pdist)
     par[colnames(xmat1)=="distance"]=pdist
#    if not all intercept or distance, fit an mronly model and fill in rest     
     notdistance=sapply(colnames(xmat1),function(x) as.numeric(length(grep("(Intercept)",x))>0)+
                                               as.numeric(length(grep("distance",x))>0)+
                                               as.numeric(length(grep("factor(distance)",x))>0))==0
     if(any(notdistance))
     {
        newpar=fit.removal(x=x,par=par,p.formula=p.formula,width=width,cutp=cutp,debug=FALSE,
              indep=TRUE,mronly=TRUE,hessian=FALSE,point=point)$fit$par
        par[notdistance]=newpar[notdistance]
     }     
     if(!indep)
     {
       if(posdep)
         par=c(par,rep(log(-log((1-pc)*pw/(pc*(1-pw)))/width),ncol(dmat)))
       else
         par=c(par,rep(-log((1-pc)*pw/(pc*(1-pw)))/width,ncol(dmat)))
     }
  } else
  {
#    Initial values can be specified as either a vector of values or a prior model
     if(class(par)[1]=="mrpt")
     {
        if(indep)
        {
          ipar=rep(0,ncol(xmat1))        
          names.match=colnames(xmat1)%in%colnames(par$dm$xmat1)
          vnames.match=colnames(par$dm$xmat1)%in%colnames(xmat1)
        } else
        {
          ipar=rep(0,ncol(xmat1)+ncol(dmat))
          names.match=c(colnames(xmat1),colnames(dmat))%in%c(colnames(par$dm$xmat1),colnames(par$dm$dmat))
          vnames.match=c(colnames(par$dm$xmat1),colnames(par$dm$dmat))%in%c(colnames(xmat1),colnames(dmat))
        }
        ipar[names.match]=par$fit$par[vnames.match]
        par=ipar       
     } else
     {
#       if vector of values specified just make sure the length is correct
        if(indep)
        {
           if(length(par)!=ncol(xmat1))
              stop("\nLength of par ",length(par)," does not match needed number of parameters ",ncol(xmat1)+ncol(dmat)," based on formula\n")
        } else
        {
            if(length(par)!=(ncol(xmat1)+ncol(dmat)))
              stop("\nLength of par ",length(par)," does not match needed number of parameters ",ncol(xmat1)+ncol(dmat)," based on formula\n")
        }
     }
   }
# Paste the 2 design matrices together and work out the number of unique values and collapse to 
# the set of unique values and counts (nfreq).  This will reduce the amount of numerical 
# integration in computation of the likelihood
  cmat=apply(cbind(xmat1,xmat2,dmat,x1$distance,x1$detected,x2$detected),1,paste,collapse="")
  nfreq=table(cmat)
  xmat1= xmat1[!duplicated(cmat),,drop=FALSE][order(unique(cmat)),,drop=FALSE]
  xmat2= xmat2[!duplicated(cmat),,drop=FALSE][order(unique(cmat)),,drop=FALSE]
  dmat=dmat[!duplicated(cmat),,drop=FALSE][order(unique(cmat)),,drop=FALSE]
  if("pair"%in%names(x))
    vars=unique(c("distance","observer","detected","pair",all.vars(p.formula),all.vars(delta.formula)))
  else
    vars=unique(c("distance","observer","detected",all.vars(p.formula),all.vars(delta.formula)))  
  x1=x1[!duplicated(cmat),,drop=FALSE][order(unique(cmat)),vars,drop=FALSE]
  x2=x2[!duplicated(cmat),,drop=FALSE][order(unique(cmat)),vars,drop=FALSE]
# Call optim to minimize the negative log likelihood for the removal configuration
  fit=optim(par,lnl.removal,x1=x1,x2=x2,nfreq=nfreq,models=list(p.formula=p.formula,delta.formula=delta.formula),
            width=width,cutp=cutp,debug=debug,indep=indep,PI=PI,use.offset=use.offset,method=method,
            hessian=FALSE,control=list(maxit=5000),posdep=posdep,point=point,mronly=mronly)
  if(hessian) 
     fit$hessian=hessian(lnl.removal,x=fit$par,method="Richardson",x1=x1,x2=x2,nfreq=nfreq,models=list(p.formula=p.formula,delta.formula=delta.formula),
            width=width,cutp=cutp,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep,point=point,mronly=mronly)
# Construct results list with fitted model, data, AICc value etc
  if(ncol(xmat1)==length(fit$par))
     gamma=0
  else
     gamma= fit$par[(ncol(xmat1)+1):length(fit$par)]
  mod=list(fit=fit,beta=fit$par[1:ncol(xmat1)],gamma=gamma,AICc=2*fit$value+2*length(fit$par)*(sum(nfreq)/(sum(nfreq)-length(fit$par)-1)),nfreq=nfreq,
          alldata=fulldata,data=rbind(x1,x2),dm=list(xmat1=xmat1,xmat2=xmat2,dmat=dmat),models=list(p.formula=p.formula,delta.formula=delta.formula),
          control=list(mronly=mronly,width=width,cutp=cutp,debug=debug,indep=indep,PI=PI,use.offset=use.offset,type="logistic",posdep=posdep,point=point))
  if(mronly)
     pp=p.removal.mr(mod$fit$par,x1,x2,mod$models)
  else
     pp=p.removal(mod$fit$par,x1,x2,mod$models,mod$control$width,mod$control$cutp,indep=mod$control$indep,
                  PI=mod$control$PI,posdep=mod$control$posdep,point=mod$control$point,use.offset=mod$control$use.offset)        
  mod$p=pp
  mod$Nhat=sum(mod$nfreq/pp$pdot)
  class(mod)=c("mrpt","removal")
  return(mod)
}
