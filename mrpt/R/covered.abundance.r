se.covered.abundance=function(model,delta=0.001)
{
  if(is.null(model$fit$hessian))
  {
    cat("\nHessian not available for computation of standard error of abundance. Returning NA\n")
    return(NA)
  }
  vcov=try(solve(model$fit$hessian))
  if(class(vcov)=="try-error")
  {
    cat("\nInversion of hessian for var-cov matrix failed. se set to NA\n")
    return(NA)
  } else
    return(sqrt(sum((1-model$p$pdot)*model$nfreq/model$p$pdot^2)+
        as.vector(DeltaMethod(model$fit$par, covered.abundance, vcov, delta=delta, model=model)$variance)))
}
covered.abundance=function(par,model)
{
  model$fit$par=par
  if(class(model)[2]=="removal")
     if(model$control$mronly)
        pp=p.removal.mr(model$fit$par,x1=model$data[model$data$observer==1,],x2=model$data[model$data$observer==2,],model$models)
     else
        pp=p.removal(model$fit$par,x1=model$data[model$data$observer==1,],x2=model$data[model$data$observer==2,],model$models,model$control$width,model$control$cutp,indep=model$control$indep,
                  PI=model$control$PI,posdep=model$control$posdep,point=model$control$point,use.offset=model$control$use.offset)
  else
    if(class(model)[2]=="single")
       pp=p.single(model$fit$par,model$data,model$models,model$control$width,model$control$cutp,
           point=model$control$point,type=model$control$type, p0=model$control$p0) 
  return(sum(model$nfreq/pp$pdot))
}
