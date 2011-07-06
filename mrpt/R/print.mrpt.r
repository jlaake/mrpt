print.mrpt=function(x,...)
{
  if(x$control$point)
     cat("\nPoint Count Summary")
  else
     cat("\nLine Transect Summary")
  if(class(x)[2]=="single")
     cat("\nObserver Configuration: Single")
  else
     cat("\nObserver Configuration: Removal")
  cat("\nSample size           : ",sum(x$nfreq))
  if(!is.null(x$control$cutp))
    cat("\nDistance bins         : ",c(x$control$cutp[,1],x$control$cutp[nrow(x$control$cutp),2]))
  else
    cat("\nDistance unbinned")
  if(x$control$type=="hn")
     cat("\nDetection function    : Half-normal")
  else
    if(x$control$type=="hr")
       cat("\nDetection function    : Hazard rate")
    else
       cat("\nDetection function    : Logistic")
  cat("\nDetection formula     :",paste(x$models$p.formula,collapse=""))
  cat("\nAICc value            : ",x$AICc)
  if(x$fit$convergence!=0)cat("\nModel fit is suspect. Optimization did not converge\n")
  cat("\nSurveyed Abundance    : ",x$Nhat, "(SE=",se.covered.abundance(x),")\n")
  se=NULL
  if(!is.null(x$fit$hessian))
  {
    se=try(diag(solve(x$fit$hessian)))
    if(class(se)=="try-error")
    {
      se=rep(NA,length(x$fit$par))
      cat("\nWarning: Inversion of hessian for var-cov matrix failed. se set to NA\n")
    }
    else
    {
      se[se<0]=NA
      se=sqrt(se)
    }
  }
  if(class(x)[2]=="single")
    xmat=x$dm$xmat
  else
    xmat=x$dm$xmat1  
  py.range=1:ncol(xmat)
  if(x$control$type=="hr")py.range=1:(ncol(xmat)+1)
  cat("\nParameter estimates for p(y)\n")
  py.estimates=data.frame(estimate=x$fit$par[py.range])
  if(!is.null(se))
    py.estimates$se=se[py.range]
  if(x$control$type=="hr")
    row.names(py.estimates)=c("Power",colnames(xmat))
  else
    row.names(py.estimates)=colnames(xmat)
  print(py.estimates)
  if(class(x)[2]!="single")
  {
     if(x$control$indep)
       cat("\nFull Independence: delta0(y)=1\n")
     else
     {
       if(x$control$PI)
         cat("\nPoint Independence: Parameter estimates for delta0(y)\n")
       else
         cat("\nLimiting Independence: Parameter estimates for delta0(y)\n")
       delta.range=(max(py.range)+1):length(x$fit$par)
       cat("\n")
       deltay.estimates=data.frame(estimate=x$fit$par[delta.range])
       if(!is.null(se))
          deltay.estimates$se=se[delta.range]
       row.names(deltay.estimates)=colnames(x$dm$dmat)
       print(deltay.estimates)
     }
  }
  invisible()
}
plot.mrpt=function(x,distance=NULL,titles="",single=TRUE,...)
{

  dummy=function(x)return(x)
######### functions for multiple ggplots on same page
  vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
  arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
   dots <- list(...)
   n <- length(dots)
   if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
   if(is.null(nrow)) { nrow = ceiling(n/ncol)}
   if(is.null(ncol)) { ncol = ceiling(n/nrow)}
        ## NOTE see n2mfrow in grDevices for possible alternative
   grid.newpage()
   pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
   ii.p <- 1
   for(ii.row in seq(1, nrow)){
   ii.table.row <- ii.row
   if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
     ii.table <- ii.p
     if(ii.p > n) break
     print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
     ii.p <- ii.p + 1
    }
   }
  }
#########
  model=x
  if(class(model)[2]=="removal")
  {
    if(is.null(model$control$cutp))
    {
      k=max(floor(.67*sqrt(nrow(model$alldata)/2)),2)
      if(!model$control$point)
          cutp=(0:k)*(model$control$width/k)
      else
          cutp=seq(0,k,2)*(model$control$width/k)  
      cutp=matrix(c(cutp[1],rep(cutp[2:(length(cutp)-1)],each=2),cutp[length(cutp)]),ncol=2,byrow=TRUE)
    }  
    else
       cutp=model$control$cutp
    if(length(grep("factor(distance)",as.character(model$models$p.formula)))==0) 
        p.formula=as.formula(gsub("distance","factor(distance)",model$models$p.formula))
    mrmodel=fit.removal(model$alldata,width=model$control$width,cutp=cutp,p.formula=p.formula,indep=TRUE,mronly=TRUE)
    pc=gcy.pt.plot(mrmodel,distance=apply(cutp,1,mean),plot=FALSE)$y
    if(model$control$mronly)
    {
       return(gcyhist.pt.plot(pc,mrmodel)+gcy.pt.plot(model,distance=distance))
    }
    else
    if(all(titles==""))
    {
    if(single)
      arrange(fyhist.pt.plot(model)+fy.pt.plot(model,distance=distance),
        gyhist.pt.plot(model,single=FALSE)+gy.pt.plot(model,distance=distance,single=FALSE), 
        gyhist.pt.plot(model,single=TRUE)+gy.pt.plot(model,distance=distance,single=TRUE), 
        gcyhist.pt.plot(pc,model)+gcy.pt.plot(model,distance=distance),
        delta.pt.plot(model),nrow=3,ncol=2)                   
    else
      arrange(fyhist.pt.plot(model)+fy.pt.plot(model,distance=distance),
        gyhist.pt.plot(model,single=FALSE)+gy.pt.plot(model,distance=distance,single=FALSE), 
        gcyhist.pt.plot(pc,model)+gcy.pt.plot(model,distance=distance),
        delta.pt.plot(model),nrow=2,ncol=2)                   
    }
    else
    {
    if(single)
      arrange(fyhist.pt.plot(model)+fy.pt.plot(model,distance=distance)+opts(title=titles[1]),
        gyhist.pt.plot(model,single=FALSE)+gy.pt.plot(model,distance=distance,single=FALSE)+opts(title=titles[2]), 
        gyhist.pt.plot(model,single=TRUE)+gy.pt.plot(model,distance=distance,single=TRUE)+opts(title=titles[3]), 
        gcyhist.pt.plot(pc,model)+gcy.pt.plot(model,distance=distance)+opts(title=titles[4]),
        delta.pt.plot(model)+opts(title=titles[5]),nrow=3,ncol=2)                   
    else     
      arrange(fyhist.pt.plot(model)+fy.pt.plot(model,distance=distance)+opts(title=titles[1]),
        gyhist.pt.plot(model,single=FALSE)+gy.pt.plot(model,distance=distance,single=FALSE)+opts(title=titles[2]), 
        gcyhist.pt.plot(pc,model)+gcy.pt.plot(model,distance=distance)+opts(title=titles[3]),
        delta.pt.plot(model)+opts(title=titles[4]),nrow=2,ncol=2)  
    }                 
  }
  invisible()
}
