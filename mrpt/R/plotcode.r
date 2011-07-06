fyhist.pt.plot=function(model,linesize=1)
{
  fyhist=as.vector(tapply(model$nfreq,model$data$distance[model$data$observer==1],sum)/(apply(model$control$cutp,1,diff)*sum(model$nfreq)))
  fyhist=c(as.vector(t(cbind(rep(0,length(fyhist)),fyhist,fyhist))),0)
  fyhist=data.frame(x=c(as.vector(t(cbind(model$control$cutp[,1],model$control$cutp))),model$control$width),y=fyhist)
  return(ggplot(fyhist,aes(x,y))+geom_step(size=linesize)+ xlab("\nDistance")+ylab("Probability Density\n"))
}

fy.pt.plot=function(model,distance=NULL,x=NULL,pdot=NULL)
{
   if(is.null(distance))distance=(0:100)*model$control$width/100
   if(!is.null(x))
   {
      if(nrow(x)>1)
        stop("\nDataframe x for design matrix can only contain one row\n")
      xlist1 <- as.list(x)
      xlist1$distance <-distance
      x1 <- expand.grid(xlist1)
      xmat=model.matrix(model$models$p.formula,x1)
      dmat=model.matrix(model$models$delta.formula,x1)
      fy=logistic.pt.remdot(distance,xmat1=xmat,xmat2=xmat,dmat=dmat,beta=model$beta,gamma=model$gamma,
              width=model$control$width,indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)/pdot     
   }
   else
   {
      vnames=names(model$data)
      vnames=vnames[!vnames%in%c("distance","detected")]
      if(length(vnames)!=0)
      {
        x=subset(model$data,select=vnames)
        row.names(x)=1:nrow(x)
        nr=row.names(unique(x))
        x$nfreq=as.vector(model$nfreq)
        xx=aggregate(x$nfreq,x[,-ncol(x),drop=FALSE],sum)
        nfreq=xx[,"x"]
        x=subset(xx,select=names(xx)[names(xx)!="x"])
        fy=rep(0,length(distance))
        rn=(1:nrow(x))[x$observer==1]
        for(i in 1:(nrow(x)/2))
        {
           xlist1 <- as.list(x[x$observer==1,,drop=FALSE][i,,drop=FALSE])
           xlist1$distance <-distance
           x1 <- expand.grid(xlist1)
           xmat1=model.matrix(model$models$p.formula,x1)
           xlist1 <- as.list(x[x$observer==2,,drop=FALSE][i,,drop=FALSE])
           xlist1$distance <-distance
           x1 <- expand.grid(xlist1)
           xmat2=model.matrix(model$models$p.formula,x1)
           dmat=model.matrix(model$models$delta.formula,x1)
           fy=fy+nfreq[rn[i]]*logistic.pt.remdot(distance,xmat1=xmat1,xmat2=xmat2,dmat=dmat,beta=model$beta,gamma=model$gamma,
              width=model$control$width,indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)/model$p$pdot[as.numeric(nr[i])]^2/model$Nhat 
        }
      } else
      {
          xmat = model.matrix(model$models$p.formula, data.frame(distance = distance))
          dmat = model.matrix(model$models$delta.formula, data.frame(distance = distance))
          fy=logistic.pt.remdot(distance,xmat1=xmat,xmat2=xmat,dmat=dmat,beta=model$beta,gamma=model$gamma,
              width=model$control$width,indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)/model$p$pdot[1]    
      }
   }   
   fy=data.frame(x=distance,y=fy)
   return(geom_line(data=fy))
}

gyhist.pt.plot=function(model,linesize=1,single=TRUE)
{
  if(single)
     gyhist=as.vector(tapply(model$nfreq[model$data$detected[model$data$observer==1]!=0],model$data$distance[model$data$observer==1][model$data$detected[model$data$observer==1]!=0],sum)/(apply(model$control$cutp^2,1,diff)/model$control$width^2*model$Nhat))
  else
     gyhist=as.vector(tapply(model$nfreq,model$data$distance[model$data$observer==1],sum)/(apply(model$control$cutp^2,1,diff)/model$control$width^2*model$Nhat))
  gyhist=c(as.vector(t(cbind(rep(0,length(gyhist)),gyhist,gyhist))),0)
  gyhist=data.frame(x=c(as.vector(t(cbind(model$control$cutp[,1],model$control$cutp))),model$control$width),y=gyhist)
  return(ggplot(gyhist,aes(x,y))+geom_step(size=linesize)+ xlab("\nDistance")+ylab("Pr(Detection)\n"))
}

gy.pt.plot=function(model,distance=NULL,x=NULL,pdot=NULL,single=TRUE)
{
   if(is.null(distance))distance=(0:100)*model$control$width/100
   if(!is.null(x))
   {
      if(nrow(x)>1)
        stop("\nDataframe x for design matrix can only contain one row\n")
      xlist1 <- as.list(x)
      xlist1$distance <-distance
      x1 <- expand.grid(xlist1)
      xmat=model.matrix(model$models$p.formula,x1)
      dmat=model.matrix(model$models$delta.formula,x1)
      if(single)
        gy=plogis(xmat%*%model$beta)
      else
        gy=rem.pdot(xmat1=xmat,xmat2=xmat,dmat=dmat,beta=model$beta,gamma=model$gamma,
                     indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)
   }
   else
   {
      vnames=names(model$data)
      vnames=vnames[!vnames%in%c("distance","detected")]
      if(length(vnames)!=0)
      {
        x=subset(model$data,select=vnames)
        row.names(x)=1:nrow(x)
        nr=row.names(unique(x))
        x$nfreq=as.vector(model$nfreq)
        xx=aggregate(x$nfreq,x[,-ncol(x),drop=FALSE],sum)
        nfreq=xx[,"x"]
        x=subset(xx,select=names(xx)[names(xx)!="x"])
        gy=rep(0,length(distance))
        rn=(1:nrow(x))[x$observer==1]
        for(i in 1:(nrow(x)/2))
        {
           xlist1 <- as.list(x[x$observer==1,,drop=FALSE][i,,drop=FALSE])
           xlist1$distance <-distance
           x1 <- expand.grid(xlist1)
           xmat1=model.matrix(model$models$p.formula,x1)
           xlist1 <- as.list(x[x$observer==2,,drop=FALSE][i,,drop=FALSE])
           xlist1$distance <-distance
           x1 <- expand.grid(xlist1)
           xmat2=model.matrix(model$models$p.formula,x1)
           dmat=model.matrix(model$models$delta.formula,x1)
           if(single)
             gy=gy+nfreq[rn[i]]*plogis(xmat1%*%model$beta)/model$p$pdot[as.numeric(nr[i])]/model$Nhat 
           else
             gy=gy+nfreq[rn[i]]*rem.pdot(xmat1=xmat1,xmat2=xmat2,dmat=dmat,beta=model$beta,gamma=model$gamma,
                     indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)/model$p$pdot[as.numeric(nr[i])]/model$Nhat 
        }
      } else
      {
          xmat = model.matrix(model$models$p.formula, data.frame(distance = distance))
          dmat = model.matrix(model$models$delta.formula, data.frame(distance = distance))
          if(single)
            gy=plogis(xmat%*%model$beta)
          else
            gy=rem.pdot(xmat1=xmat,xmat2=xmat,dmat=dmat,beta=model$beta,gamma=model$gamma,
                     indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)
      }
   }   
   gy=data.frame(x=distance,y=gy)
   return(geom_line(data=gy))
}

gcyhist.pt.plot=function(pc,model,linesize=1)
{
  gyhist=pc
  gyhist=c(as.vector(t(cbind(rep(0,length(gyhist)),gyhist,gyhist))),0)
  gyhist=data.frame(x=c(as.vector(t(cbind(model$control$cutp[,1],model$control$cutp))),model$control$width),y=gyhist)
  return(ggplot(gyhist,aes(x,y))+geom_step(size=linesize)+ xlab("\nDistance")+ylab("Conditional Pr(Detection)\n"))
}

gcy.pt.plot=function(model,distance=NULL,x=NULL,pdot=NULL,plot=TRUE)
{
   if(is.null(distance))distance=(0:100)*model$control$width/100
   if(!is.null(x))
   {
      if(nrow(x)>1)
        stop("\nDataframe x for design matrix can only contain one row\n")
      xlist1 <- as.list(x)
      xlist1$distance <-distance
      x1 <- expand.grid(xlist1)
      if(model$control$mronly&!is.null(model$control$cutp))
      {
         if(length(grep("factor(distance)",as.character(model$models$p.formula)))!=0)
         {
           cutp=model$control$cutp
           x1$distance=as.numeric(cut(x1$distance,breaks=c(cutp[,1],cutp[nrow(cutp),2]),include=TRUE))
         } else
         {
           cutp=model$control$cutp
           x1$distance=apply(cutp,1,mean)[as.numeric(cut(x1$distance,breaks=c(cutp[,1],cutp[nrow(cutp),2]),include=TRUE))]             
         }
      }
      xmat=model.matrix(model$models$p.formula,x1)
      dmat=model.matrix(model$models$delta.formula,x1)
      gy=rem.pc1(xmat1=xmat,xmat2=xmat,dmat=dmat,beta=model$beta,gamma=model$gamma,indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)
   }
   else
   {
      vnames=names(model$data)
      if(!model$control$mronly)
        vnames=vnames[!vnames%in%c("distance","detected")]
      else
        vnames=vnames[!vnames%in%c("detected")]
      if(length(vnames)!=0)
      {
        x=subset(model$data,select=vnames)
        row.names(x)=1:nrow(x)
        nr=row.names(unique(x))
        x$nfreq=as.vector(model$nfreq)
        xx=aggregate(x$nfreq,x[,-ncol(x),drop=FALSE],sum)
        nfreq=xx[,"x"]
        x=subset(xx,select=names(xx)[names(xx)!="x"])
        gy=rep(0,length(distance))
        rn=(1:nrow(x))[x$observer==1]
        for(i in 1:(nrow(x)/2))    
        {
           xlist1 <- as.list(x[x$observer==1,,drop=FALSE][i,,drop=FALSE])
           xlist1$distance <-distance
           x1 <- expand.grid(xlist1)
           if(model$control$mronly&!is.null(model$control$cutp))
           {
             if(length(grep("factor(distance)",as.character(model$models$p.formula)))!=0)
             {
               cutp=model$control$cutp
               x1$distance=as.numeric(cut(x1$distance,breaks=c(cutp[,1],cutp[nrow(cutp),2]),include=TRUE))
             } else
             {
               cutp=model$control$cutp
               x1$distance=apply(cutp,1,mean)[as.numeric(cut(x1$distance,breaks=c(cutp[,1],cutp[nrow(cutp),2]),include=TRUE))]             
             }
           }
           xmat1=model.matrix(model$models$p.formula,x1)
           xlist1 <- as.list(x[x$observer==2,,drop=FALSE][i,,drop=FALSE])
           xlist1$distance <-distance
           x1 <- expand.grid(xlist1)
           if(model$control$mronly&!is.null(model$control$cutp))
           {
             if(length(grep("factor(distance)",as.character(model$models$p.formula)))!=0)
             {
               cutp=model$control$cutp
               x1$distance=as.numeric(cut(x1$distance,breaks=c(cutp[,1],cutp[nrow(cutp),2]),include=TRUE))
             } else
             {
               cutp=model$control$cutp
               x1$distance=apply(cutp,1,mean)[as.numeric(cut(x1$distance,breaks=c(cutp[,1],cutp[nrow(cutp),2]),include=TRUE))]             
             }
           }
           xmat2=model.matrix(model$models$p.formula,x1)
           dmat=model.matrix(model$models$delta.formula,x1)
           gy=gy+nfreq[rn[i]]*rem.pc1(xmat1=xmat1,xmat2=xmat2,dmat=dmat,beta=model$beta,gamma=model$gamma,
             indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)/model$p$pdot[as.numeric(nr[i])]/model$Nhat 
        }
      } else
      {
          xmat = model.matrix(model$models$p.formula, data.frame(distance = distance))
          dmat = model.matrix(model$models$delta.formula, data.frame(distance = distance))
          gy=rem.pc1(xmat1=xmat,xmat2=xmat,dmat=dmat,beta=model$beta,gamma=model$gamma,indep=model$control$indep,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)
      }
   }                                              
   gy=data.frame(x=distance,y=gy)
   if(plot)
     return(geom_line(data=gy))
   else
     return(gy)
}


delta.pt.plot=function(model,distance=NULL,x=NULL,pdot=NULL)
{
   if(is.null(distance))distance=(0:100)*model$control$width/100
   if(!is.null(x))
   {
      if(nrow(x)>1)
        stop("\nDataframe x for design matrix can only contain one row\n")
      xlist1 <- as.list(x)
      xlist1$distance <-distance
      x1 <- expand.grid(xlist1)
      xmat=model.matrix(model$models$p.formula,x1)
      dmat=model.matrix(model$models$delta.formula,x1)
      p1=plogis(xmat%*%model$beta)
      deltay=delta(dmat,p1=p1,p2=p1,gamma=model$gamma,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)
   }
   else
   {
      vnames=names(model$data)
      vnames=vnames[!vnames%in%c("distance","detected")]
      if(length(vnames)!=0)
      {
        x=subset(model$data,select=vnames)
        row.names(x)=1:nrow(x)
        nr=row.names(unique(x))
        x$nfreq=as.vector(model$nfreq)
        xx=aggregate(x$nfreq,x[,-ncol(x),drop=FALSE],sum)
        nfreq=xx[,"x"]
        x=subset(xx,select=names(xx)[names(xx)!="x"])
        deltay=rep(0,length(distance))
        rn=(1:nrow(x))[x$observer==1]
        for(i in 1:(nrow(x)/2))
        {
           xlist1 <- as.list(x[x$observer==1,,drop=FALSE][i,,drop=FALSE])
           xlist1$distance <-distance
           x1 <- expand.grid(xlist1)
           xmat1=model.matrix(model$models$p.formula,x1)
           xlist1 <- as.list(x[x$observer==2,,drop=FALSE][i,,drop=FALSE])
           xlist1$distance <-distance
           x1 <- expand.grid(xlist1)
           xmat2=model.matrix(model$models$p.formula,x1)
           dmat=model.matrix(model$models$delta.formula,x1)
           p1=plogis(xmat1%*%model$beta)
           p2=plogis(xmat2%*%model$beta)
           deltay=deltay+nfreq[rn[i]]*delta(dmat,p1=p1,p2=p2,gamma=model$gamma,
            PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)/
            model$p$pdot[as.numeric(nr[i])]/model$Nhat 
        }
      } else
      {
          xmat = model.matrix(model$models$p.formula, data.frame(distance = distance))
          dmat = model.matrix(model$models$delta.formula, data.frame(distance = distance))
          p1=plogis(xmat%*%model$beta)
          deltay=delta(dmat,p1=p1,p2=p1,gamma=model$gamma,PI=model$control$PI,posdep=model$control$posdep,use.offset=model$control$use.offset)
      }
   }   
   deltay=data.frame(x=distance,y=deltay)
   if(model$control$indep)deltay$y=1
   return(ggplot(deltay,aes(x,y))+geom_line()+xlab("\nDistance")+ylab("Dependence\n"))
}

