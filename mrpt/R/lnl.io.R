lnl.io <-function(par,x1,x2,nfreq,models,width,cutp,debug,indep,
                       PI,use.offset,posdep,point)
{
if(debug)cat("\npar=",par,"\n")
# Compute probabilities 
p.list=p.io(par,x1,x2,models,width,cutp,indep=indep,PI=PI,use.offset=use.offset,
                    posdep=posdep,point=point)
p11=p.list$p11
p01=p.list$p01
p10=p.list$p10
p01[p01==0]=1e-6
p10[p10==0]=1e-6
# Compute log-likelihood value
lnl=sum(nfreq*(1-x1$Detected)*x2$Detected*log(p01))+ sum(nfreq*(1-x2$Detected)*x1$Detected*log(p10))+
   sum(nfreq*x1$Detected*x2$Detected*log(p11)) - sum(nfreq*log(p.list$pdot))
# If debug set, print out values
if(debug)
{
  print(summary(data.frame(p11=p11)))
  print(summary(data.frame(p01=p01)))
  print(summary(data.frame(p10=p10)))
  cat("\n-lnl=",-lnl)
}
return(-lnl)
}
