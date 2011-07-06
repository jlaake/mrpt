lnl.removal <-
function(par,x1,x2,nfreq,models,width,cutp,debug=FALSE,
              indep=FALSE,PI=TRUE,use.offset=TRUE,posdep=FALSE,point=TRUE,mronly=FALSE)
{
################################################################################
# Computes negative log likelihood for removal configuration of mrds data
#
# Arguments:
#   par            - vector of parameter values
#   xdf            - dataframe
#   xmat           - design matrix for p(y)
#   dmat           - design matrix for delta(y) 
#   nfreq          - frequency of observations 
#   width          - radius of point transect
#   cutp           - cut points for distances
#   debug          - if TRUE, output given during iteration
#   indep          - if TRUE, full independence assumed
#   PI             - if TRUE, point independence assumed; if !indep & !PI -- limiting independence 
#   use.offset     - should always be default of TRUE
#   posdep         - if TRUE forces positive dependence
#   point          - if TRUE for point counts; otherwise, line transect
#   mronly         - if TRUE, uses only mark-recapture likelihood
#
# Value: negative log-likelihood value
################################################################################
if(debug)cat("\npar=",par,"\n")
# compute probabilities 
if(mronly)
  p.list=p.removal.mr(par,x1,x2,models)
else
  p.list=p.removal(par,x1,x2,models,width,cutp,indep=indep,PI=PI,use.offset=use.offset,
                         posdep=posdep,point=point)
# if any are 0 set to a small value
p11=p.list$p11
p01=p.list$p01
p01[p01==0]=1e-6
if(debug)print(summary(data.frame(p11=p11)))
if(debug)print(summary(data.frame(p01=p01)))
if(debug)print(summary(data.frame(pdot=p.list$pdot)))
# compute negative log-likelihood value and return it
lnl=sum(nfreq*(1-x1$detected)*log(p01))+ sum(nfreq*x1$detected*log(p11))-
                                    sum(nfreq*log(p.list$pdot))
if(debug)cat("\n-lnl=",-lnl)
return(-lnl)
}

