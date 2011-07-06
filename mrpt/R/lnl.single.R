lnl.single <-
function(par,xdf,xmat,nfreq,models,width,cutp,debug,point,type,p0)
{
################################################################################
# Computes negative log likelihood for removal configuration of mrds data
#
# Arguments:
#   par            - vector of parameter values
#   xdf            - dataframe
#   xmat           - design matrix for p(y)
#   nfreq          - frequency of observations 
#   width          - radius of point transect
#   cutp           - cut points for distances
#   debug          - if TRUE, output given during iteration
#   point          - if TRUE for point counts; otherwise, line transect
#
# Value: negative log-likelihood value
################################################################################
if(debug)cat("\npar=",par,"\n")
# compute probabilities 
p.list=p.single(par,xdf,models,width,cutp,point=point,type=type,p0=p0)
# if any are 0 set to a small value
# compute negative log-likelihood value and return it
lnl=sum(nfreq*log(p.list$p1/p.list$pdot))
if(debug)cat("\n-lnl=",-lnl)
return(-lnl)
}

