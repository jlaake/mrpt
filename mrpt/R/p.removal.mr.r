################################################################################
p.removal.mr <- function(par,x1,x2,models)    
{
################################################################################
#  Computes probabilities for capture histories in removal configuration
# 
#  Arguments:
#   par            - initial values for parameters
#   x1             - dataframe of double observer data for position 1
#   x2             - dataframe of double observer data for position 2
#   models         - list containing p.formula (formula for p(y)) and delta.formula (formula for delta(y))
################################################################################
xmat1=model.matrix(models$p.formula,x1)
xmat2=model.matrix(models$p.formula,x2)
# Extract parameter vectors: beta for detection and gamma for delta
beta=par[1:ncol(xmat1)]
# Extract/compute value of gamma dependening on posdep
gamma=0
p01=rem.p01(xmat1,xmat2,dmat=NULL,beta=beta,gamma=0,indep=TRUE,PI=FALSE,use.offset=TRUE,posdep=TRUE)
p11=rem.p11(xmat1,xmat2,dmat=NULL,beta=beta,gamma=0,indep=TRUE,PI=FALSE,use.offset=TRUE)
pdot=rem.pdot(xmat1,xmat2,dmat=NULL,beta=beta,gamma=0,indep=TRUE,PI=FALSE,use.offset=TRUE,posdep=TRUE)
return(list(p11=as.vector(p11),p01=as.vector(p01),pdot=as.vector(pdot)))
}

