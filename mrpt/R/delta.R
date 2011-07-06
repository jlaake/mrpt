delta <-
function(x,p1,p2,gamma,PI,use.offset,posdep)
{
  if(all(gamma==0))
    return(rep(1,length(x))) 
  if(posdep)
  {
    gamma=exp(gamma)
    lower=rep(0,nrow(x))
  }
  else
    lower=apply(cbind(p1,p2),1,function(x)return(max(0,(x[1]+x[2]-1)/(x[1]*x[2]))))
  upper=apply(cbind(p1,p2),1,function(x)return(min(1/x[1],1/x[2])))
  if(PI | use.offset)
  {
     delta.values=rep(1,length(p1))
     offset=rep(0,length(p1))
     ix=(upper<1.0000001 & upper>.9999999) | ( lower<1.0000001 & lower>.9999999)
     offset[!ix]=log((1-lower[!ix])/(upper[!ix]-1))
     xx=(upper-lower)*plogis(x%*%gamma+offset)+lower
     delta.values[!ix]=xx[!ix]
     return(delta.values)
  }
  else
     return((upper-lower)*plogis(x%*%gamma)+lower)
}

