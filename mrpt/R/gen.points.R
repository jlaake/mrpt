gen.points <-
function(n,w,sigma1,sigma2=NULL,detfct="hn",p0_1=0.9999999,p0_2=NULL,type="single",gamma=NULL,epsilon=NULL,all=FALSE, posdep=FALSE,rotate=TRUE)
{
# Generate simulated point transect data with optional 
#  capture history for 2 observers. No observer effect allowed p1=p2.
#
# Arguments:
# n       - index to sample size (also depends on sigma etc)
# w       - radius of point transect
# sigma1  - for observer 1 sigma for half-normal detection function or sigma = distance slope for logistic
# sigma2  - for observer 1 sigma for half-normal detection function or sigma = distance slope for logistic
# detfct  - either "hn" or "logistic"
# p0_1    - prob of detection at distance 0 for observer 1
# p0_2    - prob of detection at distance 0 for observer 2
# type    - single observer "single"; dependent dual observers "removal";
#           if not those then can be for dual observers with trial generator or 
#           independent dual observers.
# gamma   - slope for distance in delta0 to create dependence except at x=0
# epsilon - if gamma=NULL, an additional form of heterogeneity is N(0,epsilon^2)
#             added to logit of detection probability if epsilon single value
#             or N(0,(epsilon[1]+epsilon[2]*distance)^2) if not
# all     - if TRUE, then it reports distances to all objects and whether it was seen or not
#
compute.p=function(distance,p0,sigma,gamma,epsilon)
{
#  Compute detection probability for single observer
   n=length(distance)
   addhet=0
   if(!is.null(epsilon)& is.null(gamma))
     if(length(epsilon)==1)
        addhet=rnorm(n,0,epsilon)
     else
        addhet=rnorm(n,0,epsilon[1]+epsilon[2]*distance)       
   if(detfct=="hn")
   {
      p1=p0*exp(-distance^2/(2*sigma^2)) 
      p1=plogis(log(p1/(1-p1))+addhet)
   }
   else
   {
      if(p0==1)p0=0.9999999
      p1=plogis(log(p0/(1-p0))+distance*sigma+addhet)
   }  
   return(p1)
}
#  Generate random x,y coordinates over 2w unit square; restrict to circle
#  by excluding those with distance >w
   if(is.null(sigma2)) sigma2=sigma1
   if(is.null(p0_2)) p0_2=p0_1
   if(type=="single")
     k=1.1/sum(compute.p((0:100)*w/100,p0_1,sigma1,gamma,epsilon)*2*(0:100)*w/100/w^2) 
   else
   {
      p1=compute.p((0:100)*w/100,p0_1,sigma1,gamma,epsilon) 
      k=1.1/sum((2*p1-p1^2)*2*(0:100)*w/100/w^2)
   }
#   x=-w+2*w*runif(k*n)
#   y=-w+2*w*runif(k*n)
#  distance=sqrt(x^2+y^2)
#   distance=distance[distance<=w]
   distance=w*sqrt(runif(k*n))
   n=length(distance)
   n1=floor(n/2)
   n2=n-n1
# Compute probability for a primary and secondary observer
   if(!rotate)
   {
     person=data.frame(person1=rep("A",n),person2=rep("B",n))
     pair=data.frame(pair1=rep("AB",n),pair2=rep("AB",n))
     p1=compute.p(distance,p0_1,sigma1,gamma,epsilon)
     if(type!="single") p2=compute.p(distance,p0_2,sigma2,gamma,epsilon)
   }
   else
   {
     person=data.frame(person1=c(rep("A",n1),rep("B",n2)),person2=c(rep("B",n1),rep("A",n2)))
     pair=data.frame(pair1=c(rep("AB",n1),rep("BA",n2)),pair2=c(rep("AB",n1),rep("BA",n2)))
     p1=compute.p(distance,p0_1,sigma1,gamma,epsilon)
     p2=compute.p(distance,p0_2,sigma2,gamma,epsilon)
     px1=c(p1[1:n1],p2[(n1+1):n])
     px2=c(p2[1:n1],p1[(n1+1):n])
     p1=px1
     p2=px2
   }
#  For dual observers, compute delta dependence for each observation
   if(type!="single")
   {
      if(is.null(gamma))   
         delta=1
      else     
         if(length(gamma)==1)  
            delta=delta(matrix(distance,ncol=1),p1,p2,gamma,PI=TRUE,posdep=posdep,use.offset=TRUE)
         else
            delta=delta(model.matrix(~distance,data.frame(distance=distance)),p1,p2,gamma,PI=TRUE,posdep=posdep,use.offset=TRUE)         
#     Next compute probability of being detected by both observers (p11), by second observer only (p01),
#     and by first observer only (p10) depending on type.
      if(type=="removal")
      {
        p10=0
        p11=p1
        p01=p2*(1-p1*delta)
      }
      else
      {
         p11=p1*p2*delta
         p01=p1-p11
         p10=p01
      }
#     If dependent second observer (removal), set p10=0
#     Sum probabilities and compute p00 for those missed by both observers
      p00=1-p11-p01-p10
#     Generate multinomial random variables for capture history; exclude p01 for removal
      if(type=="removal")
         iclass=colSums(apply(cbind(p00,p11,p01),1,function(x)(1:3)*rmultinom(1,1,x)))      
      else
         iclass=colSums(apply(cbind(p00,p11,p01,p10),1,function(x)(1:4)*rmultinom(1,1,x)))
      if(!all)
        include=iclass!=1          
      else
        include=TRUE
      iclass=iclass[include]
#     Return dataframe with distances and 0/1 value for each observer
      xdf=data.frame(distance=rep(distance[include],2),observer=c(rep(1,length(distance[include])),rep(2,length(distance[include]))),
                         detected=c(as.numeric(iclass==2 | iclass==4),as.numeric(iclass==2 | iclass==3)),
                         person=unlist(person[include,]),pair=unlist(pair[include,]))
      row.names(xdf)=1:nrow(xdf)
      return(xdf)
   }
   else
   {
#     For single observer, generate bernoulli random variable for each distance
#     and return dataframe of observed distances.
      seen=rbinom(n,1,p1)
      if(!all)
        include=seen==1
      else
        include=TRUE
      return(data.frame(distance=distance[include],observer=rep(1,length(seen[include])),detected=seen[include]))
   }
}

