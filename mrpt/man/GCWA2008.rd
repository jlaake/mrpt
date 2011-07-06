\name{GCWA2008}
\alias{GCWA2008}
\docType{data}
\title{Golden-Cheeked Warbler Point Sampling Data}
\description{
A point sampling data set with two observers in a removal configuration in which the secondary observer (observer=2) is
aware of all detections by the primary observer (observer=1) but can detect those missed by the primary observer. The data 
collection protocols are described in more detail in  Wilkins et al. (2008) and Collier et al. (2010).
}
\usage{data(GCWA2008)}
\format{
A data frame with 1406 observations on the following 12 variables. The data represent 551 sightings of
golden-cheeked warblers with a record for the primary observer and secondary observer.  Thus 2*703=1102 records, but there
are only 2*551 records of birds because 2*152 records have NA distance and represent points with no birds seen.
  \describe{
    \item{\code{observer}}{1 for primary observer and 2 for secondary observer}
    \item{\code{detected}}{1 if seen by the observer and 0 if missed; detect=1 for all observer=2 sightings}
    \item{\code{distance}}{distance from the point center to the bird. distances were collected in bins of 0-50 and 50-100 m so the distance shown is the mid-point of 25 or 75}
    \item{\code{person}}{the person in that observer position; a factor with levels \code{A} to  \code{X}}
    \item{\code{observer.pair}}{the pair of persons in order of primary-secondary that were observing; a factor with multiple levels}
    \item{\code{experience.pair}}{the experience levels of the persons paired together; a factor with levels \code{00} \code{01} \code{10} \code{11}}
    \item{\code{experience}}{1 if the observer was experienced and 0 if not}
    \item{\code{ID}}{a factor variable for the point}
    \item{\code{PatchID}}{a factor variable for the patch containing the point}
    \item{\code{CC}}{canopy cover split into 2 levels \code{[0,66]} \code{(66,100]}}
    \item{\code{Date}}{date of survey split into approximate 2-week periods with \code{(0,14]} \code{(14,28]} \code{(28,42]} \code{(42,57]}}
    \item{\code{DaysFromPeak}}{number of days from the peak territory date}
  }
}
\source{ Wilkins et al. 2008; Collier et al. 2010}
\references{
Wilkins et al. 2008; Collier et al. 2010; Laake and Collier in prep
}
\examples{
data(GCWA2008)
GCWA2008=GCWA2008[!is.na(GCWA2008$distance),]
# Fit full independence model with distance
mod2008.FI=fit.removal(GCWA2008,width=100,cutp=c(0,50,100),indep=TRUE)
mod2008.FI
dev.new()
plot(mod2008.FI)
# Fit point independence model with distance and positive dependence
mod2008.PI=fit.removal(GCWA2008,width=100,cutp=c(0,50,100),PI=TRUE)
mod2008.PI
dev.new()
plot(mod2008.PI)
# Full independence with distance and canopy cover
mod2008.FIcc=fit.removal(GCWA2008,width=100,cutp=c(0,50,100),indep=TRUE,PI=FALSE,p.formula=~distance+CC)
mod2008.FIcc
dev.new()
plot(mod2008.FIcc)
# Full independence with distance and canopy cover
mod2008.PIcc=fit.removal(GCWA2008,width=100,cutp=c(0,50,100),PI=TRUE,p.formula=~distance+CC)
mod2008.PIcc
dev.new()
plot(mod2008.PIcc)
# Fit mark-recapture only model with distance 
mod2008.MR1=fit.removal(GCWA2008,width=100,cutp=c(0,50,100),mronly=TRUE)
mod2008.MR1
dev.new()
plot(mod2008.MR1)
mod2008.MR2=fit.removal(GCWA2008,width=100,cutp=c(0,50,100),mronly=TRUE,p.formula=~factor(distance))
mod2008.MR2
dev.new()
plot(mod2008.MR2)
\dontrun{
# Model MR1 can be fitted with RMark using the following:
#Create RMark  dataframe for first (11) and second (01) observers in MRDS removal model
gcwa.rem=data.frame(ch=c(rep("01",2),rep("11",2)),
                    freq=c(sum(1-GCWA2008$detected[GCWA2008$observer==1&GCWA2008$distance==25]),sum(1-GCWA2008$detected[GCWA2008$observer==1&GCWA2008$distance==75]), sum(GCWA2008$detected[GCWA2008$observer==1&GCWA2008$distance==25]),sum(GCWA2008$detected[GCWA2008$observer==1&GCWA2008$distance==75])), 
                    distance=rep(c(25,75),2),stringsAsFactors=FALSE)
freq=with(gcwa.rem,tapply(freq,distance,sum))
#Run Huggins closed capture model on removal MRDS data with 2 distance categories
Huggins.mod=mark(gcwa.rem,model="Huggins",model.parameters=list(p=list(formula=~distance),c=list(formula=~1,fixed=1)),output=FALSE,delete=TRUE)
#Estimate detection for each bin width  - here distance is an individual covariate so need to use covariate.predictions
plist=covariate.predictions(Huggins.mod,data=data.frame(index=c(1,1),distance=c(25,75)))
p=plist$estimates$estimate
pdot=2*p-p^2
deriv=-freq*2*(1-pdot)/pdot^2
#Huggins model output
Nhat.Huggins1=sum(freq/pdot)
Nhat.Huggins1.se=sqrt(sum(freq*(1-pdot)/pdot^2)+t(deriv)\%*\%plist$vcv[1:2,1:2]\%*\%deriv )
# Model MR2 can be fitted with RMark using the following:
# change distance to factor
gcwa.rem$distance=factor(gcwa.rem$distance)
Huggins.mod=mark(gcwa.rem,model="Huggins",groups="distance",model.parameters=list(p=list(formula=~distance),c=list(formula=~1,fixed=1)),output=FALSE,realvcv=TRUE)
#Estimate detection for each bin width
p=Huggins.mod$results$real$estimate[1:2]
pdot=2*p-p^2
deriv=-freq*2*(1-pdot)/pdot^2
#Huggins model output
Nhat.Huggins2=sum(Huggins.mod$results$derived$estimate)
Nhat.Huggins2.se=sqrt(sum(freq*(1-pdot)/pdot^2)+t(deriv)\%*\%Huggins.mod$results$real.vcv[1:2,1:2]\%*\%deriv )
# These estimates will not differ much but if more than 2 distance intervals are used it is
# preferable to fit with distance as an individual continuous covariate rather than as a factor variable
# which doesn't require the function to be smooth.
#
# There are slight differences in the std error between MARK and output from this package and not certain why
}
\dontrun{
# code used to create dataframe from raw data
GCWA2008=read.delim("C:/Users/Jeff Laake/Desktop/Consulting and Students/BretCollier/PointRemoval/MRDS.GCWA.2008.csv",header=T,sep=",")
GCWA2008$Observer=as.character(GCWA2008$Observer)
GCWA2008$Surveyor1=as.character(GCWA2008$Surveyor1)
GCWA2008$Surveyor2=as.character(GCWA2008$Surveyor2)
GCWA2008$Observer[GCWA2008$Observer=="AN"]="AN "
GCWA2008$Observer[GCWA2008$Obs2==0]="AN "
GCWA2008$Surveyor2[GCWA2008$Surveyor2=="SGB"]="SGD"
names(GCWA2008)[names(GCWA2008)=="Distance"]="distance"
GCWA2008$distance[GCWA2008$distance==1]=25
GCWA2008$distance[GCWA2008$distance==2]=75
GCWA2008$CC=cut(GCWA2008$CanopyCover,c(0,66,100),include=TRUE)
GCWA2008$Date=cut(GCWA2008$Days,c(0,14,28,42,57))
GCWA2008$DaysFromPeak=abs(GCWA2008$Days-35)
GCWA2008$Observer=factor(GCWA2008$Observer)
nametab=data.frame(Observer=unique(c(names(table(factor(GCWA2008$Observer))),names(table(factor(GCWA2008$Surveyor1))),names(table(factor(GCWA2008$Surveyor2))))),person=LETTERS[1:25]) 
GCWA2008=merge(GCWA2008,nametab,by="Observer",sort=FALSE)
names(GCWA2008)[names(GCWA2008)=="person"]="person.0"
GCWA2008=merge(GCWA2008,nametab,by.x="Surveyor1",by.y="Observer",sort=FALSE)
names(GCWA2008)[names(GCWA2008)=="person"]="person.1"
GCWA2008=merge(GCWA2008,nametab,by.x="Surveyor2",by.y="Observer",sort=FALSE)
names(GCWA2008)[names(GCWA2008)=="person"]="person.2"
GCWA2008$observer.pair=paste(GCWA2008$person.1,GCWA2008$person.2,sep="")
GCWA2008$experience.pair=paste(GCWA2008$Experience1,GCWA2008$Experience2,sep="")
x1=GCWA2008
x1$detected=x1$Obs1
x1$observer=1
x1$person=GCWA2008$person.1
x1$experience=GCWA2008$Experience1
x2=GCWA2008
x2$detected=x2$Obs2
x2$observer=2
x2$person=GCWA2008$person.2
x2$experience=GCWA2008$Experience2
x1$Obs1=NULL
x2$Obs1=NULL
x1$Obs2=NULL
x2$Obs2=NULL
GCWA2008= rbind(x1,x2)
GCWA2008=subset(GCWA2008,select=c("observer","detected","distance","person","observer.pair","experience.pair","experience","StationID","PatchID","CC","Date","DaysFromPeak"))
names(GCWA2008)[names(GCWA2008)=="StationID"]="ID"
GCWA2008$observer.pair=factor(GCWA2008$observer.pair)
GCWA2008$experience.pair=factor(GCWA2008$experience.pair)
save(GCWA2008,file="GCWA2008.rda")
}

}
\keyword{datasets}
