# Version of nested sampling in which in each stage 1 point is selected 
# at some chosen distance h from all points selected in previous stages

# Only for balanced nested samples!

library(sp)
library(gstat)
library(ggplot2)

#Define fucntion for random selection of point at some chosen distance from starting point
SelectPoint<-function(start,h,area){
  dxy<-numeric(length=2)
  inArea<-NA
  while(is.na(inArea)) {
    angle<-runif(n=1,min=0,max=2*pi)
    dxy[1]=h*sin(angle)
    dxy[2]=h*cos(angle)
    xypnt<-start+dxy
    coordinates(xypnt)<-~s1+s2
    inArea<-as.numeric(over(x=xypnt,y=area))[1]
  }
  xypoint<-as.data.frame(xypnt)
  xypoint
}
#Read field of interest

load(file="Data/HunterValley4Practicals.RData")
grd <- grdHunterValley
names(grd)[c(1,2)] <- c("s1","s2")


grid<-grd #data frame needed later for selection of main station
coordinates(grd) <- c("s1","s2")
gridded(grd) <- T

#Define separation distances
lags<-c(625,125,25,5,1)

# select main station
set.seed(314)
id <- sample.int(nrow(grid),1)
mainstation <- grid[id,c(1,2)] 

#select randomly one  point at distance lag[1] from main station
newpnt<-SelectPoint(start=mainstation,h=lags[1],area=grd)
allstarts<-rbind(mainstation,newpnt)

for (j in 2:length(lags)) {
  newpnts<-NULL
  for (i in 1:nrow(allstarts)) {
    pnts<-SelectPoint(start=allstarts[i,],h=lags[j],area=grd)
    newpnts <- rbind(newpnts,pnts)
  }
  allstarts <- rbind(allstarts,newpnts)
}

nestedsample <- allstarts

#Note that code below must be adapted when number of lags is changed
nestedsample$l1<-rep(seq(1:2))
nestedsample$l2<-rep(rep(seq(1:2),each=2),times=8)
nestedsample$l3<-rep(rep(seq(1:2),each=4),times=4)
nestedsample$l4<-rep(rep(seq(1:2),each=8),times=2)

#Once data are collected, data can be analyzed as follows
library(nlme)
lmodel <- lme(z~1,data=nestedsample,random=~1|l1/l2/l3/l4)
out <- as.matrix(VarCorr(lmodel))
sigmas <- as.numeric(out[c(2,4,6,8,9),1])
ord <- sort(seq(1:length(lags)),decreasing=T)
sigmas <- sigmas[ord]
semivar <- cumsum(sigmas)

