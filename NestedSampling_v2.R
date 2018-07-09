library(sp)
library(gstat)
library(ggplot2)

# Function for random selection of 1 point at distance h from starting point
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

# Function for random selection of pair of points with separation distance h; starting point is halfway the pair of points
SelectPair<-function(start,h,area){
  dxy<-numeric(length=2)
  xypoints<-NULL
  inArea1 <- inArea2 <- NA
  while(is.na(inArea1) | is.na(inArea2)) {
    angle<-runif(n=1,min=0,max=2*pi)
    dxy[1]=h*sin(angle)/2
    dxy[2]=h*cos(angle)/2
    xypnt1<-start+dxy
    coordinates(xypnt1)<-~s1+s2
    inArea1<-as.numeric(over(x=xypnt1,y=area))[1]
    dxy[1]=-dxy[1]
    dxy[2]=-dxy[2]
    xypnt2<-start+dxy
    coordinates(xypnt2)<-~s1+s2
    inArea2<-as.numeric(over(x=xypnt2,y=area))[1]
  }
  xypoints<-rbind(as.data.frame(xypnt1),as.data.frame(xypnt2))
  xypoints
}

#Read field of interest

load(file="Data/HunterValley4Practicals.RData")
grd <- grdHunterValley
names(grd)[c(1,2)] <- c("s1","s2")

grid<-grd #data frame needed later for selection of main station
coordinates(grd) <- c("s1","s2")
gridded(grd) <- T


#Define lags in number of cells
cellsize<-25
lags<-c(120,40,20,10)

#Define for each stage the fraction of stations at which a pair of new stations is selected
fr <- c(1,1,1,1,1)
#For unbalanced designs choose 0.5 for last and/or one-but-last value of fr
#fr <- c(1,1,1,1,0.5)
#fr <- c(1,1,1,0.5,1)
#fr <- c(1,1,1,0.5,0.5)

# select main station
set.seed(614)
id <- sample.int(nrow(grid),1)
mainstation <- grid[id,c(1,2)] 

#select randomly one  point at distance lag[1] from main station
h=lags[1]*cellsize
pnt<-SelectPoint(start=mainstation,h=h,area=grd)
stations<-rbind(mainstation,pnt)
allunits <- NULL
for (j in 2:length(lags)) {
  if(fr[j]<1){
    nf <- nrow(stations)*fr[j]
    stations <- stations[1:nf,]
    #save stations not used in subsequent stages (these are units of nested sample)
    first <- nf+1
    last <- nrow(stations)
    units <- stations[first:last,]
    allunits <- rbind(allunits,units)
  }
  newstations<-NULL
  h=lags[j]*cellsize
  for (i in 1:nrow(stations)) {
    pnts<-SelectPair(start=stations[i,],h=h,area=grd)
    newstations <- rbind(newstations,pnts)
  }
  stations<-newstations
}

nestedsample <- rbind(allunits,stations)

# Assign factor levels to selected units; these factors are needed for AOV
# Number of factors is equal to number of stages - 1
# Note that code below must be adapted for unbalanced designs and when number of lags is changed

ntot <- nrow(nestedsample)
nestedsample$level1 <- rep(seq(1:2),each=ntot/2)
nestedsample$level2 <- rep(seq(1:2^2),each=ntot/2^2)
nestedsample$level3 <- rep(seq(1:2^3),each=ntot/2^3)
nestedsample$level4 <- rep(seq(1:2^4),each=ntot/2^4)

#Once data are collected, data can be analyzed as follows
library(nlme)
lmodel <- lme(z~1,data=nestedsample,random=~1|level1/level2/level3/level4)
out <- as.matrix(VarCorr(lmodel))
sigmas <- as.numeric(out[c(2,4,6,8,9),1])
ord <- sort(seq(1:length(lags)),decreasing=T)
sigmas <- sigmas[ord]
semivar <- cumsum(sigmas)
