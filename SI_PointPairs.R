# load relevant packages
library(sp)
library(ggplot2)

#Define fucntion for simple random sampling of pairs of points with separation disatnce h

SIpairs<-function(h,n,area){
    i<-1
    xy<-coordinates(area)
    y<-matrix(nrow=n,ncol=2)
    dx<-numeric(length=2)
    while (i <= n) {
        id1<-sample(x=1:length(area),1)
        y[i,1]<-area$z[id1]
        xy1st<-xy[id1,]
        angle<-runif(n=1,min=0,max=2*pi)
        dx[1]=h*sin(angle)
        dx[2]=h*cos(angle)
        xy2nd<-xy1st+dx
        xy2nd<-as.data.frame(t(xy2nd))
        coordinates(xy2nd)<-~s1+s2
        res<-as.numeric(over(x=xy2nd,y=area))[1]
        if (!is.na(res)) {y[i,2]<-res 
                          i<-i+1}
        rm(xy2nd)
    }
    ydf<-as.data.frame(y)
    names(ydf)<-c("z1","z2")
    ydf$h<-h
    return(ydf)
}

#Read simulated field

grd <- read.csv(file="Data/SimulatedField_1Exp(25).csv")
names(grd)[3] <- "z"
coordinates(grd) <- ~ s1+s2
gridded(grd) <- T

cellsize<-1

#Randomly select two points at fixed distance, and repeat this many times
set.seed(314)
samplesize<-100

#Give separation distances in number of cells
h<-c(1,3,9,27,54)
h<-h*cellsize

allpairs<-NULL
for (i in 1:length(h)){
  pairs<-SIpairs(h=h[i],n=samplesize,area=grd)
  allpairs <- rbind(allpairs,pairs)
}

gammah<-vargammah<-numeric(length=length(h))

for (i in 1:length(h)){
  ids<-which(allpairs$h==h[i])
  mysample<-allpairs[ids,]
  gammah[i]<-mean((mysample$z1-mysample$z2)^2,na.rm=TRUE)/2 
  vargammah[i]<-var((mysample$z1-mysample$z2)^2,na.rm=TRUE)/(samplesize*4)
}

#Plot sample variogram

samplevariogram<-data.frame(h,gammah,vargammah)
ggplot(data=samplevariogram) +
    geom_point(mapping = aes(x = h,y=gammah),size = 3) +
    scale_x_continuous(name = "Separation distance",limits=c(0,60)) +
    scale_y_continuous(name="Semivariance",limits=c(0,1))

#Fit model

spherical <- function(h, range, psill) {
  h <- h/range
  psill*ifelse(h < 1, (1.5 * h - 0.5 * h^3),1)
}

sphericalnugget <- function(h, range, psill, nugget) {
    h <- h/range
    nugget + psill*ifelse(h < 1, (1.5 * h - 0.5 * h^3),1)
}

exponential <- function(h, range, psill, nugget) {
  h <- h/range
  psill*(1-(exp(-h)))
}

exponentialnugget <- function(h, range, psill, nugget) {
    h <- h/range
    nugget + psill*(1-(exp(-h)))
}


fit.var <- nls(gammah~exponential(h,range,psill),
               data = samplevariogram,
               start=list(psill=0.9, range=30),
               weights=1/vargammah,
               trace=T)

allpars<-NULL
nboot<-10
for (j in 1:nboot) {
  #select bootstrap sample for each lag and compute semivariance
  gammah<-vargammah<-numeric(length=length(h))
  for (i in 1:length(h)){
    ids<-which(allpairs$h==h[i])
    pairs<-allpairs[ids,]
    mysampleids<-sample.int(samplesize,size=samplesize,replace=TRUE)
    mysample<-pairs[mysampleids,]
    gammah[i]<-mean((mysample$z1-mysample$z2)^2,na.rm=TRUE)/2 
    vargammah[i]<-var((mysample$z1-mysample$z2)^2,na.rm=TRUE)/(samplesize*4)
  }
  
  
  #fit model
  samplevariogram<-data.frame(h,gammah,vargammah)
  fittedvariogram <- nls(gammah~exponential(h,range,psill),
                         data = samplevariogram,
                         start=list(psill=0.9, range=30),
                         weights=1/vargammah) 
  pars<-coef(fittedvariogram)
  allpars<-rbind(allpars,pars)
}

#compute variance-covariance matrix of two variogram parameters
var(allpars)
