library(sampling)
library(raster)

#read data
rmap <- raster("Data/geo.tif")
grd <- as(rmap,"SpatialPixelsDataFrame")
gridded(grd)<-F
grd <- as(grd,"data.frame")

#delete rows with geo = 99 (water)
ids <- which(grd$geo==99)
grd <- grd[-ids,]

#select stratified simple random sample

#first sort units on geo
unique(grd$geo)
grd <- grd[order(grd$geo),]

#compute stratum sample sizes for proportional allocation
Nh <- tapply(grd$x,INDEX=grd$geo,FUN=length)
wh <- Nh/sum(Nh)

#set total sample size
n <- 62
nh <- round(wh*n)
sum(nh)

#minimum stratum sample size is 2, so increase nh for stratum 2 by 1, and reduce sample szie of stratum 6 by 1
nh[2] <- nh[2]+1
nh[6] <- nh[6]-1

set.seed(314)
units<-strata(grd,stratanames="geo",size=nh,method="srswr")
mysample<-getdata(grd,units)

# select random location within selected pixel
resolution <-res(rmap)
mysample$x <- jitter(mysample$x,resolution[1]/2)
mysample$y <- jitter(mysample$y,resolution[1]/2)

write.csv(mysample,file="ValidationSample.csv",row.names=F)