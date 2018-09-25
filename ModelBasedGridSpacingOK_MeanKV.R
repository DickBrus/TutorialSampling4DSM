library(gstat)
library(ggplot2)
library(sp)

# Load data on soil organic carbon
load(file="Data/DataThreeWoredasEthiopia.RData")
head(priordataEthiopia)

# Load file with discretisation grid
load("Data/CovariatesThreeWoredasEthiopia.RData")
head(grdEthiopia)

# Estimate experimental variogram

vg <- variogram(SOC~1, data = priordataEthiopia)
plot(vg)

# Fit variogram model

vgfit <- fit.variogram(vg, model = vgm(psill = 0.6, model = "Sph", range = 40, nugget = 0.6))
plot(vg,vgfit)
print(vgfit)

# Select a simple random sample of size 1000 for evaluating the square grids. 
# Add a small number to the x-coordinates and y-coordinates by drawing from a uniform distribution with lower and upper limit equal to -cellsize/2 and +cellsize/2, respectively. This can be done with function jitter.

set.seed(314)
ids<-sample.int(nrow(grdEthiopia),size=1000)
mysample<-grdEthiopia[ids,]

# Shift the randomly selected points to random points within the cells (cellsize is 1 km x 1 km)

mysample$s1 <- jitter(mysample$s1,0.5)
mysample$s2 <- jitter(mysample$s2,0.5)

coordinates(mysample)<-~s1+s2
coordinates(grdEthiopia)<-~s1+s2
gridded(grdEthiopia)<-TRUE

#Define grid spacings

spacing<-seq(from=5,to=12,by=1)

#Set number of times grid sampling of a given spacing is repeated

r<-10

MKV<-samplesize<-matrix(nrow=length(spacing),ncol=r)
for (i in 1:length(spacing)) {
  for (j in 1:r) {
    mygridxy<-spsample(x=grdEthiopia,cellsize=spacing[i],type="regular")
    mygrid<-data.frame(s1=mygridxy$x1,s2=mygridxy$x2,dummy=1)
    samplesize[i,j] <- nrow(mygrid)
    coordinates(mygrid)<-~s1+s2
    #Use gstat for ordinary kriging predictions
    predictions  <- krige(
        dummy ~ 1,
        mygrid,
        newdata = mysample,
        model = vgfit,
        nmax = 100
    )
    MKV[i,j]<-mean(predictions$var1.var)
  }
}

MMKV<-apply(MKV,MARGIN=1,FUN=mean)
Msize<-apply(samplesize,MARGIN=1,FUN=mean)

result<-data.frame(spacing,MKV,MMKV,samplesize,Msize)
#write.csv(result,file="MOKVvsGridspacing_Ethiopia.csv",row.names=F)


library(ggplot2)
ggplot(data=result)+
  geom_point(mapping=aes(x=spacing,y=MMKV),size=3)+
  scale_x_continuous(name="Spacing (km)")+
  scale_y_continuous(name="Mean Kriging Variance")

#Plot the mean kriging variance against the sample size.

ggplot(data=result)+
  geom_point(mapping=aes(x=Msize,y=MKV),size=3)+
  scale_x_continuous(name="Sample size")+
  scale_y_continuous(name="Mean Kriging Variance")