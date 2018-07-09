library(gstat)
library(sp)
library(geoR)
library(ggplot2)

#Load data on soil organic carbon

load(file="Data/DataThreeWoredasEthiopia.RData")

#Load file with discretisation grid

load("Data/CovariatesThreeWoredasEthiopia.RData")

#Estimate variogram parameters by REML
#We assume that the mean of SOC is a linear combination of several covariates such as elevation (dem), enhanced vegetation index (evi), et cetera. To optimize the grid spacing we need the variogram of the residuals. Below this variogram is estimated by REML. Initial values are obtained by the parameters of the OLS residual variogram.

#Fit linear model, compute residuals, and compute experimental variogram of residuals

vg <- variogram(SOC~dem+rfl_NIR+rfl_red+lst, data = priordataEthiopia, cutoff=20)

vgfitOLS <- fit.variogram(vg, model = vgm(model = "Sph", psill = 0.2, range = 6, nugget = 0.5))
plot(vg,vgfitOLS)

#Fit model by REML
dGeoR <- as.geodata(
    obj = as.data.frame(priordataEthiopia), 
    header=TRUE, 
    coords.col=13:14, 
    data.col=1,
    data.names=NULL, 
    covar.col=c(3,6,9,10,11,13,14)
)

lmSOC_REML<- likfit(geodata = dGeoR, trend=~dem+rfl_NIR+rfl_red+lst,cov.model="spherical", ini.cov.pars=c(vgfitOLS[2,2], vgfitOLS[2,3]),nugget=vgfitOLS[1,2], lik.method="REML")
summary(lmSOC_REML)

#Select a simple random sample of size 1000 for evaluating the square grids. Add a small number to the x-coordinates and y-coordinates by drawing from a uniform distribution with lower and upper limit equal to -cellsize/2 and +cellsize/2, respectively.

set.seed(314)
ids<-sample.int(nrow(grdEthiopia),size=1000)
mysample<-grdEthiopia[ids,]

#Shift the randomly selected grid points to random points within the cells (the cellsize is 1 km x 1 km)

mysample$s1 <- jitter(mysample$s1,0.5)
mysample$s2 <- jitter(mysample$s2,0.5)

coordinates(mysample) <- ~s1+s2

#Set variogram parameters to those estimated by REML

vgfitREML <- vgfitOLS
vgfitREML[1,2] <- lmSOC_REML$nugget
vgfitREML[2,2] <- lmSOC_REML$sigmasq
vgfitREML[2,3] <- lmSOC_REML$phi

#Specify grid spacings

spacing<-seq(from=5,to=12,by=1)

#Set number of times grid sampling of a given spacing is repeated

r<-10

coordinates(grdEthiopia) <- ~s1+s2
gridded(grdEthiopia)<-TRUE

MKV<-matrix(nrow=length(spacing),ncol=r)
for (i in 1:length(spacing)) {
    for (j in 1:r) {
        mygridxy<-spsample(x=grdEthiopia,cellsize=spacing[i],type="regular")
        #add a dummy variable for interpolation
        mygrid<-data.frame(s1=mygridxy$x1,s2=mygridxy$x2,dummy=1)
        coordinates(mygrid)<-~s1+s2
        mygrd<-data.frame(mygrid %over% grdEthiopia,mygrid)
        coordinates(mygrd)<-~s1+s2
        #Use gstat for KED predictions
        predictions  <- krige(
            dummy ~ dem+rfl_NIR+rfl_red+lst,
            mygrd,
            newdata = mysample,
#            model = vgfitREML,
            model = vgfitOLS,
            nmax = 100
        )
        MKV[i,j]<-mean(predictions$var1.var)
    }
}

#Plot the mean kriging variance against the grid spacing.

MMKV<-apply(MKV,MARGIN=1,FUN=mean)
result<-data.frame(spacing,MKV,MMKV)

ggplot(data=result)+
  geom_point(mapping=aes(x=spacing,y=MMKV),size=3)+
  scale_x_continuous(name="Spacing (km)")+
  scale_y_continuous(name="Mean Kriging Variance")

ggplot(data=result)+
    geom_point(mapping=aes(x=spacing,y=X1),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X2),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X3),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X4),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X5),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X6),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X7),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X8),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X9),size=1,alpha=0.5)+
    geom_point(mapping=aes(x=spacing,y=X10),size=1,alpha=0.5)+
    scale_x_continuous(name="Spacing (km)")+
    scale_y_continuous(name="Mean Kriging Variance")

write.csv(result,file="MKEDVvsGridspacing_Ethiopia.csv",row.names=F)

