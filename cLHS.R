# NOTE THAT, CONTRARY TO WHAT IS STATED IN THE GEODERMA PAPER, SUPPLEMENTING 
# AN EXISTING SAMPLE (LEGACY SAMPLE) WITH A CONDITIONED LATIN HYPERCUBE SAMPLE
# CAN ALSO BE DONE WITH R PACKAGE SPSANN! SEE ARGUMENT points OF FUNCTION optimCLHS


library(spcosa)
library(sp)
library(ggplot2)

# Source annealing functions
source('Functions4SSA.R')

# Read grid with covariates
load(file="Data/HunterValley4Practicals.RData")
grd<-grdHunterValley
rm(grdHunterValley)

# In which columns are the coordinates and covariates?
col.xy <- c(1,2)
col.cov <- c(3,4,5,6,7)

# Compute population correlation matrix of covariates

R<-cor(grd[,col.cov])

# Typecast grd to SpatialPixelsDataFrame
grid <- SpatialPixelsDataFrame(
    points = grd[,col.xy],
    data   = grd[,col.cov]
)

# Select legacy sample
set.seed(314)
ids <- sample.int(nrow(grd),10)
legacy <- SpatialPoints(
  coords=grd[ids,col.xy]
)

# Select spatial infill sample which is used as initial sample in annealing
set.seed(314)
samplesize<-50 #number of additional points
ntot <- samplesize+length(legacy)

myStrata <- stratify(grid,nStrata = ntot, priorPoints=legacy, equalArea=FALSE, nTry=1)
mySample <- spsample(myStrata)
plot(myStrata, mySample)

# Select the new points from mySample
ids <- which(mySample@isPriorPoint==F)
mySample <- as(mySample, "SpatialPoints")
mySample <- mySample[ids,]

# Compute lower bounds of marginal strata
by=1/ntot
probs<-seq(from=0,to=1,by=by)
lb <- apply(grd[,3:7],MARGIN=2,FUN=function(x) quantile(x,probs=probs,type=7))
lb <- lb[-length(probs),]

#set relative weight of O1 for computing the LHS criterion (O1 is for coverage of marginal strata of covariates); 1-W01  is the relative weight for O3 (for correlation)
wO1<-0.5

#now start the annealing
system.time(
annealingResult <- anneal.cLHS(
    d = mySample,
    g = grid,
    legacy = legacy,
    lb = lb,
    wO1=wO1,
    R=R,
    initialTemperature = 0.01,
    coolingRate = 0.9,
    maxAccepted = 5*length(mySample),
    maxPermuted = 5*length(mySample),
    maxNoChange=10,
    verbose = "TRUE"
    )
)

save(annealingResult,file="LHSample_50(0.5).Rdata")
load(file="LHSample_50(0.5).Rdata")

optSample<-as(annealingResult$optSample, "data.frame")
Eall<-annealingResult$Criterion

#Plot the selected points on top of one of the covariates
legacy <- as(legacy,"data.frame")
#pdf(file = "LHSample_50(05)_cti.pdf", width = 7, height = 7)
ggplot(data=grd) +
  geom_raster(mapping = aes(x = Easting/1000, y = Northing/1000, fill = cti))+  
  geom_point(data = optSample, mapping = aes(x = Easting/1000, y = Northing/1000), colour = "black") +
  geom_point(data = legacy, mapping = aes(x = Easting/1000, y = Northing/1000), colour = "red") +
  scale_x_continuous(name = "Easting (km)") +
  scale_y_continuous(name = "Northing (km)") +    
  scale_fill_gradient(name="cti",low = "darkblue", high = "red")+
  coord_fixed()
#dev.off()

#Make scatter plots
coordinates(optSample)<-~Easting+Northing
optSample <- over(optSample,grid)                

coordinates(legacy)<-~Easting+Northing
legacy <- over(legacy,grid)                

ggplot(data=grd) +
  geom_point(mapping = aes(x = ndvi, y = cti), colour = "black",size=1,alpha=0.5) +
  geom_point(data=as.data.frame(optSample), mapping = aes(x = ndvi, y = cti), colour = "red",size=2) +
  geom_point(data=as.data.frame(legacy), mapping = aes(x = ndvi, y = cti), colour = "green",size=2) +
  scale_x_continuous(name = "ndvi") +
  scale_y_continuous(name = "cti")

#Plot marginal sample sizes

mySample <-rbind(legacy,optSample)
p<-length(col.cov)
stratum<-matrix(nrow=nrow(mySample),ncol=p)
for ( i in 1:p) {
  stratum[,i]<-findInterval(mySample[,i],lb[,i],left.open = TRUE)
}

counts<-matrix(nrow=nrow(lb),ncol=p)
for (i in 1:nrow(lb)) {
  counts[i,]<-apply(stratum, MARGIN=2, function(x,i) sum(x==i), i=i)
}

index<-seq(1:ntot)
countsdf<-as.data.frame(counts)
names(countsdf)<-names(grd[col.cov])
library(reshape)
countslf<-melt(countsdf)
countslf$index<-rep(index,times=5)
  ggplot(countslf) +
  geom_point(mapping = aes(x=index,y = value), colour = "black",size=1) +
  facet_wrap(~variable) +
  scale_x_continuous(name = "Index") +
  scale_y_continuous(name = "Sample size",breaks=c(0,1,2,3,4))
