library(spcosa)
library(sp)
library(matrixcalc)

# Annealing functions
source("Functions4SSA.R")

#Read data with coordinates and other attributes of fine grid (discretization of study area)

load(file="Data/HunterValley4Practicals.RData")
grd <- grdHunterValley
coordinates(grd)<- ~Easting+Northing
gridded(grd)<-TRUE

# Select spatial coverage sample for prediction. These locations are fixed, i.e. their locations are not opimized in simulated annealing
# Note that a spatial coverage sample is not strictly needed! The alternative is to to optimize the coordinates of all locations in SSA

# Choose number of locations of spatial coverage sample
#n <-90
#set.seed(314)
#myStrata <- stratify(grd, nStrata = n, equalArea=FALSE, nTry=10)
#mySCsample <- as(spsample(myStrata),"SpatialPoints")

# Define correlation function
# Exponential
CorF <- function(D, pars) {
  C <- ifelse(D>.Machine$double.eps,pars[1]*exp(-D/pars[2]),1)
  C
}

s <- 0.8 #ratio of spatial dependence c1/(c0+c1)
range<-500
pars <- c(s,range)

# Select initial supplemental sample
nsup <- 100
ids <- sample.int(nrow(grd),nsup)
mysupsample <- as(grd[ids,],"SpatialPoints")

# Select evaluation sample
myevalsample<-spsample(x=grd,n=100,type="regular",offset=c(0.5,0.5))

# Set amount of perturbation of correlogram parameters
perturbation <- 0.01

# For variogram estimation choose one of the following minimization criterions:

# logdet: log of the determinant of inverse of Fisher information matrix
# VV: Variance of kriging Variance, see Eq. 9 in Lark (2002) Geoderma.

# For estimation of variogram and kriging choose one of the following minimization criterions:

# AV: Augmented kriging Variance, see Eq. 5 in Lark and Marchant (2018), Geoderma
# EAC: Estimation Adjusted Criterion, see Eq. 2.16 in Zhu and Stein (2006), JABES

set.seed(314)
annealingResult <- anneal.EK(
  d = mysupsample,
  g = grd,
#  spcosam=mySCsample,
  p = myevalsample,
  pars=pars,
  perturbation=perturbation,
  criterion="logdet",
  initialTemperature = 0.05, #logdet
#  initialTemperature = 0.005, #VV
#  initialTemperature = 0.002, #AV
#  initialTemperature = 0.005, #EAC
  coolingRate = 0.8,
  maxAccepted = 5*nrow(coordinates(mysupsample)),
  maxPermuted = 5*nrow(coordinates(mysupsample)),
#  maxNoChange = nrow(coordinates(mysupsample)),
  maxNoChange = 5,
verbose = "TRUE"
)

save(annealingResult,file="MBSample_VV_phi500nug02_nospcosa_HunterValley.RData")

load(file="MBSample_logdet_phi500nug02_nospcosa_HunterValley.RData")
load(file="MBSample_VV_phi500nug02_nospcosa_HunterValley.RData")
load(file="MBSample_logdet_phi500nug02_HunterValley.RData")
load(file="MBSample_VV_phi500nug02_HunterValley.RData")
load(file="MBSample_AV_phi500nug02_HunterValley.RData")
load(file="MBSample_EAC_phi500nug02_HunterValley.RData")

# compute distance of supplemental sample to nearest spatial coverage sample point
D <- spDists(mySCsample,annealingResult$optSample)
(Dmin <- apply(D,MARGIN=2,FUN=min))
hist(Dmin)

D <- spDists(annealingResult$optSample)
diag(D) <- 1E1000
(Dmin <- apply(D,MARGIN=2,FUN=min))

library(ggplot2)
mysupsampledf <- data.frame(annealingResult$optSample)
mySCsampledf <- as(mySCsample,"data.frame")


# Plot strata, spatial coverage sample and supplemental sample

pdf(file = "MB_logdet_phi500nug02_HunterValley.pdf", width = 7, height = 7)
plot(myStrata) +
  geom_point(data = mySCsampledf, mapping = aes(x= Easting, y =Northing), shape =1,size=1 )+
  geom_point(data = mysupsampledf, mapping = aes(x= Easting, y =Northing), shape =2,size=1 )
dev.off()

pdf(file = "MB_logdet_phi500nug02_nospcosa_HunterValley.pdf", width = 7, height = 7)
ggplot(grdHunterValley) +
  geom_raster(mapping = aes(x= Easting, y =Northing),fill="grey")+
  geom_point(data = mysupsampledf, mapping = aes(x= Easting, y =Northing), shape =2,size=1 )+
  coord_fixed() 
dev.off()

crit <- annealingResult$Criterion
ggplot() +
  geom_line(mapping=aes(x=1:length(crit),y=crit),colour="red")+
  scale_x_continuous(name="Chain")+
  scale_y_continuous(name="Minimization criterion")
