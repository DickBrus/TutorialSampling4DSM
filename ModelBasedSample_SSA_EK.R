library(spcosa)
library(sp)
library(matrixcalc)

# Annealing functions
source("Functions4SSA.R")

#Read data with coordinates and other attributes of fine grid (discretization of study area)

load(file="HunterValley4Practicals.RData")
grd <- grdHunterValley
coordinates(grd)<- ~Easting+Northing
gridded(grd)<-TRUE

# Define correlation function
# Exponential
CorF <- function(D, pars) {
  A <- ifelse(D>.Machine$double.eps,pars[1]*exp(-D/pars[2]),1)
  A
}

s <- 0.8 #ratio of spatial dependence c1/(c0+c1)
range<-200
thetas <- c(s,range)

# Select spatial coverage sample for prediction. These locations are fixed, i.e. their locations are not opimized in simulated annealing
# Note that a spatial coverage sample is not strictly needed! The alternative is to to optimize the coordinates of all locations in SSA

# Choose number of locations of spatial coverage sample
n <-90
set.seed(314)
myStrata <- stratify(grd, nStrata = n, equalArea=FALSE, nTry=10)
mySCsample <- as(spsample(myStrata),"SpatialPoints")

# Select initial supplemental sample
nsup <- 10
ids <- sample.int(nrow(grd),nsup)
mysupsample <- as(grd[ids,],"SpatialPoints")

# Select evaluation sample
myevalsample<-spsample(x=grd,n=100,type="regular",offset=c(0.5,0.5))

# Set amount of perturbation of correlogram parameters
perturbation <- 0.01

# Choose one of the following minimization criterions:
# AV: Augmented kriging Variance, see Eq. 5 in Lark and Marchant (2018), Geoderma
# EAC: Estimation Adjusted Criterion, see Eq. 2.16 in Zhu and Stein (2006), JABES

annealingResult <- anneal.EK(
  free = mysupsample,
  disc = grd,
  fixed = mySCsample,
  eval = myevalsample,
  thetas = thetas,
  perturbation=perturbation,
  criterion="EAC",
#  initialTemperature = 0.002, #AV
  initialTemperature = 0.005, #EAC
  coolingRate = 0.8,
  maxAccepted = 5*nrow(coordinates(mysupsample)),
  maxPermuted = 5*nrow(coordinates(mysupsample)),
  maxNoChange = 5,
verbose = "TRUE"
)

save(annealingResult,file="MBSample_EAC_phi200nug02_HunterValley.RData")

crit <- annealingResult$Criterion
ggplot() +
  geom_line(mapping=aes(x=1:length(crit),y=crit),colour="red")+
  scale_x_continuous(name="Chain")+
  scale_y_continuous(name="Minimization criterion")

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

plot(myStrata) +
  geom_point(data = mySCsampledf, mapping = aes(x= Easting, y =Northing), shape =1,size=1 )+
  geom_point(data = mysupsampledf, mapping = aes(x= Easting, y =Northing), shape =2,size=1 )
dev.off()
