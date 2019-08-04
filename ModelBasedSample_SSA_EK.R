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

s <- 0.8 #ratio of spatial dependence c1/(c0+c1)
range<-200
thetas <- c(s,range)


## SAMPLING FOR VARIOGRM ESTIMATION USING EITHER LOGDET OR VV AS CRITERION


# Select initial sample
set.seed(314)
n <- 100
ids <- sample.int(nrow(grd),n)
mysample0 <- as(grd[ids,],"SpatialPoints")

# Define sampling grid used for prediction after the second sampling phase 
xgrid <- c(0,1,2,3)
ygrid <- xgrid
grid <- expand.grid(xgrid,ygrid)
names(grid) <- c("x","y")
spacing<-300
grid$x <- grid$x*spacing
grid$y <- grid$y*spacing
coordinates(grid) <- ~x+y

# Compute evaluation point
myevalsample <- data.frame(x=mean(grid$x),y=mean(grid$y))
coordinates(myevalsample) <- ~x+y

# Set amount of perturbation of correlogram parameters
perturbation <- 0.01

# For variogram estimation choose one of the following minimization criterions:

# logdet: log of the determinant of inverse of Fisher information matrix
# VV: Variance of kriging Variance, see Eq. 9 in Lark (2002) Geoderma.

annealingResult <- anneal.EK(
  free = mysample0,
  disc = grd,
  fixed = grid,
  esample = myevalsample,
  model = "Exp",
  thetas=thetas,
  perturbation=perturbation,
  criterion="VV",
  #  initialTemperature = 0.05, #logdet
  initialTemperature = 0.00005, #VV
  coolingRate = 0.8,
  maxAccepted = 5*nrow(coordinates(mysample0)),
  maxPermuted = 5*nrow(coordinates(mysample0)),
  maxNoChange = 5,
  verbose = "TRUE"
)

save(annealingResult,file="MBSample_VV_phi200nug05_HunterValley.RData")

#load(file="MBSample_logdet_phi200nug02_HunterValley.RData")
load(file="MBSample_VV_phi200nug02_HunterValley.RData")


library(ggplot2)
mysampledf <- data.frame(annealingResult$optSample)

# Plot sample

#pdf(file = "MB_logdet_phi200nug05_HunterValley.pdf", width = 7, height = 7)
ggplot(grdHunterValley) +
  geom_raster(mapping = aes(x= Easting, y =Northing),fill="grey")+
  geom_point(data = mysampledf, mapping = aes(x= Easting, y =Northing), shape =2,size=1 )+
  coord_fixed() 
#dev.off()

crit <- annealingResult$Criterion
ggplot() +
  geom_line(mapping=aes(x=1:length(crit),y=crit),colour="red")+
  scale_x_continuous(name="Chain")+
  scale_y_continuous(name="Minimization criterion")





## SAMPLING FOR VARIOGRAM ESTIMATION AND PREDICTION USING EITHER AV OR EAC AS CRITERION






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
  esample = myevalsample,
  model = "Exp",
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
