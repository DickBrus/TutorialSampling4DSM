library(sp)
library(matrixcalc)

# Annealing functions
source("Functions4SSA.R")

#Read data with coordinates and other attributes of fine grid (discretization of study area)

load(file="HunterValley4Practicals.RData")
grd <- grdHunterValley
gridded(grd)<- ~Easting+Northing

s <- 0.5 #ratio of spatial dependence c1/(c0+c1)
phi<-200
thetas <- c(s,phi)

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
eval <- data.frame(x=mean(grid$x),y=mean(grid$y))
coordinates(eval) <- ~x+y

# Set amount of perturbation of correlogram parameters
perturbation <- 0.01

# For variogram estimation choose one of the following minimization criterions:

# logdet: log of the determinant of inverse of Fisher information matrix
# VV: Variance of kriging Variance, see Eq. 9 in Lark (2002) Geoderma.

annealingResult <- anneal.EK(
  free = mysample0,
  disc = grd,
  fixed = grid,
  eval = eval,
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
