library(spcosa)
library(gstat)
library(sp)

# Source annealing functions
source("Functions4SSA.R")

# Load data.frame with coordinates (and other attributes) of fine grid (discretization of study area)
load(file="Data/CovariatesThreeWoredasEthiopia.RData")
grid <-grdEthiopia
coordinates(grid) <- ~s1+s2
gridded(grid) <- T

# Load existing sampling points
load(file="Data/DataThreeWoredasEthiopia.RData")

# Subsample legacy sample (cellsize of grid = 1 km x 1 km)
legacy <- remove.duplicates(d,zero=1,remove.second=T)
legacy <- as(legacy,"SpatialPoints")


# Create prediction grid 'p'
p <- spsample(x = grid, n = 1000, type = "regular", offset = c(0.5, 0.5))

# set number of new sampling locations to be selected
n<-100

# Compute total sample size (existing points + new points)
ntot<-n+length(legacy)

# Grid does not have projection attributes, whereas legacy does. Remove projection attributes of legacy
proj4string(legacy)<- NA_character_

myStrata <- stratify(grid,nStrata = ntot, priorPoints=legacy,equalArea=FALSE, nTry=1)
mySample <- spsample(myStrata)
plot(myStrata, mySample)

# Select the new points from mySample
ids <- which(mySample@isPriorPoint==F)

# Change class of mySample
infill <- as(mySample, "SpatialPoints")
infill <- infill[ids,]

# Estimate the variogram from legacy sample
vg <- variogram(SOC~1,d)
plot(vg)
vgmfit <- fit.variogram(vg,model=vgm(psill=0.6, "Sph", range=40,nugget=0.6))
print(vgmfit)

# Start the optimization
annealingResult <- anneal.K(
    d = infill,
    g = grid,
    p = p,
    legacy=legacy,
    model=vgmfit,
    nmax=20,
    initialTemperature = 0.0005,
    coolingRate = 0.9,
    maxAccepted = 2*nrow(coordinates(infill)),
    maxPermuted = 2*nrow(coordinates(infill)),
    maxNoChange = 2*nrow(coordinates(infill)),
    verbose = "TRUE"
    )

save(annealingResult,file="ModelBasedSample_OK_Ethiopia.RData")
load("d:/UserData/Sampling4DSM/Rscripts/ModelBasedSample_OK_Ethiopia.RData")

library(ggplot2)
infillSampledf <-data.frame(annealingResult$optSample)
legacy <- as(legacy,"data.frame")

pdf(file = "ModelBasedInfillSample_OK_Ethiopia.pdf", width = 5, height = 5)
ggplot() +
  geom_tile(grdEthiopia,mapping=aes(x=s1,y=s2),fill="grey")+
  geom_point(data = infillSampledf, mapping = aes(x = s1, y = s2)) +
  geom_point(data = legacy, mapping = aes(x = s1, y = s2), shape=2) +
  coord_fixed()
dev.off()

traceMOKV <- annealingResult$Criterion
pdf(file = "TraceMOKV_Ethiopia.pdf", width = 7, height = 7)
ggplot() +
  geom_line(mapping=aes(x=1:length(traceMOKV),y=traceMOKV),colour="red")+
  scale_x_continuous(name="Chain")+
  scale_y_continuous(name="Mean ordinary kriging variance")
dev.off()
