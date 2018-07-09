library(spcosa)
library(gstat)
library(sp)

# Source annealing functions
source("Functions4SSA.R")

# Load data.frame with coordinates and other attributes of fine grid (discretization of study area)
load(file="Data/CovariatesThreeWoredasEthiopia.RData")
covariates <- c("dem","rfl_NIR","rfl_red","lst")
grid <-grdEthiopia[,c("s1","s2",covariates)]
coordinates(grid) <- ~s1+s2
gridded(grid) <- T

# Load existing sampling points
load(file="Data/DataThreeWoredasEthiopia.RData")

# Subsample legacy sample (cellsize of grid = 1 km x 1 km)
legacy <- remove.duplicates(priordataEthiopia,zero=1,remove.second=T)
legacy <- as(legacy,"data.frame")
legacy <- legacy[,c("s1","s2",covariates)]
coordinates(legacy) <- ~s1+s2

# Create prediction grid 'p'
p <- spsample(x = grid, n = 1000, type = "regular", offset = c(0.5, 0.5))
overgrd <- (p %over%grid)
p <- data.frame(coordinates(p),overgrd)
names(p)[c(1,2)] <- c("s1","s2")
coordinates(p) <- ~s1+s2

# Set number of new sampling locations to be selected
n<-100

# Select spatial coverage sample as initial infill sample
# Compute total sample size (existing points + new points)
ntot<-n+length(legacy)

myStrata <- stratify(grid,nStrata = ntot, priorPoints=as(legacy,"SpatialPoints"),equalArea=FALSE, nTry=1)
mySample <- spsample(myStrata)
plot(myStrata, mySample)

# Select the new points from mySample
ids <- which(mySample@isPriorPoint==F)

# Change class of mySample
infill <- as(mySample, "SpatialPoints")
infill <- infill[ids,]

# Overlay with grid
overgrd <-  (infill%over%grid)
infill <- data.frame(coordinates(infill),overgrd)
coordinates(infill) <- ~s1+s2

# Estimate the variogram from legacy sample
vg <- variogram(SOC~dem+rfl_NIR+rfl_red+lst, data = priordataEthiopia, cutoff=20)
vgfitOLS <- fit.variogram(vg, model = vgm(model = "Sph", psill = 0.2, range = 6, nugget = 0.5))
plot(vg,vgfitOLS)

# Fit model by REML
library(geoR)
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

# Set variogram parameters to those estimated by REML
vgfitREML <- vgfitOLS
vgfitREML[1,2] <- lmSOC_REML$nugget
vgfitREML[2,2] <- lmSOC_REML$sigmasq
vgfitREML[2,3] <- lmSOC_REML$phi

# Start the optimization
annealingResult <- anneal.K(
    d = infill,
    g = grid,
    p = p,
    legacy=legacy,
    model=vgfitREML,
    nmax=40,
    initialTemperature = 0.0005,
    coolingRate = 0.9,
    maxAccepted = 2*nrow(coordinates(infill)),
    maxPermuted = 2*nrow(coordinates(infill)),
    maxNoChange = 2*nrow(coordinates(infill)),
    verbose = "TRUE"
    )

save(annealingResult,file="ModelBasedSample_KED_Ethiopia.RData")

library(ggplot2)
infillSampledf <-data.frame(annealingResult$optSample)
legacy <- as(legacy,"data.frame")

pdf(file = "ModelBasedInfillSample_KED_Ethiopia.pdf", width = 5, height = 5)
ggplot() +
  geom_tile(grdEthiopia,mapping=aes(x=s1,y=s2,fill=rfl_NIR))+
  geom_point(data = infillSampledf, mapping = aes(x = s1, y = s2)) +
  geom_point(data = legacy, mapping = aes(x = s1, y = s2), shape=2) +
  scale_fill_continuous(low="darkblue",high="red",name="NIR")+
  coord_fixed()
dev.off()

traceMOKV <- annealingResult$Criterion
pdf(file = "TraceMKEDV_Ethiopia.pdf", width = 7, height = 7)
ggplot() +
  geom_line(mapping=aes(x=1:length(traceMOKV),y=traceMOKV),colour="red")+
  scale_x_continuous(name="Chain")+
  scale_y_continuous(name="Mean ordinary kriging variance")
dev.off()
