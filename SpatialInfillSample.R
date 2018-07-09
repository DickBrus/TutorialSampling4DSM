# Loading R packages

library(sp)
library(spcosa)
library(ggplot2)

#Load data.frame with coordinates (and other attributes) of fine grid (discretization of study area)
load("Data/CovariatesThreeWoredasEthiopia.RData")

#Load existing sampling points
load(file="Data/DataThreeWoredasEthiopia.RData")

#Plot prior points
ggplot(data = as.data.frame(priordataEthiopia)) +
    geom_raster(data=grdEthiopia,mapping = aes(x = s1,y = s2),fill="grey") +
    geom_point(mapping = aes(x = s1,y = s2),size = 2) +
    scale_x_continuous(name = "Easting (km)") +
    scale_y_continuous(name = "Northing (km)") +
    coord_equal(ratio = 1)

#Construct compact geostrata for spatial infill sampling of 100 points, and select an infill sample of 100 points.

#Change class of grdEthiopia from data.frame to SpatialPixelsDataFrame
coordinates(grdEthiopia)<-~s1+s2
gridded(grdEthiopia)<-TRUE

#Set number of new sampling locations to be selected
n<-100

#Compute total sample size (existing points + new points)
ntot<-n+length(priordataEthiopia)

#Change class of d (existing points) from SpatialPointsDataFrame to SpatialPoints
priordataEthiopia<-as(priordataEthiopia,"SpatialPoints")

#grdEthiopia does not have projection attributes, whereas d does. Remove projection attributes of d
proj4string(grdEthiopia)
proj4string(priordataEthiopia)<- NA_character_

#Compute geostrata with option priorPoints=priordataEthiopia
set.seed(314)
myStrata <- stratify(grdEthiopia, nStrata = ntot, priorPoints=priordataEthiopia, nTry=10)

#Select sampling points of infill sample (centres of geostrata)
mySample <- spsample(myStrata)

#Plot geostrata and sampling points (centres of geostrata)
plot(myStrata, mySample)


#In mySample both prior points (listed first) and new points are listed. Whether a point is prior or new is stored in the slot isPriorPoint. Determine how many new points are sampled. 

#Select the new points from mySample
ids <- which(mySample@isPriorPoint==F)

#Change class of mySample to data.frame
mySample <- as(mySample,"data.frame")
mySamplenew <- mySample[ids,]