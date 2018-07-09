library(sp)
library(spcosa)

#Load data

load("Data/CovariatesThreeWoredasEthiopia.RData")

#Change class into SpatialPointsDataFrame and next into SpatialPixelsDataFrame

coordinates(grdEthiopia)<-~s1+s2
gridded(grdEthiopia)<-TRUE

#Set number of sampling locations to be selected

n<-100

#Compute compact geostrata

set.seed(314)
myStrata <- stratify(grdEthiopia, nStrata = n, nTry=10)
plot(myStrata)
getObjectiveFunctionValue(myStrata)

#Select the centres of the geostrata with function spsample

mySample <- spsample(myStrata)

#Plot geostrata and sampling points (centres of geostrata)

plot(myStrata, mySample)

#Write sampling points to csv file

#first change class of mySample to data.frame
mySample<-as(mySample,"data.frame")
write.csv(mySample,file="mySampleEthiopia.csv")