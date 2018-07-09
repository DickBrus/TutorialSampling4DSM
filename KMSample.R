library(sp)
library(fields)
library(ggplot2)

#Read data with coordinates and other attributes of fine grid (discretization of study area)

load(file="Data/HunterValley4Practicals.RData")
summary(grdHunterValley)
cor(grdHunterValley[,c(3,4,5,6,7)])

#Set number of sampling locations to be selected

n<-20

#Compute clusters

set.seed(314)
myClusters <- kmeans(scale(grdHunterValley[,c(3,4,5,6,7)]), centers=n, iter.max=100,nstart=10)
grdHunterValley$clusters <- myClusters$cluster

#Select locations closest to the centers of the clusters

rdist.out <- rdist(x1=myClusters$centers,x2=scale(grdHunterValley[,c(3,4,5,6,7)]))
ids.mindist <- apply(rdist.out,MARGIN=1,which.min)
mySample <- grdHunterValley[ids.mindist,]

#Plot clusters and sampling points

pdf(file = "KMSample_HunterValley.pdf", width = 7, height = 7)
ggplot(grdHunterValley) +
  geom_tile(mapping = aes(x = Easting, y = Northing, fill = factor(clusters))) +
  scale_fill_discrete(name = "cluster") +
  geom_point(data=mySample,mapping=aes(x=Easting,y=Northing),size=2) +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "") +
  coord_fixed() +
  theme(legend.position="none")
dev.off()

pdf(file = "Scatterplot_KMSample_HunterValley.pdf", width = 7, height = 7)
ggplot(grdHunterValley) +
  geom_point(mapping=aes(y=elevation_m,x=cti,colour=factor(clusters))) +
  geom_point(data=mySample,mapping=aes(y=elevation_m,x=cti),size=2) +
  scale_y_continuous(name = "Elevation") +
  scale_x_continuous(name = "CTI") +
  theme(legend.position="none")
dev.off()
