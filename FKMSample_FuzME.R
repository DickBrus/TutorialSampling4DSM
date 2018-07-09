# read data
load(file="Data/HunterValley4Practicals.RData")

# read memberships computed with FuzMe
m <- read.csv(file="Data/20_class.txt",sep="")
m <- m[,-c(1,2,3)]

# defuzzify, i.e. compute for each gridcell the cluster wth largest membership
grdHunterValley$cluster <- apply(m,MARGIN=1,which.max)

n <- ncol(m)-3

#select locations with largest membership in cluster 1...k
units <- apply(m,MARGIN=2,FUN=which.max)
myFKMSample <- grdHunterValley[units,]

#plot clusters and sampling points
#pdf(file = "FKMSample_phi13_HunterValley.pdf", width = 7, height = 7)
ggplot(grdHunterValley) +
  geom_raster(mapping = aes(x = Easting, y = Northing, fill = factor(cluster))) +
  scale_fill_discrete(name = "cluster") +
  geom_point(data=myFKMSample,mapping=aes(x=Easting,y=Northing),size=2) +
  scale_x_continuous(name = "") +
  scale_y_continuous(name = "") +
  coord_fixed() +
  theme(legend.position="none")
#dev.off()

#pdf(file = "Scatterplot_FKMSample_phi13_HunterValley.pdf", width = 7, height = 7)
ggplot(grdHunterValley) +
    geom_point(mapping=aes(y=elevation_m,x=cti,colour=factor(cluster))) +
  geom_point(data=myFKMSample,mapping=aes(y=elevation_m,x=cti),size=2) +
  scale_y_continuous(name = "Elevation") +
  scale_x_continuous(name = "CTI") +
  theme(legend.position="none")
#dev.off()
