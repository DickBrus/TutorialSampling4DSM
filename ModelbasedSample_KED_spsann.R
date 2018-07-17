library(sp)
library(gstat)
library(spsann)

#Load grid with EMdata; the nodes of this grid are the potential sampling locations 
grid <- read.csv(file="Data/EMUzbekistan_25m.csv")
coordinates(grid) <- ~x1+x2
gridded(grid)<-TRUE
boundary <- rgeos::gUnaryUnion(as(grid, "SpatialPolygons"))

# residual variogram
variogram <- vgm(nugget=0.1, psill=0.075, "Exp", range=100)

gridded(grid)<-FALSE
candi <- data.frame(x=coordinates(grid)[,1],y=coordinates(grid)[,2])
covars <- as(grid,"data.frame")
covars <- covars[,c(3,4,1)]
names(covars)[c(1,2)] <- c("x","y")

schedule <- scheduleSPSANN(initial.acceptance = 0.8,initial.temperature = 0.002,
                           temperature.decrease=0.95,
                           chains=500,
                           chain.length=5,
                           stopping=10,
                           x.min=25,y.min=25,
                           cellsize=25)

set.seed(321)
res <- optimMKV(
  points = 50, candi = candi, covars=covars, vgm = variogram,
  eqn = z ~ lnEM1m, plotit = TRUE, schedule = schedule,boundary=boundary)
save(res,file="ModelbasedSample_KED_Uzbekistan_spsann.RData")


#load(file="d:/UserData/SamplingBookandCourses/Figures/Rscripts/Sampling4Mapping/ModelbasedSample_KED_Uzbekistan_spsann.RData")

#candi <- data.frame(x=coordinates(grid)[,1],y=coordinates(grid)[,2])
#covars <- as(grid,"data.frame")
#covars <- covars[,c(2,3,1)]

sample<-candi[res$points$id,]
#add covariate to sample
ids <- as.integer(rownames(sample))
sample$EM <- covars[ids,3]

library(ggplot2)
pdf(file = "KEDSample_Uzbekistan.pdf", width = 6, height = 4)
ggplot(data = covars) +
  geom_raster(mapping = aes(x = x, y = y, fill = EM)) +
  geom_point(data = sample, mapping = aes(x = x, y = y), size=2,colour = "black") +
#  scale_fill_gradient(name="x",low = "skyblue", high = "darkblue") +
  scale_fill_continuous(name = "lnEMv1m", low = rgb(0, 0.2, 1), high = rgb(1, 0.2, 0)) +
  scale_y_continuous(name = "Northing") +
  scale_x_continuous(name = "Easting") +
  coord_equal()
dev.off()

pdf(file="TraceMKV_MBSample_SSA_OK_Square.pdf",width=6,height=4)
ggplot(res$objective$energy) +
  geom_line(mapping = aes(x=1:nrow(res$objective$energy),y = obj),colour="red") +
  scale_y_continuous(name = "Mean Kriging Variance",limits=c(22.5,25)) +
  scale_x_continuous(name = "Chain")
dev.off()

pdf(file = "HistogramEMUzbekistan_Sample.pdf", width = 4, height = 3)
ggplot(data = sample) +
  geom_histogram(mapping = aes(EM),breaks=seq(from=2.75,to=5,by=0.25),color="orange")
dev.off()

griddf <- as(grid,"data.frame")

pdf(file = "HistogramEMUzbekistan_Population.pdf", width = 4, height = 3)
ggplot(data = griddf) +
  geom_histogram(mapping = aes(EM),breaks=seq(from=2.75,to=5,by=0.25),color="orange")
dev.off()
