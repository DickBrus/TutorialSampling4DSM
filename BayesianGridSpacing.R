library(BayesianTools)
library(sp)
library(gstat) #for MoM estimation of variogram
library(ggplot2)

exponential <- function(D, pars) {
  C <- pars[1]*pars[2]*(exp(-D/pars[3]))
  diag(C) <- pars[1]
  C
}

spherical <- function(D, pars) {
  C <- ifelse(D<pars[3],pars[1]*pars[2]*(1 - (1.5*D/pars[3]) + 0.5*(D/pars[3])^3),0)
  diag(C) <- pars[1]
  C
}
# Define loglikelihood
ll <-function(pars) {
  C<-spherical(D,pars)
  cholC <- try(chol(C),silent=TRUE)
  if (!is.character(cholC)){
    C_inv <- chol2inv(chol(C))
    ### estimate the unknown mean
    XC <- crossprod(X, C_inv)
    XCy <- XC %*% y
    XCX <- XC %*% X
    betahat <- solve(XCX , XCy)
    mu <- as.numeric(t(betahat)%*%t(X))
    e <- y - mu
    eC <- crossprod(e, C_inv)
    eCe <- eC %*% e
    logDetC <- determinant(x = C, logarithm = TRUE)$modulus
    logLik<- as.numeric(-0.5 * (logDetC + eCe +length(y)*log(2*pi)))
    return(logLik)
  }else{return(-Inf)}
}

# Load existing sampling points
load(file="Data/DataThreeWoredasEthiopia.RData")

# compute matrix with Euclidian distances between sampling points
D <- spDists(priordataEthiopia)

# Estimate variogram by method-of-moments
vg <- variogram(SOC~1,data=priordataEthiopia)
plot(vg)

# vgfit <- fit.variogram(vg,model=vgm("Exp",psill=0.6,range=30,nugget=0.6))
vgfit <- fit.variogram(vg,model=vgm("Sph",psill=0.6,range=30,nugget=0.6))
plot(vg,vgfit)
print(vgfit)

# estimate variogram by maximum likelihood. Use MoM variogram parameters as initial values
sigma2.ini <- vgfit$psill[1]+vgfit$psill[2]
xi.ini <- vgfit$psill[2]/sigma2.ini #proportion spatially structured variance 
phi.ini <- vgfit$range[2]
pars <- c(sigma2.ini,xi.ini,phi.ini)

X <- matrix(1 , nrow(d) , 1)
y <- d$SOC
vgML <- optim(pars, ll, control = list(fnscale = -1), hessian = TRUE)
#invfisher_info <- solve(-vgML$hessian)
#diag(invfisher_info)

sigma2.ini <- vgML$par[1]
xi.ini <- vgML$par[2]
phi.ini <- vgML$par[3]

priors <- createUniformPrior(lower = c(0,0,0), upper = c(5,1,100))
setup <- createBayesianSetup(likelihood=ll,prior=priors,best=c(sigma2.ini,xi.ini,phi.ini), names=c("sigma2","xi","phi"))

set.seed(314)

DEzs.out <- runMCMC(setup,sampler="DEzs")
save(DEzs.out,file="DEzs_Ethiopia.RData")

#load(file="DEzs_Ethiopia.RData")
summary(DEzs.out)
plot(DEzs.out)
correlationPlot(DEzs.out)
marginalPlot(DEzs.out)
mcmcsam <- getSample(DEzs.out,start=1000,numSamples=1000)
mcmcsample <-data.frame(mcmcsam)

# Select a simple random sample of size 1000 for evaluating the square grids.
# Add a small number to the x-coordinates and y-coordinates by drawing from a uniform distribution with lower and upper limit equal to -cellsize/2 and +cellsize/2, respectively. This can be done with function jitter.

# Load file with discretisation grid
load("Data/CovariatesThreeWoredasEthiopia.RData")
set.seed(314)
ids<-sample.int(nrow(grdEthiopia),size=1000)
mysample<-grdEthiopia[ids,]

# Shift the randomly selected points to random points within the cells (cellsize is 1 km x 1 km)

mysample$s1 <- jitter(mysample$s1,0.5)
mysample$s2 <- jitter(mysample$s2,0.5)

coordinates(mysample)<-~s1+s2
coordinates(grdEthiopia)<-~s1+s2
gridded(grdEthiopia)<-TRUE

# Define grid spacings
spacing<-seq(from=1,to=12,by=1)

MKV<-matrix(nrow=length(spacing),ncol=nrow(mcmcsample))
for (i in 1:length(spacing)) {
    mygridxy<-spsample(x=grdEthiopia,cellsize=spacing[i],type="regular",offset=c(0.5,0.5))
    mygrid<-data.frame(s1=mygridxy$x1,s2=mygridxy$x2,dummy=1)
    coordinates(mygrid)<-~s1+s2
  for (j in 1:nrow(mcmcsample)) {
  #Use gstat for ordinary kriging predictions
    vgfit$psill[1] <- (1-mcmcsample$xi[j])*mcmcsample$sigma2[j]
    vgfit$psill[2] <- mcmcsample$xi[j]*mcmcsample$sigma2[j]
    vgfit$range[2] <- mcmcsample$phi[j]    
    predictions  <- krige(
      dummy ~ 1,
      mygrid,
      newdata = mysample,
      model = vgfit,
      nmax = 100
    )
    MKV[i,j]<-mean(predictions$var1.var)
  }
}

save(MKV,file="MOKV_Bayesian_Ethiopia.RData")

(MMKV <- apply(MKV,MARGIN=1,FUN=mean))

# Set target for MKV
t.MKV <- 0.8

# Compute for each variogram parameter vector the spacing for which the MKV equals the target variance
spacing.tol <- numeric(length=ncol(MKV))
for (i in 1:ncol(MKV)) {
  spacing.tol[i] <- approx(x=MKV[,i],y=spacing,xout=t.MKV)$y
}
summary(spacing.tol)
pdf("Histogram_TolerableGridSpacing_Ethiopia.pdf",width=5,height=5)
ggplot()+
  geom_histogram(mapping=aes(spacing.tol),binwidth=1,colour="orange")+
  scale_x_continuous(name="Tolerable grid spacing",breaks=seq(2:12))+
  scale_y_continuous(name="# MCMC samples")
dev.off()


# Compute for each grid spacing the porportion of MCMC samples with MKV smaller or equal to target MKV
F <- numeric(length=length(spacing))
for (i in 1:length(spacing)) {
  F[i] <- sum(MKV[i,] < t.MKV)
}
F <- F/ncol(MKV)
df <- data.frame(spacing,F)
pdf("ProportionvsGridspacing_Ethiopia.pdf",width=5,height=5)
ggplot(df)+
  geom_line(mapping=aes(x=spacing,y=F),se=FALSE,colour="red")+
  scale_x_continuous(breaks=spacing,name="Grid spacing")+
  scale_y_continuous(limits=c(0,1),breaks=seq(from=0,to=1,by=0.2),name="Proportion of MCMC samples")
dev.off()
