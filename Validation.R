library(ggplot2)

#read validation sample data
d<-read.csv(file="Data/Stratifiedrandomsample.csv",header=T)

#compute squared errors
d$e2KED <- (d$SOM_A_hori - d$KED)^2
d$e2RF <- (d$SOM_A_hori - d$RF)^2


#compute stratum sample sizes
(nh<-tapply(d$SOM_A_hori,INDEX=d$stratum,FUN=length))

#read sizes of strata
strata<-read.csv(file="Data/StrataSize.csv",header=T)
names(strata)[1] <- "stratum"

#compute stratum weights
Nh <- strata$Number.of.Pixels
wh <- Nh/sum(Nh)

#estimate MSE
me2KED_h<-tapply(d$e2KED,INDEX=d$stratum,FUN=mean)
(mseKED<-sum(wh*me2KED_h))

me2RF_h<-tapply(d$e2RF,INDEX=d$stratum,FUN=mean)
(mseRF<-sum(wh*me2RF_h))

#estimate standard error of estimated mean squared error

#Note that in stratum 2 there is only 1 point
#A solution is to collapse strata 2 and 1 (strata 1 and 2 are similar geological units)

#collapse strata with 1 point only
levels<-sort(unique(d$stratum))
collapsedstrata<-c(1,1,2,3,4,5,6,7)
lut <- data.frame(stratum = levels, collapsedstrata)
d <- merge(x = d, y = lut)
(strata <- merge(x = strata, y = lut))

#estimate sample sizes and total number of pixels per collapsed stratum
nhc<-tapply(d$SOM_A_hori,INDEX=d$collapsedstrata,FUN=length)
Nhc<-tapply(strata$Number.of.Pixels,INDEX=strata$collapsedstrata,FUN=sum)

#compute collapsed stratum weights and estimate standard errors
whc<-Nhc/sum(Nhc)

vare2KED_h<-tapply(d$e2KED,INDEX=d$collapsedstrata,FUN=var)
(sdmseKED<-sqrt(sum(whc^2*vare2KED_h/nhc)))

vare2RF_h<-tapply(d$e2RF,INDEX=d$collapsedstrata,FUN=var)
(sdmseRF<-sqrt(sum(whc^2*vare2RF_h/nhc)))