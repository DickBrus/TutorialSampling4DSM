#Read validation sample data
d<-read.csv(file="Data/Stratifiedrandomsample.csv",header=T)

#Compute errors and squared errors
d$eKED <- d$SOM_A_hori - d$KED
d$eRF <- d$SOM_A_hori - d$RF
d$e2KED <- d$eKED^2
d$e2RF <- d$eRF^2

#Compute stratum sample sizes
(nh<-tapply(d$SOM_A_hori,INDEX=d$stratum,FUN=length))

#Read sizes of strata
strata<-read.csv(file="Data/StrataSize.csv",header=T)
names(strata)[1] <- "stratum"

#Compute stratum weights
Nh <- strata$Number.of.Pixels
wh <- Nh/sum(Nh)

#Estimate population Mean Error (ME) for RF
meRF_h<-tapply(d$eRF,INDEX=d$stratum,FUN=mean)
(meRF<-sum(wh*meRF_h))

#Estimate standard error of estimated mean error

#Note that in stratum 2 there is only 1 point
#A solution is to collapse strata 2 and 1 (strata 1 and 2 are similar geological units)

#Collapse strata with 1 point only
levels<-sort(unique(d$stratum))
collapsedstrata<-c(1,1,2,3,4,5,6,7)
lut <- data.frame(stratum = levels, collapsedstrata)
d <- merge(x = d, y = lut)
(strata <- merge(x = strata, y = lut))

#Estimate sample sizes and total number of pixels per collapsed stratum
nhc<-tapply(d$SOM_A_hori,INDEX=d$collapsedstrata,FUN=length)
Nhc<-tapply(strata$Number.of.Pixels,INDEX=strata$collapsedstrata,FUN=sum)

#Compute collapsed stratum weights and estimate standard error of estimated ME
whc<-Nhc/sum(Nhc)

vareRF_h<-tapply(d$eRF,INDEX=d$collapsedstrata,FUN=var)
(sdmeRF<-sqrt(sum(whc^2*vareRF_h/nhc)))

#t test of hypothesis ME = 0 (no bias, systematic error = 0)
t <- meRF/sdmeRF
lowertail <- (t<0)
df <- nrow(d) - length(unique(d$collapsedstrata))
(p <- 2*pt(t,df=df,lower.tail=lowertail))


#Paired t test of hypothesis MSE(KED) = MSE(RF)

#Compute for each validation point differences of squared errors (paired differences)
d$dife2 <- d$e2KED-d$e2RF
mdife2_h<-tapply(d$dife2,INDEX=d$stratum,FUN=mean)
(mdife2<-sum(wh*mdife2_h))

vardife2_h<-tapply(d$dife2,INDEX=d$collapsedstrata,FUN=var)
(sdmdife2<-sqrt(sum(whc^2*vardife2_h/nhc)))

t <- mdife2/sdmdife2
lowertail <- (t<0)
df <- nrow(d) - length(unique(d$collapsedstrata))
(p <- 2*pt(t,df=df,lower.tail=lowertail))