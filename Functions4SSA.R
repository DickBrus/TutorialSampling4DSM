#Author: Dick Brus and Dennis Walvoort, November 29, 2017


# Function for generating a series of spatial samples

permute<-function(d, g)  {
  # extract coordinates of observation points 'd' and grid cells 'g'
  s_d <- coordinates(d)
  s_g <- coordinates(g)
  
  # randomly select one location in 'd'
  i_d <- sample(x = seq_len(nrow(s_d)), size = 1)
  
  # compute squared Euclidean distances 'd2' between the selected location and all grid cells
  d2 <- (s_g[, 1] - s_d[i_d, 1])^2 +
    (s_g[, 2] - s_d[i_d, 2])^2
  
  # randomly select a grid cell with a probability inverse to squared distance (p ~ 1/distance^2)
  i_g <- sample(x = seq_len(nrow(s_g)), size = 1, prob = 1/(d2 + 1))
  
  # replace randomly selected location in actual sample (s_d[i_d, ]) by a new location within the randomly selected grid cell (g[i_g, ])
  gridTopology <- as(getGridTopology(g), "data.frame")
  s_d[i_d, ] <- s_g[i_g, ] + runif(n = 2, min = -0.5, max = 0.5) * gridTopology$cellsize
  
  # return result
  SpatialPoints(coords = s_d)
}


# Annealing function for OK and KED
#d: SpatialPoints (OK) or SpatialPointsDataFrame (KED) of sampling points
#g: SpatialPixelsDataFrame (discretisation of study area)
#p: SpatialPoints (OK) or SpatialPointsDataFrame (KED) of prediction points
#legacy: SpatialPoints (OK) or SpatialPointsDataFrame (KED) of legacy sample
#model: semivariogram (gstat object)
#nmax: maxumum number of sampling points used in kriging

anneal.K<-function(d, g, p, legacy, model, nmax = 50,
                 initialTemperature = 1, coolingRate = 0.9, maxAccepted = 10 * nrow(coordinates(d)),
                 maxPermuted=10* nrow(coordinates(d)), maxNoChange=nrow(coordinates(d)), verbose = getOption("verbose")) {
  if(!(class(d) %in% c("SpatialPoints","SpatialPointsDataFrame")))
    stop("Error: d must be SpatialPoints(DataFrame)")
  if(!(class(g) %in% c("SpatialPixels","SpatialPixelsDataFrame"))) 
    stop("Error: g must be SpatialPixels(DataFrame)")
  if(!(class(p) %in% c("SpatialPoints","SpatialPointsDataFrame")))
    stop("Error: p must be SpatialPoints(DataFrame)")
  
  stopifnot(is.na(proj4string(legacy)))

  # set initial temperature
  T <- initialTemperature
  
  # merge infill sample and legacy sample
  dall <- d
  if(!missing(legacy)){
    dall <- rbind(d,legacy)
  }
  
  # compute the criterion (mean kriging variance)
  E <- getCriterion.K(dall, p, model, nmax)
  
  # store criterion
  E_prv <- E
  
  # Define structure for storing time series of criterion
  Eall<-NULL

  # initialize number of zero changes of objective function
  nNoChange <-0
  
  # start cooling loop
  repeat{
    
    # initialize number of accepted configurations
    nAccepted <- 0
    
    # initialize number of permuted configurations
    nPermuted <- 0
    
    # initialize number of improved configurations
    nImproved <- 0
    
    # start permutation loop
    repeat {
      
      # increase the number of permutations
      nPermuted <- nPermuted + 1
      
      # propose new sample by making use of function permute
      d_p <- permute(d, g)
      
      # for KED overlay new sample with grid
      if(length(names(p))>0) {
      d_p <- SpatialPointsDataFrame(
        coords = d_p,
        data = d_p %over% g                
      )}
      
      #merge infill sample and legacy sample
      dall_p <- d_p
      if(!missing(legacy)){
        dall_p <- rbind(d_p,legacy)
      }
      
      # compute the criterion of this new sample by using function getCriterion
      E_p <- getCriterion.K(dall_p, p, model, nmax)
      
      # accept/reject proposal by means of Metropolis criterion
      dE <- E_p - E
      if (dE < 0) {
        nImproved <- nImproved + 1
        prob <- 1 # always accept improvements
      } else {
        prob <- exp(-dE / T) # use Boltzmann to judge if deteriorations should be accepted
      }
      u <- runif(n = 1) # draw uniform deviate
      if (u < prob) { # accept proposal
        nAccepted <- nAccepted + 1
        d <- d_p
        E <- E_p
      }
      # are conditions met to lower temperature?
      lowerTemperature <- (nPermuted == maxPermuted) |
        (nAccepted == maxAccepted)
      if (lowerTemperature) {
        if (nImproved==0)
        {nNoChange<-nNoChange+1}
        else
        {nNoChange<-0}
        Eall<-rbind(Eall,E)
        break  
      }
    }
    
    if (verbose) {
      cat(
        format(Sys.time()), "|",
        sprintf("T = %e  E = %e  permuted = %d  accepted = %d  improved = %d  acceptance rate = %f \n",
                T, E, nPermuted, nAccepted, nImproved, nAccepted / nPermuted)
      )
    }
    
    # check on convergence
    if (nNoChange == maxNoChange) {
        break
    }
    E_prv <- E
    
    # lower temperature
    T <- coolingRate * T
  }
  
  # return result
  list(
    optSample=d,Criterion=Eall
  )
}


getCriterion.K<-function(d,p,model,nmax) {
  
  # add dummy variable
  if(class(d)=="SpatialPoints") {
    d <- SpatialPointsDataFrame(
      coords = d,
      data = data.frame(dum = rep(1, times = length(d)))
    )
  } else {
    d$dum=1
  }
  
  if(length(names(p))>0) {
    formul <- as.formula(paste("dum", paste(names(p), collapse = "+"), sep = "~"))} else {
      formul <- as.formula(paste("dum", paste(1, collapse = "+"), sep = "~"))
    }
  
  # compute variance of prediction error
  result <- krige(
    formula=formul,
    locations = d,
    newdata = p,
    model = model,
    nmax=nmax,
    debug.level = 0
  )
  mean(result$var1.var)
}



# Annealing function for estimation of variogram and kriging

#spcosam: SpatialPoints of spatial coverage sample
anneal.EK<-function(d, g, spcosam, p, pars, perturbation=0.01, criterion,
                  initialTemperature = 1, coolingRate = 0.9, maxAccepted = 10 * nrow(coordinates(d)),
                  maxPermuted=10* nrow(coordinates(d)), maxNoChange=nrow(coordinates(d)), verbose = getOption("verbose")) {

  if(!(class(d) %in% c("SpatialPoints","SpatialPointsDataFrame")))
    stop("Error: d must be SpatialPoints(DataFrame)")
  if(!(class(g) %in% c("SpatialPixels","SpatialPixelsDataFrame"))) 
    stop("Error: g must be SpatialPixels(DataFrame)")
  if(!(class(p) %in% c("SpatialPoints","SpatialPointsDataFrame")))
    stop("Error: p must be SpatialPoints(DataFrame)")
  if(!(criterion %in% c("logdet","VV","AV","EAC"))) 
    stop("Error: criterion must be one of logdet, VV, AV or EAC")
  
  # set initial temperature
  T <- initialTemperature
  
  # merge supplemental sample and spatial coverage sample (if present)
  dall <- d
  if(!missing(spcosam)){
    dall <- rbind(d,spcosam)
  }
  
  # compute the criterion (mean kriging variance)
  E <- getCriterion.EK(dall,p,pars,perturbation,criterion)
  
  # store criterion
  E_prv <- E
  
  # Define structure for storing time series of criterion
  Eall<-NULL

  # initialize number of zero changes of objective function
  nNoChange <-0
  
  # start cooling loop
  repeat{
    
    # initialize number of accepted configurations
    nAccepted <- 0
    
    # initialize number of permuted configurations
    nPermuted <- 0
    
    # initialize number of improved configurations
    nImproved <- 0
    
    # start permutation loop
    repeat {
      
      # increase the number of permutations
      nPermuted <- nPermuted + 1
      
      # propose new sample by making use of function permute
      d_p <- permute(d, g)

      # for KED overlay new sample with grid
      if(length(names(p))>0) {
      d_p <- SpatialPointsDataFrame(
        coords = d_p,
        data = d_p %over% g                
      )}

      #merge supplemental sample and spatial coverage sample
      dall_p <- d_p
      if(!missing(spcosam)){
        dall_p <- rbind(d_p,spcosam)
      }
      
      # compute the criterion of this new sample by using function getCriterion
      E_p <- getCriterion.EK(dall_p,p,pars,perturbation,criterion)

      # accept/reject proposal by means of Metropolis criterion
      dE <- E_p - E
    
      if (dE < 0) {
        nImproved <- nImproved + 1
        prob <- 1 # always accept improvements
      } else {
        prob <- exp(-dE / T) # use Boltzmann to judge if deteriorations should be accepted
      }
      u <- runif(n = 1) # draw uniform deviate
      if (u < prob) { # accept proposal
        nAccepted <- nAccepted + 1
        d <- d_p
        E <- E_p
      }
      
      # are conditions met to lower temperature?
      lowerTemperature <- (nPermuted == maxPermuted) |
        (nAccepted == maxAccepted)
      if (lowerTemperature) {
        if (nImproved==0)
        {nNoChange<-nNoChange+1}
        else
        {nNoChange<-0}
        Eall<-rbind(Eall,E)
        break  
      }
    }
    
    if (verbose) {
      cat(
        format(Sys.time()), "|",
        sprintf("T = %e  E = %e  permuted = %d  accepted = %d  improved = %d  acceptance rate = %f \n",
                T, E, nPermuted, nAccepted, nImproved, nAccepted / nPermuted)
      )
    }

    # check on convergence
    if (nNoChange == maxNoChange) {
      break
    }
    
    E_prv <- E
    
    # lower temperature
    T <- coolingRate * T
  }
  
  # return result
  list(
    optSample=d,Criterion=Eall
  )
}

getCriterion.EK<-function(d,p,pars,perturbation,criterion)  {
  nobs <- length(d)
  D <- spDists(d)
  
  C <- CorF(D,pars)
  pars.pert <- pars
  pC <- dC <- list()
  for (i in 1:length(pars)) {
    pars.pert[i] <- (1+perturbation)*pars[i]
    pC[[i]] <- CorF(D,pars.pert)
    dC[[i]] <- (pC[[i]]-C)/(pars[i]*perturbation)
    pars.pert <- pars
  }
  
  cholC <- try(chol(C),silent=TRUE)
  if (is.character(cholC)){
    return(1E20)} else {
      # inverse of the covariance matrix
      invC <- chol2inv(chol(C))
      # compute Fisher information matrix, see Eq. 7 Geoderma paper Lark, 2002
      F <- matrix(0,length(pars),length(pars))
      for (i in 1:length(pars)){
        for (j in i:length(pars)){
          F[i,j]=F[j,i]=0.5*matrix.trace(invC%*%dC[[i]]%*%invC%*%dC[[j]])
        }
      }
      
      cholF <- try(chol(F),silent=TRUE)
      if (is.character(cholF)){
        return(1E20)} else {
          
        # inverse of the Fisher information matrix
        invF <- chol2inv(chol(F))
        
        if(criterion=="logdet"){
          logdet <- -1*determinant(F,logarithm=TRUE)$modulus #This is equal to determinant(invF,logarithm=T)$modulus
          return(logdet)} else {          
        # add dummy variable
        if(class(d)=="SpatialPoints") {
          d <- SpatialPointsDataFrame(
            coords = d,
            data = data.frame(dum = rep(1, times = length(d)))
          )
        } else {
          d$dum=1
        }
        
        if(length(names(p))>0) {
          formul <- as.formula(paste("dum", paste(names(p), collapse = "+"), sep = "~"))} else {
          formul <- as.formula(paste("dum", paste(1, collapse = "+"), sep = "~"))
        }
        
        m = model.frame(terms(formul), as(d, "data.frame"), na.action = na.fail)
        term = attr(m, "terms")
        X = model.matrix(term, m)
        
        terms.f = delete.response(terms(formul))
        mf.f = model.frame(terms.f, as(p,"data.frame")) #, na.action = na.action)
        x0 = model.matrix(terms.f, mf.f)
        
        nrowA <- nobs + ncol(X)
        A <- matrix(nrow=nrowA,ncol=nrowA)
        A[,] <- 0                      
        A[1:nobs,1:nobs] <- C
        A[1:nobs,(nobs+1):nrowA] <- X
        A[(nobs+1):nrowA,1:nobs] <- t(X)

        #compute matrix with covariances between prediction nodes and sampling points
        D0 <- spDists(x=p,y=d)
        C0 <- CorF(D0,pars)

        pars.pert <- pars
        pA <- pC0 <- pb <-list()
        for (i in 1:length(pars)) {
          pA[[i]] <- A
          pA[[i]][1:nobs,1:nobs] <- pC[[i]]
          
          pars.pert[i] <- (1+perturbation)*pars[i]
          pC0[[i]] <- CorF(D0,pars.pert)
          pb[[i]] <- cbind(pC0[[i]],x0)
          pars.pert <- pars
        }
        
        L <- matrix(nrow=length(p),ncol=nobs) #matrix with kriging weights
        pL <- array(dim=c(length(p),length(d),length(pars))) #array with perturbed kriging weights
        var <- numeric(length=length(p)) #kriging variance
        pvar <- matrix(nrow=length(p),ncol=length(pars)) #matrix with perturebed kriging variances
        for (i in 1:length(p)) {
          b <- c(C0[i,],x0[i,])
          l <- solve(A,b)
          L[i,] <- l[1:nobs]
          var[i] <- 1 - l[1:nobs] %*% C0[i,] - crossprod(l[1:nobs],X) %*% l[-(1:nobs)]
          for (j in 1:length(pars)){
            l <- solve(pA[[j]],pb[[j]][i,])
            pL[i,,j] <- l[1:nobs]
            pvar[i,j] <- 1 - l[1:nobs] %*% pC0[[j]][i,] - crossprod(l[1:nobs],X) %*% l[-(1:nobs)]
          }
        }
                
        dvar <- dL <- list()
        for (j in 1:length(pars)) {
          dvar[[j]] <- (pvar[,j]-var)/(pars[j]*perturbation)
          dL[[j]] <- (pL[,,j] - L)/(pars[j]*perturbation)
        }         

        # compute standard deviations and correlation of variogram parameters from invF
        sigma <- matrix(nrow=nrow(invF),ncol=ncol(invF))
        diag(sigma) <- sqrt(diag(invF))
        for (i in 1:(length(pars)-1)){
          for (j in (i+1):length(pars)){
            sigma[i,j]<-sigma[j,i]<-invF[i,j]/(sigma[i,i]*sigma[j,j])
          }
        }

        #tausq: expectation of additional variance due to uncertainty in ML estimates of variogram parameters, see Eq. 5 Lark and Marchant 2018
        tausq <- numeric(length=length(p))
        tausqk <- 0
        for (k in 1:length(p)) {
          for (i in 1:length(dL)){
            for (j in 1:length(dL)){
              tausqijk <- sigma[i,j]*t(dL[[i]][k,])%*%C%*%dL[[j]][k,]
              tausqk <- tausqk+tausqijk
            }
          }
          tausq[k] <- tausqk
          tausqk<-0
        }
        varplus <- var+tausq
        MVar <- mean(varplus)}
        if (criterion=="AV"){
          return(MVar)
        }  else {
        #VV: variance of kriging variance, see Eq. 9 Lark (2002) Geoderma. This variance is computed per evaluation point
        VV <- numeric(length=length(var))
        for (i in 1:length(dvar)){
          for (j in 1:length(dvar)){
            VVij <- sigma[i,j]*dvar[[i]]*dvar[[j]]
            VV <- VV+VVij
          }
        }
        
        MVarVar <- mean(VV)}
        if (criterion=="VV"){
          return(MVarVar)
        } else {
        
        EAC <- mean(varplus+1/(2*varplus)*VV) #Estimation Adjusted Criterion of Zhu and Stein (2006), see Eq. 2.16
        return(EAC)}
      }
    }
}


# Annealing function for cLHS

anneal.cLHS<-function(d, g, legacy, breaks, pp, wO1, R,
                      initialTemperature = 1, coolingRate = 0.9, maxAccepted = 10 * nrow(coordinates(d)),
                      maxPermuted=10* nrow(coordinates(d)),maxNoChange=nrow(coordinates(d)),verbose = getOption("verbose")) {
  
  # set initial temperature
  T <- initialTemperature
  
  # merge infill sample and legacy sample
  dall <- d
  if(!missing(legacy)) {
    
    #    if(class(legacy) != "SpatialPointsDataFrame") {
    #      stop("legacy should be SpatialPointsDataFrame")
    #    }
    
    #    if(proj4string(d) != proj4string(legacy)) {
    #      stop("projections don't match")
    #    }
    dall <- rbind(d,legacy)
  }
  
  # compute the criterion
  criterion <- getCriterion.cLHS(dall, g, breaks, pp, wO1, R)
  
  # store criterion
  criterion_prv <- criterion
  
  # Define structure for storing time series of criterion
  Eall<-NULL
  
  # initialize number of zero changes of objective function
  nNoChange <-0
  
  # start cooling loop
  repeat{
    
    # initialize number of accepted configurations
    nAccepted <- 0
    
    # initialize number of permuted configurations
    nPermuted <- 0
    
    # initialize number of improved configurations
    nImproved <- 0
    
    # start permutation loop
    repeat {
      
      # increase the number of permutations
      nPermuted <- nPermuted + 1
      
      # propose new sample by making use of function permute
      d_p <- permute(d, g)
      
      #merge infill sample and legacy sample
      dall_p <- d_p
      if(!missing(legacy)){
        dall_p <- rbind(d_p,legacy)
      }
      
      # compute the criterion of this new sample by using function getCriterion
      criterion_p <- getCriterion.cLHS(dall_p, g, breaks, pp, wO1, R)
      
      # accept/reject proposal by means of Metropolis criterion
      dE <- criterion_p["E"] - criterion["E"]
      if (dE < 0) {
        nImproved <- nImproved + 1
        p <- 1 # always accept improvements
      } else {
        p <- exp(-dE / T) # use Boltzmann to judge if deteriorations should be accepted
      }
      u <- runif(n = 1) # draw uniform deviate
      if (u < p) { # accept proposal
        nAccepted <- nAccepted + 1
        d <- d_p
        criterion <- criterion_p
      }
      
      # are conditions met to lower temperature?
      lowerTemperature <- (nPermuted == maxPermuted) |
        (nAccepted == maxAccepted)
      if (lowerTemperature) {
        if (nImproved==0)
        {nNoChange<-nNoChange+1}
        else
        {nNoChange<-0}
        Eall<-rbind(Eall,criterion)
        break  
      }
    }
    
    if (verbose) {
      cat(
        format(Sys.time()), "|",
        sprintf("T = %e  E = %e  permuted = %d  accepted = %d  improved = %d  acceptance rate = %f  \n",
                T, criterion["E"], nPermuted, nAccepted, nImproved, nAccepted / nPermuted)
      )
    }
    
    # check on convergence
    if (nNoChange == maxNoChange) {
      break
    }
    criterion_prv <- criterion
    
    # lower temperature
    T <- coolingRate * T
  }
  
  # return result
  list(
    optSample=d,Criterion=Eall
  )
}
# Function for computing minimization criterion of cLHS

getCriterion.cLHS<-function(d,g,breaks,pp,wO1,R)  {
  #determine values of covariates at locations in d
  d <- SpatialPointsDataFrame(
    coords = d,
    data = over(d,g)                
  )

  counts <- lapply(1:ncol(d), function (i) 
    hist(as.data.frame(d[,i])[,1], breaks[,i], plot = FALSE)$counts
  )
  countslf <- data.frame(counts=unlist(counts))
  sampleprop <- countslf/length(d)
  
  
  O1<-sum(abs(sampleprop-popprop))

  #compute sum of absolute differences of correlations
  r<-cor(as.data.frame(d)[1:ncol(d)])
  dr <- abs(R-r)
  offdiagonal <- (!row(dr)==col(dr))
  O3<-mean(dr[offdiagonal])
  
  #compute LHS criterion
  E<-wO1*O1+(1-wO1)*O3
  
  # return result
  c(E = E, O1 = O1, O3=O3)
}

