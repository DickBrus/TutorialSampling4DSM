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

anneal.K<-function(d, g, p, legacy=NULL, model, nmax = 50, prob=0.50,
                 initialTemperature = 1, coolingRate = 0.9, maxAccepted = 10 * nrow(coordinates(d)),
                 maxPermuted=10*nrow(coordinates(d)), maxNoChange=nrow(coordinates(d)), verbose = getOption("verbose")) {
  if(!(class(d) %in% c("SpatialPoints","SpatialPointsDataFrame")))
    stop("Error: d must be SpatialPoints(DataFrame)")
  if(!(class(g) %in% c("SpatialPixels","SpatialPixelsDataFrame"))) 
    stop("Error: g must be SpatialPixels(DataFrame)")
  if(!(class(p) %in% c("SpatialPoints","SpatialPointsDataFrame")))
    stop("Error: p must be SpatialPoints(DataFrame)")
  if(!is.null(legacy)){
  stopifnot(is.na(proj4string(legacy)))}
  if(prob <0 | prob > 1){
    stop("Error: prob must be in open interval (0,1)")
  }

  # set initial temperature
  T <- initialTemperature
  
  # merge infill sample and legacy sample
  dall <- d
  if(!is.null(legacy)){
    dall <- rbind(d,legacy)
  }
  
  # compute the criterion (mean kriging variance)
  E <- getCriterion.K(dall, p, model, nmax, prob)
  
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
      E_p <- getCriterion.K(dall_p, p, model, nmax, prob)
      
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


getCriterion.K<-function(d,p,model,nmax,prob) {
  
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
  quantile(result$var1.var,probs=prob)
}


# Annealing function for estimation of variogram and kriging
anneal.EK<-function(free, disc, fixed, esample, model, thetas, perturbation=0.01, criterion,
                    initialTemperature = 1, coolingRate = 0.9, maxAccepted = 10 * nrow(coordinates(free)),
                    maxPermuted=10* nrow(coordinates(free)), maxNoChange=nrow(coordinates(free)), verbose = getOption("verbose")) {
  
  if(!(class(free) %in% c("SpatialPoints","SpatialPointsDataFrame")))
    stop("Error: free must be SpatialPoints(DataFrame)")
  if(!(class(disc) %in% c("SpatialPixels","SpatialPixelsDataFrame"))) 
    stop("Error: disc must be SpatialPixels(DataFrame)")
  if(!(class(esample) %in% c("SpatialPoints","SpatialPointsDataFrame")))
    stop("Error: esample must be SpatialPoints(DataFrame)")
  if(!(criterion %in% c("logdet","VV","AV","EAC"))) 
    stop("Error: criterion must be one of logdet, VV, AV or EAC")
  
  # set initial temperature
  T <- initialTemperature
  
  # merge free and fixed sample, if present and only for criterion AV and EAC
  sample <- free
  if(!missing(fixed) & criterion %in% c("AV","EAC")){
    sample <- rbind(free,fixed)
  }
  
  # compute the criterion (mean kriging variance)
  if (criterion %in% c("logdet","VV")) {
    E <- getCriterion.E(sample=sample,grid=fixed,esample=esample,model,thetas,perturbation,criterion)} else {
    E <- getCriterion.EK(sample=sample,esample=esample,model,thetas,perturbation,criterion)
  }
  
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
      free_p <- permute(free, disc)
      
      # for KED overlay new sample with grid
      if(length(names(esample))>0) {
        free_p <- SpatialPointsDataFrame(
          coords = free_p,
          data = free_p %over% disc                
        )}
      
      #merge proposed free sample and fixed sample when present
      sample_p <- free_p
      if(!missing(fixed) & criterion %in% c("AV","EAC")){
        sample_p <- rbind(free_p,fixed)
      }
      
      # compute the criterion of this new sample by using function getCriterion
      if (criterion %in% c("logdet","VV")) {
        E_p <- getCriterion.E(sample=sample_p,grid=fixed,esample=esample,model,thetas,perturbation,criterion)} else {
        E_p <- getCriterion.EK(sample=sample_p,esample=esample,model,thetas,perturbation,criterion)
        }

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
        free <- free_p
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
    optSample=free,Criterion=Eall
  )
}

getCriterion.E<-function(sample,grid,esample,model,thetas,perturbation,criterion)  {
  nobs <- length(sample)
  #compute distance matrix of sample for variogram estimation
  D <- spDists(sample)
  A <- variogramLine(vgm(model=model,psill=thetas[1],range=thetas[2],nugget=1-thetas[1]),
                     dist_vector=D,covariance=TRUE)
  thetas.pert <- thetas
  pA <- dA <- list()
  for (i in 1:length(thetas)) {
    thetas.pert[i] <- (1+perturbation)*thetas[i]
    pA[[i]] <- variogramLine(vgm(model=model,psill=thetas.pert[1],range=thetas.pert[2],nugget=1-thetas.pert[1]),
                             dist_vector=D,covariance=TRUE)
    dA[[i]] <- (pA[[i]]-A)/(thetas[i]*perturbation)
    thetas.pert <- thetas
  }
  
  cholA <- try(chol(A),silent=TRUE)
  if (is.character(cholA)){
    return(1E20)} else {
      # inverse of the correlation matrix
      invA <- chol2inv(chol(A))
      # compute Fisher information matrix, see Eq. 7 Geoderma paper Lark, 2002
      I <- matrix(0,length(thetas),length(thetas))
      for (i in 1:length(thetas)){
        for (j in i:length(thetas)){
          I[i,j]=I[j,i]=0.5*matrix.trace(invA%*%dA[[i]]%*%invA%*%dA[[j]])
        }
      }
      
      cholI <- try(chol(I),silent=TRUE)
      if (is.character(cholI)){
        return(1E20)} else {
          
          # inverse of the Fisher information matrix
          invI <- chol2inv(chol(I))
          
          if(criterion=="logdet"){
            logdet <- determinant(invI,logarithm=TRUE)$modulus
            return(logdet)} else {          
              
              #compute distance matrix and correlation matrix of grid nodes
              D <- spDists(grid)
              A <- variogramLine(vgm(model=model,psill=thetas[1],range=thetas[2],nugget=1-thetas[1]),
                                 dist_vector=D,covariance=TRUE)
              #extend correlation matrix A with a column and row with ones (ordinary kriging)
              nobs<-length(grid)
              B <- matrix(data=0,nrow=nobs+1,ncol=nobs+1)
              B[1:nobs,1:nobs] <- A
              B[1:nobs,nobs+1] <- 1
              B[nobs+1,1:nobs] <- 1
              #compute matrix with correlations between evaluation node and sampling points
              D0 <- spDists(x=esample,y=grid)
              A0 <- variogramLine(vgm(model=model,psill=thetas[1],range=thetas[2],nugget=1-thetas[1]),
                                  dist_vector=D0,covariance=TRUE) 
              b <- cbind(A0,1)
              #compute perturbed correlation matrix (pA)
              thetas.pert <- thetas
              pA <- pA0 <- list()
              for (i in 1:length(thetas)) {
                thetas.pert[i] <- (1+perturbation)*thetas[i]
                pA[[i]] <- variogramLine(vgm(model=model,psill=thetas.pert[1],range=thetas.pert[2],nugget=1-thetas.pert[1]),
                                         dist_vector=D,covariance=TRUE) 
                pA0[[i]] <- variogramLine(vgm(model=model,psill=thetas.pert[1],range=thetas.pert[2],nugget=1-thetas.pert[1]),
                                          dist_vector=D0,covariance=TRUE)
                thetas.pert <- thetas
              }
              #extend pA and pA0 with ones
              pB <- pb <-list()
              for (i in 1:length(thetas)) {
                pB[[i]] <- matrix(data=0,nrow=nobs+1,ncol=nobs+1)
                pB[[i]][1:nobs,1:nobs] <-pA[[i]]
                pB[[i]][1:nobs,nobs+1] <- 1
                pB[[i]][nobs+1,1:nobs] <- 1
                pb[[i]] <- cbind(pA0[[i]],1)
              }
              
              #compute perturbed kriging variances (pvar)
              var <- numeric(length=length(esample)) #kriging variance
              pvar <- matrix(nrow=length(esample),ncol=length(thetas)) #matrix with perturbed kriging variances
              for (i in 1:length(esample)) {
                b <- c(A0[i,],1)
                l <- solve(B,b)
                var[i] <- 1 - l[1:nobs] %*% A0[i,] - l[nobs+1]
                for (j in 1:length(thetas)){
                  pl <- solve(pB[[j]],pb[[j]][i,])
                  pvar[i,j] <- 1 - pl[1:nobs] %*% pA0[[j]][i,] - pl[nobs+1]
                }
              }
              
              #approximate partial derivatives of kriging variance to correlogram parameters
              dvar <- list()
              for (i in 1:length(thetas)) {
                dvar[[i]] <- (pvar[,i]-var)/(thetas[i]*perturbation)
              }         
              #compute variance of kriging variance for evaluation points.
              VV <- numeric(length=length(var))
              for (i in 1:length(thetas)){
                for (j in 1:length(thetas)){
                  VVij <- invI[i,j]*dvar[[i]]*dvar[[j]]
                  VV <- VV+VVij
                }
              }
              MVV <- mean(VV)
              return(MVV)
            }
        }
    }
}


getCriterion.EK<-function(sample,esample,model,thetas,perturbation,criterion)  {
  nobs <- length(sample)
  D <- spDists(sample)
  A <- variogramLine(vgm(model=model,psill=thetas[1],range=thetas[2],nugget=1-thetas[1]),
                     dist_vector=D,covariance=TRUE)
  thetas.pert <- thetas
  pA <- dA <- list()
  for (i in 1:length(thetas)) {
    thetas.pert[i] <- (1+perturbation)*thetas[i]
    pA[[i]] <- variogramLine(vgm(model=model,psill=thetas.pert[1],range=thetas.pert[2],nugget=1-thetas.pert[1]),
                             dist_vector=D,covariance=TRUE)
    dA[[i]] <- (pA[[i]]-A)/(thetas[i]*perturbation)
    thetas.pert <- thetas
  }
  
  cholA <- try(chol(A),silent=TRUE)
  if (is.character(cholA)){
    return(1E20)} else {
      # inverse of the covariance matrix
      invA <- chol2inv(chol(A))
      # compute Fisher information matrix, see Eq. 7 Geoderma paper Lark, 2002
      I <- matrix(0,length(thetas),length(thetas))
      for (i in 1:length(thetas)){
        for (j in i:length(thetas)){
          I[i,j]=I[j,i]=0.5*matrix.trace(invA%*%dA[[i]]%*%invA%*%dA[[j]])
        }
      }
      
      cholI <- try(chol(I),silent=TRUE)
      if (is.character(cholI)){
        return(1E20)} else {
          
        # inverse of the Fisher information matrix
        invI <- chol2inv(chol(I))
        
        # add dummy variable
        if(class(sample)=="SpatialPoints") {
          sample <- SpatialPointsDataFrame(
            coords = sample,
            data = data.frame(dum = rep(1, times = length(sample)))
          )
        } else {
          sample$dum=1
        }
        
        if(length(names(esample))>0) {
          formul <- as.formula(paste("dum", paste(names(esample), collapse = "+"), sep = "~"))} else {
          formul <- as.formula(paste("dum", paste(1, collapse = "+"), sep = "~"))
        }
        
        m = model.frame(terms(formul), as(sample, "data.frame"), na.action = na.fail)
        term = attr(m, "terms")
        X = model.matrix(term, m)
        
        terms.f = delete.response(terms(formul))
        mf.f = model.frame(terms.f, as(esample,"data.frame"))
        x0 = model.matrix(terms.f, mf.f)
        
        nrowB <- nobs + ncol(X)
        B <- matrix(data=0,nrow=nrowB,ncol=nrowB)
        B[1:nobs,1:nobs] <- A
        B[1:nobs,(nobs+1):nrowB] <- X
        B[(nobs+1):nrowB,1:nobs] <- t(X)

        #compute matrix with covariances between prediction nodes and sampling points
        D0 <- spDists(x=esample,y=sample)
        A0 <- variogramLine(vgm(model=model,psill=thetas[1],range=thetas[2],nugget=1-thetas[1]),
                            dist_vector=D0,covariance=TRUE) 
        #compute pB and pb by extending pA and pA0 with X
        thetas.pert <- thetas
        pB <- pA0 <- pb <-list()
        for (i in 1:length(thetas)) {
          pB[[i]] <- B
          pB[[i]][1:nobs,1:nobs] <- pA[[i]]
          
          thetas.pert[i] <- (1+perturbation)*thetas[i]
          pA0[[i]] <- variogramLine(vgm(model=model,psill=thetas.pert[1],range=thetas.pert[2],nugget=1-thetas.pert[1]),
                                    dist_vector=D0,covariance=TRUE)
          pb[[i]] <- cbind(pA0[[i]],x0)
          thetas.pert <- thetas
        }
        
        L <- matrix(nrow=length(esample),ncol=nobs) #matrix with kriging weights
        pL <- array(dim=c(length(esample),length(sample),length(thetas))) #array with perturbed kriging weights
        var <- numeric(length=length(esample)) #kriging variance
        pvar <- matrix(nrow=length(esample),ncol=length(thetas)) #matrix with perturbed kriging variances
        for (i in 1:length(esample)) {
          b <- c(A0[i,],x0[i,])
          l <- solve(B,b)
          L[i,] <- l[1:nobs]
          var[i] <- 1 - l[1:nobs] %*% A0[i,] - x0[i,] %*% l[-(1:nobs)]
          for (j in 1:length(thetas)){
            l <- solve(pB[[j]],pb[[j]][i,])
            pL[i,,j] <- l[1:nobs]
            pvar[i,j] <- 1 - l[1:nobs] %*% pA0[[j]][i,] - x0[i,] %*% l[-(1:nobs)]
          }
        }
                
        dvar <- dL <- list()
        for (i in 1:length(thetas)) {
          dvar[[i]] <- (pvar[,i]-var)/(thetas[i]*perturbation)
          dL[[i]] <- (pL[,,i] - L)/(thetas[i]*perturbation)
        }         

        #tausq: expectation of additional variance due to uncertainty in ML estimates of variogram parameters, see Eq. 5 Lark and Marchant 2018
        tausq <- numeric(length=length(esample))
        tausqk <- 0
        for (k in 1:length(esample)) {
          for (i in 1:length(dL)){
            for (j in 1:length(dL)){
              tausqijk <- invI[i,j]*t(dL[[i]][k,])%*%A%*%dL[[j]][k,]
              tausqk <- tausqk+tausqijk
            }
          }
          tausq[k] <- tausqk
          tausqk<-0
        }
        augmentedvar <- var+tausq
        MVar <- mean(augmentedvar)
        if (criterion=="AV"){
          return(MVar)
        }  else {
        #VV: variance of kriging variance, see Eq. 9 Lark (2002) Geoderma. This variance is computed per evaluation point
        VV <- numeric(length=length(var))
        for (i in 1:length(dvar)){
          for (j in 1:length(dvar)){
            VVij <- invI[i,j]*dvar[[i]]*dvar[[j]]
            VV <- VV+VVij
          }
        }
        EAC <- mean(augmentedvar+VV/(2*var)) #Estimation Adjusted Criterion of Zhu and Stein (2006), see Eq. 2.16
        return(EAC)
      }
    }
  }
}


# Annealing function for cLHS

anneal.cLHS<-function(d, g, legacy, lb, wO1, R,
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
  criterion <- getCriterion.cLHS(dall, g, lb, wO1,R)
  
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
      criterion_p <- getCriterion.cLHS(dall_p, g, lb, wO1, R)
      
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

getCriterion.cLHS<-function(d,g,lb,wO1,R)  {
  #determine values of covariates at locations in d
  d <- SpatialPointsDataFrame(
    coords = d,
    data = over(d,g)                
  )
  
  #Determine in which stratum the sampling locations are
  stratum<-matrix(nrow=length(d),ncol=ncol(d))
  for ( i in 1:ncol(d) ) {
    stratum[,i]<-findInterval(as.data.frame(d[,i])[,1],lb[,i])
  }
  
  #count number of points in marginal strata
  counts<-matrix(nrow=nrow(lb),ncol=ncol(d))
  for (i in 1:nrow(lb)) {
    counts[i,]<-apply(stratum, MARGIN=2, function(x,i) sum(x==i), i=i)
  }
  O1<-mean(abs(counts-1))
  
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

