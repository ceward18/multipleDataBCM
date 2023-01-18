################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, smoothC, smoothD,
                          hospData, deathData) {
    
    # constants that are the same for all models
    N <- 1e6
    I0 <- 5
    R0 <- 0
    S0 <- N - I0 - R0
    tau <- length(incData)
    
    ### constants
    n <- 50
    
    minC <- floor(min(smoothC))
    minD <- floor(min(smoothD))
    maxC <- ceiling(max(smoothC))
    maxD <- ceiling(max(smoothD))
    
    xC <- seq(0, maxC, length.out = n)
    xD <- seq(0, maxD, length.out = n)
    
    rho <- -0.6
    corMat <- matrix(c(1, rho, 
                       rho, 1), ncol = 2, byrow = T)
    sds <- rep(0.7, 2)
    Sigma <- diag(sds) %*% corMat %*% diag(sds)
    
    constantsList <- list(tau = tau,
                          N = N,
                          S0 = S0,
                          I0 = I0,
                          minC = minC,
                          minD = minD,
                          maxC = maxC,
                          maxD = maxD,
                          n = n,
                          xC = xC,
                          xD = xD,
                          maxInf = 10,
                          zeros = rep(0, 2),
                          Sigma = Sigma)
    
    ### data
    dataList <- list(Istar = incData,
                     smoothC = smoothC,
                     smoothD = smoothD,
                     constrain_deltas = 1)
    
    if (modelType == 'simple') {
        
        repeat {
            ### inits 
            initsList <- list(beta = runif(1, 1/7, 1),
                              nuC = runif(1, 1, 10),
                              nuD = runif(1, 1, 10),
                              x0C = runif(1, maxC/10, maxC/5),
                              x0D = runif(1, maxD/10, maxD/5),
                              Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                              w0 = rnorm(1, 3, 0.5),
                              k = rgamma(1, 100, 100))
            
            
            delta <- multiBeta(initsList$Z)
            
            if (sum(delta) <= 1) break
        }
        
        # end modeltype == 'simple'
        
    } else if (modelType == 'simpleThresh') {
      
      repeat {
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          HC = runif(1, 0, maxC/N/5),
                          HD = runif(1, 0, maxD/N/5),
                          Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                          w0 = rnorm(1, 3, 0.5),
                          k = rgamma(1, 100, 100))
        
        
        delta <- multiBeta(initsList$Z)
        
        if (sum(delta) <= 1) break
      }
      
      
      # end modeltype == 'simpleThresh'
      
    } else if (modelType == 'full') {
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(beta = runif(1, 1/7, 1),
                              gamma1 = runif(1), # IR
                              gamma2 = runif(1), # HR
                              lambda = runif(1), # IH
                              phi = runif(1) ,   # HD
                              nuC = runif(1, 1, 10),
                              nuD = runif(1, 1, 10),
                              x0C = runif(1, maxC/20, maxC/5),
                              x0D = runif(1, maxD/20, maxD/5),
                              Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                              RstarI = round(0.3 * c(rep(0, 3), I0, dataList$Istar[1:(tau-4)])),
                              RstarH = round(0.3 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            
            delta <- multiBeta(initsList$Z)
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1) & sum(delta) <= 1) break
            
        } 
        
        names(initsList$RstarI) <- paste0('RstarI[', 1:tau, ']')
        names(initsList$RstarH) <- paste0('RstarH[', 1:tau, ']')
        
        # end modeltype == 'full'
        
    }  else if (modelType == 'fullThresh') {
      
      dataList$Hstar <- hospData
      dataList$Dstar <- deathData
      
      repeat {
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          gamma1 = runif(1), # IR
                          gamma2 = runif(1), # HR
                          lambda = runif(1), # IH
                          phi = runif(1) ,   # HD
                          HC = runif(1, 0, maxC/N/5),
                          HD = runif(1, 0, maxD/N/5),
                          Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                          RstarI = round(0.3 * c(rep(0, 3), I0, dataList$Istar[1:(tau-4)])),
                          RstarH = round(0.3 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
        
        probIH <- 1 - exp(-initsList$lambda)
        probIR <- 1 - exp(-initsList$gamma1)
        
        probHR <- 1 - exp(-initsList$gamma2)
        probHD <- 1 - exp(-initsList$phi)
        
        
        delta <- multiBeta(initsList$Z)
        
        if ((probIH + probIR < 1) & (probHR + probHD < 1) & sum(delta) <= 1) break
        
      } 
      
      names(initsList$RstarI) <- paste0('RstarI[', 1:tau, ']')
      names(initsList$RstarH) <- paste0('RstarH[', 1:tau, ']')
      
      # end modeltype == 'fullThresh'
      
    } else if (modelType == 'inc') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              minC = minC,
                              maxC = maxC,
                              n = n,
                              xC = xC,
                              maxInf = 10,
                              zeros = rep(0, 2),
                              Sigma = Sigma)
        
        ### data
        dataList <- list(Istar = incData,
                         smoothC = smoothC)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          nuC = runif(1, 1, 10),
                          x0C = runif(1, maxC/20, maxC/5),
                          deltaC = runif(1, 0, 1),
                          w0 = rnorm(1, 3, 0.5),
                          k = rgamma(1, 100, 100))
        
        
        
        # end modeltype == 'inc'
        
    } 
    
    ### MCMC specifications
    niter <- 4e5
    nburn <- 2e5
    nthin <- 10
    
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xC = xC,
         xD = xD)
    
}
