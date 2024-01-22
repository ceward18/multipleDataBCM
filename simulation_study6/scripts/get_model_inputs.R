################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, assumeType, prior,
                          smoothC, smoothD, hospData, deathData) {
    
    # constants that are the same for all models
    N <- 1e6
    I0 <- 5
    R0 <- 0
    S0 <- N - I0 - R0
    tau <- length(incData)
    maxInf <- 10
    
    ### initial conditions probability
    if (modelType %in% c('full', 'fullNoAlarm')) {
        
        # SIHRD model - start of epidemic so no H, R, D
        initProb <- c(S0, I0, 0, 0, 0)/N
        comp_init <- rmulti(1, N, initProb)
        
    } else if (modelType %in% c('simple', 'inc', 'simpleNoAlarm')) {
        
        # SIR models - start of epidemic so R[0] = 0
        initProb <- c(S0, rep(I0/maxInf, maxInf), R0)/N
        comp_init <- rmulti(1, N, initProb)
        
    }
    
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
    
    # different priors on mean of Gaussian copula
    if (prior == 1) {
        
        # equal importance on deaths and cases
        zMean <- rep(0, 2)
        
    } else if (prior == 2) {
        
        # cases more important than deaths
        zMean <- c(1, -1) * 0.8
        
    } else if (prior == 3) {
        
        # deaths more important than cases
        zMean <- c(-1, 1) * 0.8
        
    }
    
    
    
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          minC = minC,
                          minD = minD,
                          maxC = maxC,
                          maxD = maxD,
                          n = n,
                          xC = xC,
                          xD = xD,
                          maxInf = maxInf,
                          zMean = zMean,
                          Sigma = Sigma)
    
    ### data
    dataList <- list(detectIstar = incData,
                     smoothC = smoothC,
                     smoothD = smoothD,
                     constrain_deltas = 1)
    
    
    if (modelType == 'simple') {
        
        repeat {
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1, 250, 750),
                              beta = runif(1, 1/7, 1),
                              nuC = rinvgamma(1, 11, 40),
                              nuD = rinvgamma(1, 11, 40),
                              x0C = runif(1, minC + 1, maxC/4),
                              x0D = runif(1, minD + 1, maxD/4),
                              Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                              w0 = rnorm(1, 5, 0.5),
                              k = rgamma(1, 100, 100))
            
            initsList$Istar <- round(incData / 0.25)
            
            delta <- multiBeta(initsList$Z)
            
            if (sum(delta) <= 1) break
        }
        
        # end modeltype == 'simple'
        
    } else if (modelType == 'full') {
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1, 250, 750),
                              beta = runif(1, 1/7, 1),
                              gamma1 = rgamma(1, 20, 100), # IR
                              gamma2 = rgamma(1, 20, 100), # HR
                              lambda = rgamma(1, 3, 100), # IH
                              phi = rgamma(1, 10, 100) ,   # HD
                              nuC = rinvgamma(1, 11, 40),
                              nuD = rinvgamma(1, 11, 40),
                              x0C = runif(1, minC + 1, maxC/4),
                              x0D = runif(1, minD + 1, maxD/4),
                              Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                              RstarI = round(0.3 * c(rep(0, 3), I0, dataList$detectIstar[1:(tau-4)])),
                              RstarH = round(0.3 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
            
            initsList$Istar <- round(incData / 0.25)
            
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
        
    } else if (modelType == 'inc') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              initProb = initProb,
                              minC = minC,
                              maxC = maxC,
                              n = n,
                              xC = xC,
                              maxInf = maxInf)
        
        ### data
        dataList <- list(detectIstar = incData,
                         smoothC = smoothC)
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          probDetect = rbeta(1, 250, 750),
                          beta = runif(1, 1/7, 1),
                          nuC = runif(1, 1, 10),
                          x0C = runif(1, minC + 1, maxC/4),
                          deltaC = runif(1, 0, 1),
                          w0 = rnorm(1, 3, 0.5),
                          k = rgamma(1, 100, 100))
        
        initsList$Istar <- round(incData / 0.25)
        
        
        # end modeltype == 'inc'
        
    } else if (modelType == 'fullNoAlarm') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              initProb = initProb,
                              maxInf = maxInf)
        
        ### data
        dataList <- list(detectIstar = incData,
                         Hstar = hospData,
                         Dstar = deathData)
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1, 250, 750),
                              beta = runif(1, 1/7, 1),
                              gamma1 = rgamma(1, 20, 100), # IR
                              gamma2 = rgamma(1, 20, 100), # HR
                              lambda = rgamma(1, 3, 100), # IH
                              phi = rgamma(1, 10, 100) ,   # HD
                              RstarI = round(0.3 * c(rep(0, 3), I0, dataList$detectIstar[1:(tau-4)])),
                              RstarH = round(0.3 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
            
            initsList$Istar <- round(incData / 0.25)
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1)) break
        }
        
        # end modeltype == 'fullNoAlarm'
        
    } else if (modelType == 'simpleNoAlarm') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              initProb = initProb,
                              maxInf = maxInf)
        
        ### data
        dataList <- list(detectIstar = incData)
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          probDetect = rbeta(1, 250, 750),
                          beta = runif(1, 1/7, 1),
                          w0 = rnorm(1, 5, 0.5),
                          k = rgamma(1, 100, 100))
        
        initsList$Istar <- round(incData / 0.25)
        
        # end modeltype == 'simpleNoAlarm'
        
    }
    
    if (assumeType == 'casesOnly') {
        initsList$probDetect <- NULL
        initsList$Istar <- NULL
        
        # change data
        dataList$detectIstar <- NULL
        dataList$Istar <- incData

    }
    
    ### MCMC specifications
    niter <- 3e5
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
