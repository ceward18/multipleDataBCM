################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, peak,
                          smoothC, smoothD,
                          hospData, deathData,
                          N, S0, I0, H0, D0, R0) {
    
    maxInf <- 14
    tau <- length(incData)
    
    #### initial conditions probability
    if (grepl('SIHRD', modelType)) {
        
        # SIHRD model
        initProb <- c(S0, I0, H0, R0, D0)/N
        comp_init <- rmulti(1, N, initProb)
        
    } else if (grepl('SIR', modelType)) {
        
        # SIR models 
        initProb <- c(S0, rep(I0/maxInf, maxInf), R0)/N
        comp_init <- rmulti(1, N, initProb)
        
    }
  
    ### constants
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          maxInf = maxInf)
    
    ### data
    dataList <- list(Istar = incData,
                     smoothC = smoothC,
                     smoothD = smoothD)
    
    
    
    # only relevant for SIHRD models
    # distribute initial infections to be removed anytime in first maxInf days
    # ~2% infected go to hospital, so use 90% to be conservative
    if (peak == 1) {
        RstarI <- round(0.9 * c(rep(0, maxInf),
                          dataList$Istar[1:(tau-maxInf)]))
    } else{
        RstarI <- round(0.9 * c(rmulti(1, I0, rep(I0/maxInf, maxInf)),
                          dataList$Istar[1:(tau-maxInf)]))
    }
    
    
    # RstarH is removals from hospitalizations
    if (peak == 1) {
        
        missing_idx <- which(is.na(hospData))
        hospData[missing_idx[1:10]] <- 0
        hospData<- round(approxfun(1:length(hospData), hospData)(1:length(hospData)))
        
    } else if (peak == 2) {
        # use midpoints
        hospData<- round(approxfun(1:length(hospData), hospData)(1:length(hospData)))
        
    }
    
    
    # with no removals, we would have this many in H at any time
    totalH_noR <- cumsum(c(H0, hospData)) - cumsum(c(0, deathData))
    RstarH <- round(0.1 * diff(totalH_noR))
    RstarH <- pmax(rep(0, length(RstarH)), RstarH)


    
    if (modelType == 'SIHRD_full') {
        
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              beta = runif(1, 0.01, 0.5),
                              gamma1 = rgamma(1, 1429, 10000),  # IR
                              gamma2 = rgamma(1, 6700, 100000), # HR
                              lambda = rgamma(1, 2000, 100000), # IH
                              phi = rgamma(1, 6700, 100000),    # HD
                              k = runif(1, 0.0001, 0.01),
                              alpha = rbeta(1, 1, 1),
                              RstarI = RstarI,
                              RstarH = RstarH,
                              Hstar = hospData)
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1)) break
            
        } 
        
        names(initsList$RstarI) <- paste0('RstarI[', 1:tau, ']')
        names(initsList$RstarH) <- paste0('RstarH[', 1:tau, ']')
        
        # end modeltype == 'SIHRD_full'
        
    } else if (modelType == 'SIHRD_inc') {
        
        dataList$Dstar <- deathData
        dataList$smoothD <- NULL
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              beta = runif(1, 0.01, 0.5),
                              gamma1 = rgamma(1, 1429, 10000),  # IR
                              gamma2 = rgamma(1, 6700, 100000), # HR
                              lambda = rgamma(1, 2000, 100000), # IH
                              phi = rgamma(1, 6700, 100000),    # HD
                              k = runif(1, 0.0001, 0.01),
                              RstarI = RstarI,
                              RstarH = RstarH,
                              Hstar = hospData)
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1)) break
            
        } 
        
        names(initsList$RstarI) <- paste0('RstarI[', 1:tau, ']')
        names(initsList$RstarH) <- paste0('RstarH[', 1:tau, ']')
        
        # end modeltype == 'SIHRD_inc'
        
    } else if (modelType == 'SIR_full') {
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          beta = runif(1, 0.01, 0.5),
                          k = runif(1, 0.0001, 0.01),
                          alpha = rbeta(1, 1, 1),
                          w0 = runif(1, 6.8, 7.2),
                          nu = rgamma(1, 1000, 1000))
        
        # end modeltype == 'SIR_full'
        
    } else if (modelType == 'SIR_inc') {
        
        ### data (remove smoothD)
        dataList$smoothD <- NULL
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          beta = runif(1, 0.01, 0.5),
                          k = runif(1, 0.0001, 0.01),
                          alpha = rbeta(1, 1, 1),
                          w0 = runif(1, 6.8, 7.2),
                          nu = rgamma(1, 1000, 1000))
        
        # end modeltype == 'SIR_inc'
        
    } else if (modelType == 'SIHRD_noAlarm') {
        
        ### data
        dataList$Dstar <- deathData
        dataList$smoothC <- NULL
        dataList$smoothD <- NULL
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              beta = runif(1, 0.01, 0.5),
                              gamma1 = rgamma(1, 1429, 10000),  # IR
                              gamma2 = rgamma(1, 6700, 100000), # HR
                              lambda = rgamma(1, 2000, 100000), # IH
                              phi = rgamma(1, 6700, 100000),    # HD
                              RstarI = RstarI,
                              RstarH = RstarH,
                              Hstar = hospData)
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1)) break
        }
        
        # end modeltype == 'SIHRD_noAlarm'
        
    } else if (modelType == 'SIR_noAlarm') {
        
        ### data
        dataList$smoothC <- NULL
        dataList$smoothD <- NULL
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          beta = runif(1, 0.01, 0.5),
                          w0 = runif(1, 6.8, 7.2),
                          nu = rgamma(1, 1000, 1000))
        
        # end modeltype == 'SIR_noAlarm'
        
    }
    
    ### MCMC specifications
    niter <- 8e5
    nburn <- 3e5
    nthin <- 10
    
    ### MCMC specifications
    niter <- 100
    nburn <- 1
    nthin <- 1
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin)
    
}
