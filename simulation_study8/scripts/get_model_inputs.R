################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, assumeType, 
                          smoothC, smoothD, hospData, deathData) {
    
    # constants that are the same for all models
    N <- 1e6
    I0 <- 5
    R0 <- 0
    S0 <- N - I0 - R0
    tau <- length(incData)
    maxInf <- 10
    
    ### initial conditions probability
    if (modelType %in% c('SIHRD_full', 'SIHRD_noAlarm')) {
        
        # SIHRD model - start of epidemic so no H, R, D
        initProb <- c(S0, I0, 0, 0, 0)/N
        comp_init <- rmulti(1, N, initProb)
        
    } else if (modelType %in% c('SIR_full', 'SIR_inc', 'SIR_noAlarm')) {
        
        # SIR models - start of epidemic so R[0] = 0
        initProb <- c(S0, rep(I0/maxInf, maxInf), R0)/N
        comp_init <- rmulti(1, N, initProb)
        
    }
    
    ### constants
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          maxInf = maxInf)
    
    # conditional inference - fix unobserved Istar at 4 times what was observed
    # (for models that allow undetected infections)
    Istar <- round(incData / 0.25)
    
    ### data
    dataList <- list(detectIstar = incData,
                     Istar = Istar,
                     smoothC = smoothC,
                     smoothD = smoothD)

    if (modelType == 'SIHRD_full') {
        
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
                              k = runif(1, 0, 0.02),
                              alpha = rbeta(1, 1, 1),
                              RstarI = round(0.3 * c(rep(0, 3), I0, dataList$detectIstar[1:(tau-4)])),
                              RstarH = round(0.3 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)

            if ((probIH + probIR < 1) & (probHR + probHD < 1)) break
            
        } 
        
        names(initsList$RstarI) <- paste0('RstarI[', 1:tau, ']')
        names(initsList$RstarH) <- paste0('RstarH[', 1:tau, ']')
        
        # end modeltype == 'SIHRD_full'
        
    } else if (modelType == 'SIR_full') {
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          probDetect = rbeta(1, 250, 750),
                          beta = runif(1, 1/7, 1),
                          k = runif(1, 0, 0.02),
                          alpha = rbeta(1, 1, 1),
                          w0 = rnorm(1, 5, 0.25),
                          nu = rgamma(1, 100, 100))
        
        # end modeltype == 'SIR_full'
        
    } else if (modelType == 'SIR_inc') {

        ### data (remove smoothD)
        dataList$smoothD <- NULL
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          probDetect = rbeta(1, 250, 750),
                          beta = runif(1, 1/7, 1),
                          k = runif(1, 0, 0.02),
                          alpha = rbeta(1, 1, 1),
                          w0 = rnorm(1, 5, 0.25),
                          nu = rgamma(1, 100, 100))
        
        # end modeltype == 'SIR_inc'
        
    } else if (modelType == 'SIHRD_noAlarm') {
        
        ### data
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        dataList$smoothC <- NULL
        dataList$smoothD <- NULL
        
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
                          probDetect = rbeta(1, 250, 750),
                          beta = runif(1, 1/7, 1),
                          w0 = rnorm(1, 5, 0.25),
                          nu = rgamma(1, 100, 100))
    
        
        # end modeltype == 'SIR_noAlarm'
        
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
    nburn <- 1e5
    nthin <- 10
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin)
    
}
