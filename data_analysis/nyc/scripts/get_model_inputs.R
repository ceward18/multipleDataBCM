################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, assumeType, peak,
                          smoothC, smoothD,
                          hospData, deathData,
                          N, S0, I0, H0, D0, R0) {
    
    maxInf <- 14
    tau <- length(incData)
    
    #### initial conditions probability
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIHRD_noAlarm')) {
        
        # SIHRD model
        initProb <- c(S0, I0, H0, R0, D0)/N
        comp_init <- rmulti(1, N, initProb)
        
        #  make sure I[1] has enough people to cover Hstar[1] and RstarI[1]
        # for undetected model only (otherwise I is known)
        if (assumeType == 'undetected') {
            
            if (I0 < H0) {
                initProb <- c(S0, I0 + H0, H0, R0, D0)/N
                comp_init <- rmulti(1, N, initProb)
            }
            
            
        }
        
    } else if (modelType %in% c('SIR_full', 'SIR_inc', 'SIR_noAlarm')) {
        
        # SIR models 
        initProb <- c(S0, rep(I0/maxInf, maxInf), R0)/N
        comp_init <- rmulti(1, N, initProb)
        
    }
    
    # prior for probDetect depends on peak 
    #   (https://covid19.healthdata.org/united-states-of-america/new-york?view=infections-testing&tab=trend&test=infections)
    # wave 1: Feb 25 - 11 July 2020          17% detection
    # wave 2: Aug 23, 2020 - March 20, 2021  55% detection
    # wave 3: March 21 - July 17, 2021       
    # wave 4: July 18 - Dec 4, 2021          20% detection
    
    # prior for hospitalization and death rates also depends on peak
    if (peak == 1) {
        probDetectMean <- 0.17
        detectA <- 170
        detectB <- 930
        
    } else if (peak == 2) {
        probDetectMean <- 0.55
        detectA <- 550
        detectB <- 450
        
    } else if (peak == 4) {
        probDetectMean <- 0.2
        detectA <- 200
        detectB <- 800
    }
    
    ### constants
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          maxInf = maxInf,
                          detectA = detectA,
                          detectB = detectB)
    
    # conditional inference - fix unobserved Istar at 4 times what was observed
    # (for models that allow undetected infections)
    Istar <- round(incData / probDetectMean)
    
    
    ### data
    dataList <- list(detectIstar = incData,
                     Istar = Istar,
                     smoothC = smoothC,
                     smoothD = smoothD)
    
    # only relevant for SIHRD models
    if (assumeType == 'casesOnly') {
        RstarI <- round(0.3 * c(rmulti(1, I0, rep(I0/7, 7)), 
                                dataList$detectIstar[1:(tau-7)]))
    } else {
        RstarI <- round(0.3 * c(rmulti(1, I0, rep(I0/7, 7)),
                                dataList$Istar[1:(tau-7)]))
    }
    
    if (H0 == 0) {
        
        RstarH <- round(0.1 * c(rep(0, 15),
                                hospData[1:(tau-15)]))
    } else {
        RstarH <- round(0.1 * c(rmulti(1, H0, rep(H0/15, 15)),
                                hospData[1:(tau-15)]))
    }
    
    
    
    if (modelType == 'SIHRD_full') {
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1,  constantsList$detectA,
                                                 constantsList$detectB),
                              beta = runif(1, 1/7, 2),
                              gamma1 = rgamma(1, 14, 100), # IR
                              gamma2 = rgamma(1, 67, 1000), # HR
                              lambda = rgamma(1, 1, 10), # IH
                              phi = rgamma(1, 67, 1000),    # HD
                              k = runif(1, 0, 0.02),
                              alpha = rbeta(1, 1, 1),
                              RstarI = RstarI,
                              RstarH = RstarH)
            
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
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1,  constantsList$detectA,
                                                 constantsList$detectB),
                              beta = runif(1, 1/7, 2),
                              gamma1 = rgamma(1, 14, 100), # IR
                              gamma2 = rgamma(1, 67, 1000), # HR
                              lambda = rgamma(1, 1, 10), # IH
                              phi = rgamma(1, 67, 1000),    # HD
                              k = runif(1, 0, 0.02),
                              RstarI = RstarI,
                              RstarH = RstarH)
            
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
                          probDetect = rbeta(1, constantsList$detectA,
                                             constantsList$detectB),
                          beta = runif(1, 1/7, 2),
                          k = runif(1, 0, 0.02),
                          alpha = rbeta(1, 1, 1),
                          w0 = rnorm(1, 7, 0.25),
                          nu = rgamma(1, 100, 100))
        
        # end modeltype == 'SIR_full'
        
    } else if (modelType == 'SIR_inc') {
        
        ### data (remove smoothD)
        dataList$smoothD <- NULL
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          probDetect = rbeta(1, constantsList$detectA, 
                                             constantsList$detectB),
                          beta = runif(1, 1/7, 2),
                          k = runif(1, 0, 0.02),
                          alpha = rbeta(1, 1, 1),
                          w0 = rnorm(1, 7, 0.25),
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
                              probDetect = rbeta(1, constantsList$detectA, constantsList$detectB),
                              beta = runif(1, 1/7, 2),
                              gamma1 = rgamma(1, 14, 100), # IR
                              gamma2 = rgamma(1, 67, 1000), # HR
                              lambda = rgamma(1, 1, 10), # IH
                              phi = rgamma(1, 67, 1000),    # HD
                              RstarI = RstarI,
                              RstarH = RstarH)
            
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
                          probDetect = rbeta(1, constantsList$detectA, 
                                             constantsList$detectB),
                          beta = runif(1, 1/7, 2),
                          w0 = rnorm(1, 7, 0.25),
                          nu = rgamma(1, 100, 100))
        
        # end modeltype == 'SIR_noAlarm'
        
    }
    
    
    if (assumeType == 'casesOnly') {
        constantsList$detectA <- NULL
        constantsList$detectB <- NULL
        initsList$probDetect <- NULL
        initsList$Istar <- NULL
        
        # change data
        dataList$detectIstar <- NULL
        dataList$Istar <- incData
        
    }
    
    ### MCMC specifications
    niter <- 8e5
    nburn <- 3e5
    nthin <- 10
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin)
    
}
