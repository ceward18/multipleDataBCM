################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, city, modelType, peak,
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
    
    # prior for probDetect depends on peak 
    # NYC:
    #   (https://covid19.healthdata.org/united-states-of-america/new-york?view=infections-testing&tab=trend&test=infections)
    # wave 1: Feb 25 - 11 July 2020          17% detection
    # wave 2: Aug 23, 2020 - March 20, 2021  55% detection
    # wave 4: July 18 - Dec 4, 2021          20% detection
    # Montreal:
    #   (https://www.healthdata.org/sites/default/files/covid_briefs/101_briefing_Canada.pdf)
    # wave 1: Feb 25 - 11 July 2020          25% detection
    # wave 2: Aug 23, 2020 - March 20, 2021  40% detection
    # wave 4: July 18 - Dec 4, 2021          20% detection
    
    # prior for hospitalization and death rates also depends on peak
    if (city == 'nyc') {
        probDetectMean <- switch(peak_i, 
                                 '1' = 0.17,
                                 '2' = 0.55,
                                 '3' = 0.2)
    } else if (city == 'montreal') {
        probDetectMean <- switch(peak_i, 
                                 '1' = 0.25,
                                 '2' = 0.4,
                                 '3' = 0.2)
    }
    
    detectA <- probDetectMean*1000
    detectB <- 1000 - detectA
    
    ### constants
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          maxInf = maxInf,
                          detectA = detectA,
                          detectB = detectB)
    
    # conditional inference - fix unobserved Istar at 1/probDetect times what was observed
    Istar <- round(incData / probDetectMean)
    
    ### data
    dataList <- list(detectIstar = incData,
                     Istar = Istar,
                     smoothC = smoothC,
                     smoothD = smoothD)
    
    # only relevant for SIHRD models
    # distribute initial infections to be removed anytime in first maxInf days
    # ~2% infected go to hospital, so use 90% to be conservative
    RstarI <- round(0.9 * c(rmulti(1, I0, rep(I0/maxInf, maxInf)),
                            dataList$Istar[1:(tau-maxInf)]))
    
    # RstarH is removals from hospitalizations
    
    # with no removals, we would have this many in H at any time
    totalH_noR <- cumsum(c(H0, hospData)) - cumsum(c(0, deathData))
    RstarH <- round(0.25 * diff(totalH_noR))
    RstarH <- pmax(rep(0, length(RstarH)), RstarH)


    
    if (modelType == 'SIHRD_full') {
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1,  constantsList$detectA,
                                                 constantsList$detectB),
                              beta = runif(1, 0.01, 1),
                              gamma1 = rgamma(1, 1429, 10000),  # IR
                              gamma2 = rgamma(1, 6700, 100000), # HR
                              lambda = rgamma(1, 2000, 100000), # IH
                              phi = rgamma(1, 6700, 100000),    # HD
                              k = runif(1, 0.0001, 0.01),
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
        dataList$smoothD <- NULL
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1,  constantsList$detectA,
                                                 constantsList$detectB),
                              beta = runif(1, 0.01, 1),
                              gamma1 = rgamma(1, 1429, 10000),  # IR
                              gamma2 = rgamma(1, 6700, 100000), # HR
                              lambda = rgamma(1, 2000, 100000), # IH
                              phi = rgamma(1, 6700, 100000),    # HD
                              k = runif(1, 0.0001, 0.01),
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
                          beta = runif(1, 0.01, 1),
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
                          probDetect = rbeta(1, constantsList$detectA, 
                                             constantsList$detectB),
                          beta = runif(1, 0.01, 1),
                          k = runif(1, 0.0001, 0.01),
                          alpha = rbeta(1, 1, 1),
                          w0 = runif(1, 6.8, 7.2),
                          nu = rgamma(1, 1000, 1000))
        
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
                              beta = runif(1, 0.01, 1),
                              gamma1 = rgamma(1, 1429, 10000),  # IR
                              gamma2 = rgamma(1, 6700, 100000), # HR
                              lambda = rgamma(1, 2000, 100000), # IH
                              phi = rgamma(1, 6700, 100000),    # HD
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
                          beta = runif(1, 0.01, 1),
                          w0 = runif(1, 6.8, 7.2),
                          nu = rgamma(1, 1000, 1000))
        
        # end modeltype == 'SIR_noAlarm'
        
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
