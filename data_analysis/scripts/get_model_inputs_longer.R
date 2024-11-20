################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput_longer <- function(incData, city, modelType, peak, timePeriod,
                          smoothC, smoothD,
                          hospData, deathData,
                          N, S0, I0, H0, D0, R0, chainIdx) {
    
    maxInf <- 14
    tau <- length(incData)
    
    # read in previous saved chains
    lastRun <- readRDS(paste0('./output/chains_', modelType, '_', city,
                               '_peak', peak, '_weeks', timePeriod, '.rds'))
    
    chainSamples <- lastRun[[chainIdx]]
    lastChainIdx <- nrow(chainSamples)
    
    #### initial conditions probability
    if (grepl('SIHRD', modelType)) {
        
        # SIHRD model
        initProb <- c(S0, I0, H0, R0, D0)/N
        
    } else if (grepl('SIR', modelType)) {
        
        # SIR models 
        initProb <- c(S0, rep(I0/maxInf, maxInf), R0)/N
        
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
    
    # use last iteration of previous chains as initial values
    
    RstarI <- as.numeric(chainSamples[lastChainIdx,grep('RstarI', colnames(chainSamples))])
    RstarH <- as.numeric(chainSamples[lastChainIdx,grep('RstarH', colnames(chainSamples))])
    
    
    if (modelType == 'SIHRD_full') {
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = as.numeric(chainSamples[lastChainIdx,grep('comp_init', colnames(chainSamples))]),
                              probDetect = as.numeric(chainSamples[lastChainIdx,grep('probDetect', colnames(chainSamples))]),
                              beta = as.numeric(chainSamples[lastChainIdx,grep('beta', colnames(chainSamples))]),
                              gamma1 = as.numeric(chainSamples[lastChainIdx,grep('gamma1', colnames(chainSamples))]),  # IR
                              gamma2 = as.numeric(chainSamples[lastChainIdx,grep('gamma2', colnames(chainSamples))]), # HR
                              lambda = as.numeric(chainSamples[lastChainIdx,grep('lambda', colnames(chainSamples))]), # IH
                              phi = as.numeric(chainSamples[lastChainIdx,grep('phi', colnames(chainSamples))]),    # HD
                              k = as.numeric(chainSamples[lastChainIdx,grep('k', colnames(chainSamples))]),
                              alpha = as.numeric(chainSamples[lastChainIdx,grep('alpha', colnames(chainSamples))]),
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
            initsList <- list(comp_init = as.numeric(chainSamples[lastChainIdx,grep('comp_init', colnames(chainSamples))]),
                              probDetect = as.numeric(chainSamples[lastChainIdx,grep('probDetect', colnames(chainSamples))]),
                              beta = as.numeric(chainSamples[lastChainIdx,grep('beta', colnames(chainSamples))]),
                              gamma1 = as.numeric(chainSamples[lastChainIdx,grep('gamma1', colnames(chainSamples))]),  # IR
                              gamma2 = as.numeric(chainSamples[lastChainIdx,grep('gamma2', colnames(chainSamples))]), # HR
                              lambda = as.numeric(chainSamples[lastChainIdx,grep('lambda', colnames(chainSamples))]), # IH
                              phi = as.numeric(chainSamples[lastChainIdx,grep('phi', colnames(chainSamples))]),    # HD
                              k = as.numeric(chainSamples[lastChainIdx,grep('k', colnames(chainSamples))]),
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
        initsList <- list(comp_init = as.numeric(chainSamples[lastChainIdx,grep('comp_init', colnames(chainSamples))]),
                          probDetect = as.numeric(chainSamples[lastChainIdx,grep('probDetect', colnames(chainSamples))]),
                          beta = as.numeric(chainSamples[lastChainIdx,grep('beta', colnames(chainSamples))]),
                          k = as.numeric(chainSamples[lastChainIdx,grep('k', colnames(chainSamples))]),
                          alpha = as.numeric(chainSamples[lastChainIdx,grep('alpha', colnames(chainSamples))]),
                          w0 = as.numeric(chainSamples[lastChainIdx,grep('w0', colnames(chainSamples))]),
                          nu = as.numeric(chainSamples[lastChainIdx,grep('nu', colnames(chainSamples))]))
        
        # end modeltype == 'SIR_full'
        
    } else if (modelType == 'SIR_inc') {
        
        ### data (remove smoothD)
        dataList$smoothD <- NULL
        
        ### inits 
        initsList <- list(comp_init = as.numeric(chainSamples[lastChainIdx,grep('comp_init', colnames(chainSamples))]),
                          probDetect = as.numeric(chainSamples[lastChainIdx,grep('probDetect', colnames(chainSamples))]),
                          beta = as.numeric(chainSamples[lastChainIdx,grep('beta', colnames(chainSamples))]),
                          k = as.numeric(chainSamples[lastChainIdx,grep('k', colnames(chainSamples))]),
                          w0 = as.numeric(chainSamples[lastChainIdx,grep('w0', colnames(chainSamples))]),
                          nu = as.numeric(chainSamples[lastChainIdx,grep('nu', colnames(chainSamples))]))
        
        # end modeltype == 'SIR_inc'
        
    } else if (modelType == 'SIHRD_noAlarm') {
        
        ### data
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        dataList$smoothC <- NULL
        dataList$smoothD <- NULL
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = as.numeric(chainSamples[lastChainIdx,grep('comp_init', colnames(chainSamples))]),
                              probDetect = as.numeric(chainSamples[lastChainIdx,grep('probDetect', colnames(chainSamples))]),
                              beta = as.numeric(chainSamples[lastChainIdx,grep('beta', colnames(chainSamples))]),
                              gamma1 = as.numeric(chainSamples[lastChainIdx,grep('gamma1', colnames(chainSamples))]),  # IR
                              gamma2 = as.numeric(chainSamples[lastChainIdx,grep('gamma2', colnames(chainSamples))]), # HR
                              lambda = as.numeric(chainSamples[lastChainIdx,grep('lambda', colnames(chainSamples))]), # IH
                              phi = as.numeric(chainSamples[lastChainIdx,grep('phi', colnames(chainSamples))]),    # HD
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
        initsList <- list(comp_init = as.numeric(chainSamples[lastChainIdx,grep('comp_init', colnames(chainSamples))]),
                          probDetect = as.numeric(chainSamples[lastChainIdx,grep('probDetect', colnames(chainSamples))]),
                          beta = as.numeric(chainSamples[lastChainIdx,grep('beta', colnames(chainSamples))]),
                          w0 = as.numeric(chainSamples[lastChainIdx,grep('w0', colnames(chainSamples))]),
                          nu = as.numeric(chainSamples[lastChainIdx,grep('nu', colnames(chainSamples))]))
        
        # end modeltype == 'SIR_noAlarm'
        
    }
    
    ### MCMC specifications
    niter <- 8e5
    nburn <- 0
    nthin <- 10
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin)
    
}
