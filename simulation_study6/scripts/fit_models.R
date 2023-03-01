################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, modelType, alarmBase, 
                          smoothC, smoothD, hospData, deathData, seed) {
    
    source('./scripts/model_code.R')
    source('./scripts/get_model_inputs.R')
    
    # for reproducibility so inits are always the same
    set.seed(seed + 3)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, smoothC, smoothD,
                                 hospData, deathData)
    
    ### MCMC specifications
    niter <- modelInputs$niter
    nburn <- modelInputs$nburn
    nthin <- modelInputs$nthin
    
    ### get appropriate model code
    modelCode <- get(paste0('SIHRD_', modelType))
    
    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    myConfig <- configureMCMC(myModel)

    myConfig$removeSampler('probDetect')
    myConfig$addSampler(target = 'probDetect', type = "slice")
    
    if (modelType == 'full') {
        
        myConfig$addMonitors(c('yAlarmC', 'yAlarmD', 'alarmC', 'alarmD', 
                               'R0', 'delta'))
        
        paramsForBlock <- c('beta', 'Z', 'x0C', 'x0D', 'nuC', 'nuD')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "RW_block",
                            control = list(propCov = diag(c(0.3,
                                                            0.7, 0.7, 
                                                            150, 50, 
                                                            3, 3)))) 
        
        # samplers for RstarI and RstarH
        myConfig$removeSamplers('RstarI') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarI'),
                            type = "RstarUpdate")
        
        myConfig$removeSamplers('RstarH') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarH'),
                            type = "RstarUpdate")
        
        myConfig$addMonitors(c('RstarI', 'RstarH'))
        
        # use slice sampling for rate parameters
        paramsForSlice <- c('gamma1', 'gamma2', 'lambda','phi')
        myConfig$removeSampler(paramsForSlice)
        for (j in 1:length(paramsForSlice)) {
            myConfig$addSampler(target = paramsForSlice[j], type = "slice")
        }
    } else if (modelType == 'simple') {
        
        myConfig$addMonitors(c('yAlarmC', 'yAlarmD', 'alarmC', 'alarmD', 
                               'R0', 'delta'))
        
        # use slice sampling for hill parameters
        paramsForBlock <- c('Z', 'x0C', 'x0D', 'nuC', 'nuD')
        myConfig$removeSampler(paramsForBlock)
        myConfig$addSampler(target = paramsForBlock, 
                            type = "RW_block",
                            control = list(propCov = diag(c(1.7, 1.7, 
                                                            30^2, 10^2, 
                                                            3, 3))))
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('Z', 'x0C', 'x0D', 'nuC', 'nuD')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    } else { # if model == 'inc'
        
        myConfig$addMonitors(c('yAlarmC','alarmC', 'R0', 'deltaC'))
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('deltaC', 'x0C', 'nuC')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    }
    
    myConfig$removeSamplers('Istar') # Nodes will be expanded
    myConfig$addSampler(target = c('Istar'),
                        type = "RstarUpdate",
                        control = list(nUpdates = 500))
    myConfig$addMonitors(c('Istar'))
    
    print(myConfig)
    
    # browser()
    nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
    myMCMC <- buildMCMC(myConfig)
    compiled <- compileNimble(myModel, myMCMC) 
    
    runMCMC(compiled$myMCMC, 
            niter = niter, 
            nburnin = nburn,
            thin = nthin,
            setSeed = seed)
    
}


