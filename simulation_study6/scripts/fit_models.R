################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, modelType, assumeType, alarmBase, 
                          smoothC, smoothD, hospData, deathData, seed) {
    
    source('./scripts/model_code.R')
    source('./scripts/get_model_inputs.R')
    
    # for reproducibility so inits are always the same
    set.seed(seed + 3)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, assumeType, smoothC, smoothD,
                                 hospData, deathData)
    
    ### MCMC specifications
    niter <- modelInputs$niter
    nburn <- modelInputs$nburn
    nthin <- modelInputs$nthin
    
    ### get appropriate model code
    modelCode <- get(paste0('SIHRD_', modelType, '_', assumeType))
    
    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    myConfig <- configureMCMC(myModel)

    
    myConfig$addMonitors('R0')
    
    if (modelType == 'full') {
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('beta', 'Z', 'x0C', 'x0D', 'nuC', 'nuD')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
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
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('Z', 'x0C', 'x0D', 'nuC', 'nuD')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    } else if (modelType == 'inc') { 
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('deltaC', 'x0C', 'nuC')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    } else if (modelType == 'fullNoAlarm') {
        
        # samplers for RstarI and RstarH
        myConfig$removeSamplers('RstarI') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarI'),
                            type = "RstarUpdate")
        
        myConfig$removeSamplers('RstarH') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarH'),
                            type = "RstarUpdate")
        
        myConfig$addMonitors(c('RstarI', 'RstarH'))
        
        # use slice sampling for rate parameters
        paramsForSlice <- c('beta', 'gamma1', 'gamma2', 'lambda','phi')
        myConfig$removeSampler(paramsForSlice)
        for (j in 1:length(paramsForSlice)) {
            myConfig$addSampler(target = paramsForSlice[j], type = "slice")
        }
        
    } else if (modelType == 'simpleNoAlarm') {
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    } 
    
    # monitor alarm functions when present
    if (modelType %in% c('simple', 'full')) {
        myConfig$addMonitors(c('yAlarmC', 'yAlarmD', 'alarmC', 'alarmD', 'delta'))
    } else if (modelType == 'inc') {
        myConfig$addMonitors(c('yAlarmC', 'alarmC', 'deltaC'))
    }
    
    if (assumeType == 'undetected') {
        
        myConfig$removeSampler('probDetect')
        myConfig$addSampler(target = 'probDetect', type = "slice")
        
        myConfig$removeSamplers('Istar') # Nodes will be expanded
        myConfig$addSampler(target = c('Istar'),
                            type = "RstarUpdate",
                            control = list(nUpdates = 1000))
        myConfig$addMonitors(c('Istar'))
        
    }
    
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


