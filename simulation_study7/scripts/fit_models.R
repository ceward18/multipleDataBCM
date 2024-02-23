################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, modelType, assumeType,  
                          smoothC, smoothD, hospData, deathData, seed) {
    
    source('./scripts/model_code.R')
    source('./scripts/get_model_inputs.R')
    
    # for reproducibility so inits are always the same
    set.seed(seed + 3)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, assumeType,
                                 smoothC, smoothD, hospData, deathData)
    
    ### MCMC specifications
    niter <- modelInputs$niter
    nburn <- modelInputs$nburn
    nthin <- modelInputs$nthin
    
    ### get appropriate model code
    modelCode <- get(paste0(modelType, '_', assumeType))
    
    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)

    myConfig <- configureMCMC(myModel)

    myConfig$addMonitors('R0')
    
    if (modelType == 'SIHRD_full') {
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('beta', 'delta', 'x0', 'nu')
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
        paramsForSlice <- c('gamma1', 'gamma2', 'lambda', 'phi')
        myConfig$removeSampler(paramsForSlice)
        for (j in 1:length(paramsForSlice)) {
            myConfig$addSampler(target = paramsForSlice[j], type = "slice")
        }
        
    } else if (modelType == 'SIR_full') {
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('delta', 'x0', 'nu')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    } else if (modelType == 'SIR_inc') { 
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('delta', 'x0', 'nu')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    } else if (modelType == 'SIHRD_noAlarm') {
        
        # samplers for RstarI and RstarH
        myConfig$removeSamplers('RstarI') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarI'),
                            type = "RstarUpdate")
        
        myConfig$removeSamplers('RstarH') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarH'),
                            type = "RstarUpdate")
        
        myConfig$addMonitors(c('RstarI', 'RstarH'))
        
        # use slice sampling for rate parameters
        paramsForSlice <- c('beta', 'gamma1', 'gamma2', 'lambda', 'phi')
        myConfig$removeSampler(paramsForSlice)
        for (j in 1:length(paramsForSlice)) {
            myConfig$addSampler(target = paramsForSlice[j], type = "slice")
        }
        
    } else if (modelType == 'SIR_noAlarm') {
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    } 
    
    # monitor alarm functions when present
    if (modelType %in% c('SIHRD_full', 'SIR_full', 'SIR_inc')) {
        myConfig$addMonitors(c('alarm'))
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


