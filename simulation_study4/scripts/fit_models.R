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
    
    # need to ensure all stochastic nodes are monitored for WAIC calculation
    myConfig$addMonitors(c('yAlarmC', 'yAlarmD', 'alarmC', 'alarmD', 
                           'R0', 'delta'))
    
    # # use slice sampling for hill parameters
    paramsForSlice <- c('Z', 'x0C', 'x0D', 'nuC', 'nuD')
    myConfig$removeSampler(paramsForSlice)
    myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
    
    # joint slice sampler for Z
    # myConfig$removeSampler(c('Z'))
    # myConfig$addSampler(target = c('Z'), type = "AF_slice")
    
    if (modelType == 'simple') {
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
    }
    
    if (modelType == 'full') {
        # samplers for RstarI and RstarH
        myConfig$removeSamplers('RstarI') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarI'),
                            type = "RstarUpdate")
        
        myConfig$removeSamplers('RstarH') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarH'),
                            type = "RstarUpdate")
        
        myConfig$addMonitors(c('RstarI', 'RstarH'))
        
        # use slice sampling for hill parameters
        paramsForSlice <- c('gamma1', 'gamma2', 'lambda','phi')
        myConfig$removeSampler(paramsForSlice)
        for (j in 1:length(paramsForSlice)) {
            myConfig$addSampler(target = paramsForSlice[j], type = "slice")
        }
        
    }
    
    # browser()
    
    myMCMC <- buildMCMC(myConfig)
    compiled <- compileNimble(myModel, myMCMC) 
    
    runMCMC(compiled$myMCMC, 
            niter = niter, 
            nburnin = nburn,
            thin = nthin,
            setSeed = seed)
    
}


