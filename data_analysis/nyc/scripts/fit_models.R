################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, modelType, assumeType, peak,
                          smoothC, smoothD, hospData, deathData, 
                          N, S0, I0, H0, D0, R0, seed) {
    
    
    source('./scripts/model_code.R')
    source('./scripts/get_model_inputs.R')
    
    # for reproducibility, inits vary across peaks/models/assumptions/chains
    seed <- seed + peak + 
        as.numeric(factor(modelType, levels = c('SIHRD_full', 'SIHRD_inc',
                                                'SIR_full', 'SIR_inc',
                                                'SIHRD_noAlarm', 'SIR_noAlarm'))) + 
        as.numeric(factor(assumeType, levels = c('undetected', 'casesOnly')))
    
    set.seed(seed)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, assumeType, peak,
                                 smoothC, smoothD,
                                 hospData, deathData,
                                 N, S0, I0, H0, D0, R0)
    
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
    peak# browser()
    myConfig <- configureMCMC(myModel)
    
    # track reproductive number
    myConfig$addMonitors('R0')
    
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc')) {

        # samplers for RstarI and RstarH
        myConfig$removeSamplers('RstarI') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarI'),
                            type = "RstarUpdate")
        
        myConfig$removeSamplers('RstarH') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarH'),
                            type = "RstarUpdate")
        
        # need to ensure all stochastic nodes are monitored for WAIC calculation
        myConfig$addMonitors(c('RstarI', 'RstarH'))
        
    } else if (modelType == 'SIR_full') {
        
        # use slice sampling for transmission parameters
        paramsForSlice <- c('beta', 'k', 'w0')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
    } else if (modelType == 'SIR_inc') { 
        
        # use slice sampling for transmission parameters
        paramsForSlice <- c('beta', 'k', 'w0')
        myConfig$removeSampler(paramsForSlice)
        myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
        
    } else if (modelType == 'SIHRD_noAlarm') {
        
        # samplers for RstarI and RstarH
        myConfig$removeSamplers('RstarI') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarI'),
                            type = "RstarUpdate")
        
        myConfig$removeSamplers('RstarH') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarH'),
                            type = "RstarUpdate")
        
        # need to ensure all stochastic nodes are monitored for WAIC calculation
        myConfig$addMonitors(c('RstarI', 'RstarH'))
        
    } else if (modelType == 'SIR_noAlarm') {
        
        # joint sampler for beta and w0
        myConfig$removeSampler(c('beta', 'w0'))
        myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
        
    } 
   
    # monitor alarm functions when present
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIR_full', 'SIR_inc')) {
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


