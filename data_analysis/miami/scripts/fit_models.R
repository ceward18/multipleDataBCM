################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, city, modelType, peak, probDetectMean,
                          smoothC, smoothD, hospData, deathData, 
                          N, S0, I0, H0, D0, R0, seed) {
    
    
    source('./scripts/model_code.R')
    source('./scripts/get_model_inputs.R')
    
    # for reproducibility, inits vary across peaks/models/chains/city/n_weeks
    seed <- seed + 
        as.numeric(peak) + 
        as.numeric(factor(modelType, levels = c('SIHRD_full', 'SIHRD_inc',
                                                'SIR_full', 'SIR_inc',
                                                'SIHRD_noAlarm', 'SIR_noAlarm'))) + 
        as.numeric(factor(city, levels = c('montreal', 'miami'))) + 
        length(incData)
    
    set.seed(seed)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, peak, probDetectMean,
                                 smoothC, smoothD,
                                 hospData, deathData,
                                 N, S0, I0, H0, D0, R0)
    
    ### MCMC specifications
    niter <- modelInputs$niter
    nburn <- modelInputs$nburn
    nthin <- modelInputs$nthin
    
    ### get appropriate model code
    modelCode <- get(modelType)
    
    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    
    # browser()
    myConfig <- configureMCMC(myModel)
    
    # track reproductive number
    myConfig$addMonitors('R0')
    
    if (grepl('SIHRD', modelType)) {
        
        # samplers for RstarI and RstarH
        myConfig$removeSamplers('RstarI') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarI'),
                            type = "RstarUpdate")
        
        myConfig$removeSamplers('RstarH') # Nodes will be expanded
        myConfig$addSampler(target = c('RstarH'),
                            type = "RstarUpdate")
        
        # impute missing hospitalizations during peak 1
        if (peak == 1) {
            
            myConfig$removeSamplers('Hstar') # Nodes will be expanded
            myConfig$addSampler(target = c('Hstar'),
                                type = "HstarUpdate",
                                control = list(fixed = which(is.na(hospData))))
        }
        
        # need to ensure all stochastic nodes are monitored for WAIC calculation
        myConfig$addMonitors(c('RstarI', 'RstarH', 'Hstar'))
    
        
        # use joint sampling for transmission parameters
        jointParams <- c('beta', 'gamma1', 'lambda', 'gamma2', 'phi')
        # add k and alpha where appropriate
        if (!grepl('noAlarm', modelType)) {
            jointParams <- c(jointParams, 'k')
        } 
        if (grepl('full', modelType)) {
            jointParams <- c(jointParams, 'alpha')
        } 
        myConfig$removeSampler(jointParams)
        myConfig$addSampler(target = jointParams, type = "AF_slice")
        
    } else if (grepl('SIR', modelType)) {
        
        # use slice sampling for transmission parameters
        jointParams <- c('beta', 'w0')
        if (!grepl('noAlarm', modelType)) {
            jointParams <- c(jointParams, 'k')
        } 
        if (grepl('full', modelType)) {
            jointParams <- c(jointParams, 'alpha')
        } 
        myConfig$removeSampler(jointParams)
        myConfig$addSampler(target = jointParams, type = "AF_slice")
        
    } 
    
    # monitor alarm functions when present
    if (!grepl('noAlarm', modelType)) {
        myConfig$addMonitors('alarm')
    } 
    
    # data augmentation for undetected infections
    myConfig$removeSamplers('Istar') # Nodes will be expanded
    myConfig$addSampler(target = c('Istar'),
                        type = "RstarUpdate",
                        control = list(nUpdates = 1500))
    myConfig$addMonitors('Istar')
    
    print(myConfig)
    
    # browser()
    nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
    myMCMC <- buildMCMC(myConfig)
    compiled <- compileNimble(myModel, myMCMC) 
    
    runMCMC(compiled$myMCMC, 
            niter = niter, 
            nburnin = nburn,
            thin = nthin)
    
}


