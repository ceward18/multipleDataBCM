################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, alarmBase, N, I0, R0, seed) {
    
    source('./scripts/model_codes.R')
    source('./scripts/get_model_inputs_uni.R')
    
    # for reproducibility so inits are always the same
    set.seed(seed + 3)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, alarmBase, N, I0, R0)
    
    ### MCMC specifications
    niter <- modelInputs$niter
    nburn <- modelInputs$nburn
    nthin <- modelInputs$nthin
    
    ### create nimble model
    myModel <- nimbleModel(SIR_gp_uni, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    myConfig <- configureMCMC(myModel)
    
    # need to ensure all stochastic nodes are monitored for WAIC calculation
    myConfig$addMonitors(c('yAlarm', 'alarm', 'R0'))
    
    # use slice sampling for GP parameters
    paramsForSlice <- c('l', 'sigma')
    myConfig$removeSampler(paramsForSlice)
    for (j in 1:length(paramsForSlice)) {
        myConfig$addSampler(target = paramsForSlice[j], type = "slice")
    }
    
    # joint sampler for beta and w0
    myConfig$removeSampler(c('beta', 'w0'))
    myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
    
    myMCMC <- buildMCMC(myConfig)
    compiled <- compileNimble(myModel, myMCMC) 
    
    runMCMC(compiled$myMCMC, 
            niter = niter, 
            nburnin = nburn,
            thin = nthin,
            setSeed  = seed)
    
}




