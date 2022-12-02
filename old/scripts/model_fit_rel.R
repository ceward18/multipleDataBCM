################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, smoothC, smoothH, smoothD, N, I0, R0, seed) {
    
    source('./scripts/model_codes.R')
    source('./scripts/get_model_inputs_rel.R')
    
    # for reproducibility so inits are always the same
    set.seed(seed + 3)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, smoothC, smoothH, smoothD, N, I0, R0)
    
    ### MCMC specifications
    niter <- modelInputs$niter
    nburn <- modelInputs$nburn
    nthin <- modelInputs$nthin
    
    ### create nimble model
    myModel <- nimbleModel(SIR_gp_rel, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    myConfig <- configureMCMC(myModel)
    
    # need to ensure all stochastic nodes are monitored for WAIC calculation
    myConfig$addMonitors(c('yC', 'yH', 'yD', 'alarm', 'R0'))
    
    # if gaussian process model, use slice sampling
    paramsForSlice <- c('sigma', 'lC', 'lH', 'lD')
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


