################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, smoothC, smoothH, smoothD, N, I0, R0, seed) {
    
    source('./scripts/model_codes.R')
    source('./scripts/get_model_inputs_thresh.R')
    
    # for reproducibility so inits are always the same
    set.seed(seed + 3)
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, smoothC, smoothH, smoothD, N, I0, R0)
    
    ### MCMC specifications
    niter <- modelInputs$niter
    nburn <- modelInputs$nburn
    nthin <- modelInputs$nthin
    
    ### create nimble model
    myModel <- nimbleModel(SIR_thresh_multi, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    myConfig <- configureMCMC(myModel)
    
    # need to ensure all stochastic nodes are monitored for WAIC calculation
    myConfig$addMonitors(c('yAlarm', 'alarm', 'R0'))
 
    
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


