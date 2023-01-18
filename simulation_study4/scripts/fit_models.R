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
  if (modelType != 'inc') {
    
    myConfig$addMonitors(c('yAlarmC', 'yAlarmD', 'alarmC', 'alarmD', 
                           'R0', 'delta'))
    
    if (modelType %in% c('simple', 'full')) {
      
      # use slice sampling for hill parameters
      paramsForSlice <- c('Z', 'x0C', 'x0D', 'nuC', 'nuD')
      myConfig$removeSampler(paramsForSlice)
      myConfig$addSampler(target = paramsForSlice, 
                          type = "RW_block",
                          control = list(propCov = diag(c(0.5, 0.5, 
                                                          100, 10, 
                                                          2, 2))))
      
    } else if (modelType %in% c('simpleThresh', 'fullThresh')) {

      # use slice sampling for threshold parameters
      paramsForSlice <- c('Z', 'HC', 'HD')
      myConfig$removeSampler(paramsForSlice)
      myConfig$addSampler(target = paramsForSlice, 
                          type = "RW_block",
                          control = list(propCov = diag(c(0.5, 0.5, 
                                                          50/modelInputs$constantsList$N,
                                                          15/modelInputs$constantsList$N))))
      
    }
    
  } else { # if model == 'inc'
    
    myConfig$addMonitors(c('yAlarmC','alarmC', 'R0', 'deltaC'))
    
    # use slice sampling for hill parameters
    paramsForSlice <- c('deltaC', 'x0C', 'nuC')
    myConfig$removeSampler(paramsForSlice)
    myConfig$addSampler(target = paramsForSlice, type = "AF_slice")
    
  }
  
  if (modelType %in% c('simple', 'inc')) {
    
    # joint sampler for beta and w0
    myConfig$removeSampler(c('beta', 'w0'))
    myConfig$addSampler(target = c('beta', 'w0'), type = "AF_slice")
  }
  
  if (modelType %in% c('full', 'fullThresh')) {
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
    
  }
  
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


