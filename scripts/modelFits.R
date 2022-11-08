################################################################################
# function to fit models 
################################################################################

fitAlarmModel <- function(incData, smoothC, smoothH, smoothD, alarmFit,
                          N, I0, R0, Rstar0, lengthI,seed) {
  
  source('./scripts/modelCodes.R')
  source('./scripts/getModelInputs.R')
  
  # get appropriate model code
  modelCode <- get(paste0('SIR_', alarmFit, '_fixed'))
  
  # for reproducibility so inits are always the same
  set.seed(seed + 3)

  # model-specific constants, data, and inits
  modelInputs <- getModelInput(alarmFit, incData, smoothC, smoothH, smoothD,
                               infPeriod = 'fixed', 
                               N, I0, R0, Rstar0, lengthI)
  
  ### MCMC specifications
  niter <- modelInputs$niter
  nburn <- modelInputs$nburn
  nthin <- modelInputs$nthin
  
  ### create nimble model
  myModel <- nimbleModel(modelCode, 
                         data = modelInputs$dataList, 
                         constants = modelInputs$constantsList,
                         inits = modelInputs$initsList)
  myConfig <- configureMCMC(myModel)
  
  # need to ensure all stochastic nodes are monitored for WAIC calculation
  myConfig$addMonitors(c('yC', 'yH', 'yD', 
                         'fCases', 'fHosp', 'fDeath', 'alarm', 'R0'))
  
  # if gaussian process model, use slice sampling
  if (alarmFit == 'gp') {

    paramsForSlice <- c('beta', 'lC', 'sigmaC', 'lH', 'sigmaH', 'lD', 'sigmaD')
    myConfig$removeSampler(paramsForSlice)
    for (j in 1:length(paramsForSlice)) {
        myConfig$addSampler(target = paramsForSlice[j], type = "slice")
    }

  }
  
  myMCMC <- buildMCMC(myConfig)
  compiled <- compileNimble(myModel, myMCMC) 
  
  runMCMC(compiled$myMCMC, 
          niter = niter, 
          nburnin = nburn,
          thin = nthin,
          setSeed  = seed)
  
}




