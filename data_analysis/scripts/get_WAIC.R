################################################################################
# Calculate WAIC from matrix of posterior samples and model object
# need to recompile model outside of the parallel environment that was used 
#   to run multiple chains
# called from summarize post after samples have been combined across chains
################################################################################

getWAIC <- function(samples, city, modelType, peak,
                    smoothC, smoothD,
                    hospData, deathData,
                    N, S0, I0, H0, D0, R0) {
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, city, modelType, peak,
                                 smoothC, smoothD,
                                 hospData, deathData,
                                 N, S0, I0, H0, D0, R0)
    
    ### get appropriate model code
    modelCode <- get(modelType)
    
    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    
    compiled <- compileNimble(myModel) 
    
    Istar <- modelInputs$dataList$Istar
    incDataSamples <- matrix(rep(Istar), NROW(samples),
                             ncol = length(Istar), byrow = T)
    colnames(incDataSamples) <- paste0('Istar[', 1:ncol(incDataSamples), ']')
    samples <- cbind(samples, incDataSamples)
    
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIHRD_noAlarm')) {
        
        hospDataSamples <- matrix(rep(hospData), NROW(samples),
                                 ncol = length(hospData), byrow = T)
        colnames(hospDataSamples) <- paste0('Hstar[', 1:ncol(hospDataSamples), ']')
        samplesHstar <- cbind(samples, hospDataSamples)
        
        deathDataSamples <- matrix(rep(deathData), NROW(samples),
                                 ncol = length(deathData), byrow = T)
        colnames(deathDataSamples) <- paste0('Dstar[', 1:ncol(deathDataSamples), ']')
        samplesDstar <- cbind(samplesHstar, deathDataSamples)
        
        waicList <- calculateWAIC(samplesDstar, compiled)
        
    } else {
        
        waicList <- calculateWAIC(samples, compiled)
    }
    
    
    
    data.frame(waic = waicList$WAIC,
               lppd = waicList$lppd,
               pWAIC = waicList$pWAIC)
    
}
