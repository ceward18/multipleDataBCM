################################################################################
# Calculate WAIC from matrix of posterior samples and model object
# need to recompile model outside of the parallel environment that was used 
#   to run multiple chains
# called from summarize post after samples have been combined across chains
################################################################################

getWAIC <- function(samples, modelType, incData, smoothC, smoothD,
                    hospData, deathData) {
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, smoothC, smoothD,
                                 hospData, deathData)
    
    ### get appropriate model code
    modelCode <- get(paste0('SIHRD_', modelType))
    
    ### create nimble model
    myModel <- nimbleModel(modelCode, 
                           data = modelInputs$dataList, 
                           constants = modelInputs$constantsList,
                           inits = modelInputs$initsList)
    
    compiled <- compileNimble(myModel) 
    
    
    if (modelType == 'full') {
        
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
