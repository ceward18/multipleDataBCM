################################################################################
# function to do posterior predictive distribution of observed epidemic curve 
################################################################################

# can only run if modelType == 'full' or 'fullThresh' or 'inc'

postPredFit <- function(incData, modelType, assumeType,
                        smoothC, smoothD, hospData, deathData, 
                        paramsSamples) {
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, assumeType, smoothC, smoothD,
                                 hospData, deathData)
    
    modelInputs$constantsList$bw <- 30
    
    # get model code
    if (modelType %in% c('full', 'inc')) {
        modelCode <- get(paste0('SIHRD_', modelType, '_', assumeType, '_sim')) 
    } else if (modelType %in% c('fullNoAlarm', 'simpleNoAlarm')){
        # don't need separate code to simulate from these models, as they 
        #   have an alarm function which depends on epidemic trajectory
        modelCode <- get(paste0('SIHRD_', modelType, '_', assumeType))
    }
    
    
    # compile model and simulator
    myModelPred <- nimbleModel(modelCode, 
                               constants = modelInputs$constantsList)
    
    compiledPred  <- compileNimble(myModelPred) 
    
    tau <- modelInputs$constantsList$tau
    
    if (modelType %in% c('full', 'fullNoAlarm')) {
        dataNodes <- c('Istar', 'Hstar', 'Dstar', 'RstarH', 'RstarI')
    } else {
        dataNodes <- 'Istar'
    }
    if (assumeType == 'undetected') dataNodes <- c(dataNodes, 'detectIstar')
    
    dataNodes <- myModelPred$expandNodeNames(dataNodes)
    
    sim_R <- simulator(myModelPred, dataNodes)
    sim_C <- compileNimble(sim_R)
    
    # get order of parameters
    parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
    parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
    parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
    
    nPost <- 10000
    postPredInc <- matrix(NA, nrow = tau, ncol = nPost)
    postPredCases <- matrix(NA, nrow = tau, ncol = nPost)
    postPredHosp <- matrix(NA, nrow = tau, ncol = nPost)
    postPredDeath <- matrix(NA, nrow = tau, ncol = nPost)
    set.seed(1)
    for (j in 1:nPost) {
        
        postIdx <- sample(1:nrow(paramsSamples), 1)
        
        trueVals <- paramsSamples[postIdx, parentNodes]
        
        postPredAll <- apply(sim_C$run(trueVals, 10), 2, median)
        
        
        postPredInc[,j] <- postPredAll[grep('^Istar', dataNodes)]
        
        if (modelType %in% c('full', 'fullNoAlarm')) {
            postPredHosp[,j] <- postPredAll[grep('Hstar', dataNodes)]
            postPredDeath[,j] <- postPredAll[grep('Dstar', dataNodes)]
        }
        
        if (assumeType == 'undetected') {
            postPredCases[,j] <- postPredAll[grep('detectIstar', dataNodes)]
        }
    }
    
    rbind(postPredInc, postPredCases, postPredHosp, postPredDeath)
    
}
