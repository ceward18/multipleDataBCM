################################################################################
# function to do posterior predictive distribution of observed epidemic curve 
################################################################################

# can only run if modelType == 'full' or 'inc'

postPredFit <- function(incData, modelType, 
                        smoothC, smoothD, hospData, deathData, 
                        paramsSamples, N, S0, I0, H0, D0, R0,
                        Istar0, Dstar0) {
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, smoothC, smoothD,
                                 hospData, deathData,
                                 N, S0, I0, H0, D0, R0)
    
    modelInputs$constantsList$bw <- 30
    
    modelInputs$constantsList$smoothC0 <- smoothC[1]
    modelInputs$constantsList$Istar0 <- Istar0
    modelInputs$constantsList$Istar0Length <- length(Istar0)
    
    if (modelType == 'full') {
        modelInputs$constantsList$smoothI0 <- smoothD[1]
        modelInputs$constantsList$Dstar0 <- Dstar0
        modelInputs$constantsList$Dstar0Length <- length(Dstar0)
    }
   
    # get model code
    modelCode <- get(paste0('SIHRD_', modelType, '_sim'))
    
    # compile model and simulator
    myModelPred <- nimbleModel(modelCode, 
                               constants = modelInputs$constantsList)
    
    compiledPred  <- compileNimble(myModelPred) 
    
    tau <- modelInputs$constantsList$tau
    
    if (modelType %in% c('full')) {
        dataNodes <- c('Istar', 'Hstar', 'Dstar', 'RstarH', 'RstarI')
    } else {
        dataNodes <- 'Istar'
    }
    dataNodes <- myModelPred$expandNodeNames(dataNodes)
    
    sim_R <- simulator(myModelPred, dataNodes)
    sim_C <- compileNimble(sim_R)
    
    # get order of parameters
    parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
    parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
    parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
    
    # add detectIStar to paramsSamples
    detectIstarSamples <- matrix(incData, ncol = tau, nrow = nrow(paramsSamples), byrow = T)
    colnames(detectIstarSamples) <- paste0('detectIstar[', 1:tau, ']')
    paramsSamples <- cbind(paramsSamples, detectIstarSamples)
    
    nPost <- 10000
    postPredInc <- matrix(NA, nrow = tau, ncol = nPost)
    postPredHosp <- matrix(NA, nrow = tau, ncol = nPost)
    postPredDeath <- matrix(NA, nrow = tau, ncol = nPost)
    set.seed(1)
    for (j in 1:nPost) {
        
        postIdx <- sample(1:nrow(paramsSamples), 1)
        
        trueVals <- paramsSamples[postIdx, parentNodes]
        
        postPredAll <- apply(sim_C$run(trueVals, 10), 2, median)
        
        
        postPredInc[,j] <- postPredAll[grep('Istar', dataNodes)]
        
        if (modelType %in% c('full')) {
            postPredHosp[,j] <- postPredAll[grep('Hstar', dataNodes)]
            postPredDeath[,j] <- postPredAll[grep('Dstar', dataNodes)]
        }
    }
    
    rbind(postPredInc, postPredHosp, postPredDeath)
    
}
