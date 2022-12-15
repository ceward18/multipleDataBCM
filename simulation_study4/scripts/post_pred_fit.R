################################################################################
# function to do posterior predictive distribution of observed epidemic curve 
################################################################################

# can only run if modelType == 'full'

postPredFit <- function(incData, smoothC, smoothD, hospData, deathData, 
                        paramsSamples) {
    
    modelType <- 'full'
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, smoothC, smoothD,
                                 hospData, deathData)
    
    modelInputs$constantsList$bw <- 30
    
    # compile model and simulator
    myModelPred <- nimbleModel(SIHRD_full_sim, 
                               constants = modelInputs$constantsList)
    
    compiledPred  <- compileNimble(myModelPred) 
    
    tau <- modelInputs$constantsList$tau
    dataNodes <- c('Istar', 'Hstar', 'Dstar', 'RstarH', 'RstarI')
    dataNodes <- myModelPred$expandNodeNames(dataNodes)
    
    sim_R <- simulator(myModelPred, dataNodes)
    sim_C <- compileNimble(sim_R)
    
    # get order of parameters
    parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
    parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
    parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
    
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
        postPredHosp[,j] <- postPredAll[grep('Hstar', dataNodes)]
        postPredDeath[,j] <- postPredAll[grep('Dstar', dataNodes)]
    }
    
    rbind(postPredInc, postPredHosp, postPredDeath)

}
