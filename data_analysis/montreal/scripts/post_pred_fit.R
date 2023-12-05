################################################################################
# function to do posterior predictive distribution of observed epidemic curve 
################################################################################

# can only run if modelType %in% 'SIHRD_full', 'SIR_inc', 'SIHRD_noAlarm', 'SIR_noAlarm')

postPredFit <- function(incData, modelType, assumeType, peak,
                        smoothC, smoothD, hospData, deathData, 
                        paramsSamples, N, S0, I0, H0, D0, R0,
                        Istar0, Dstar0) {
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, assumeType, peak,
                                 smoothC, smoothD,
                                 hospData, deathData,
                                 N, S0, I0, H0, D0, R0)
    
    modelInputs$constantsList$bw <- 30
    
    modelInputs$constantsList$smoothC0 <- smoothC[1]
    modelInputs$constantsList$Istar0 <- Istar0
    modelInputs$constantsList$Istar0Length <- length(Istar0)
    
    if (modelType %in% c('SIHRD_full', 'SIHRD_noAlarm')) {
        modelInputs$constantsList$smoothI0 <- smoothD[1]
        modelInputs$constantsList$Dstar0 <- Dstar0
        modelInputs$constantsList$Dstar0Length <- length(Dstar0)
    }
   
    # get model code
    if (modelType %in% c('SIHRD_full', 'SIR_inc')) {
        modelCode <- get(paste0(modelType, '_', assumeType, '_sim'))
    } else if (modelType %in% c('SIHRD_noAlarm', 'SIR_noAlarm')){
        # don't need separate code to simulate from these models, as they 
        #   have an alarm function which depends on epidemic trajectory
        modelCode <- get(paste0(modelType, '_', assumeType))
    }
    
    # compile model and simulator
    myModelPred <- nimbleModel(modelCode, 
                               constants = modelInputs$constantsList)
    
    compiledPred  <- compileNimble(myModelPred) 
    
    tau <- modelInputs$constantsList$tau
    
    if (modelType %in% c('SIHRD_full', 'SIHRD_noAlarm')) {
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
        
        if (modelType %in% c('SIHRD_full', 'SIHRD_noAlarm')) {
            postPredHosp[,j] <- postPredAll[grep('Hstar', dataNodes)]
            postPredDeath[,j] <- postPredAll[grep('Dstar', dataNodes)]
        }
        
        
        if (assumeType == 'undetected') {
            postPredCases[,j] <- postPredAll[grep('detectIstar', dataNodes)]
        }
    }
    
    rbind(postPredInc, postPredCases, postPredHosp, postPredDeath)
    
}
