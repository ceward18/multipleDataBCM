################################################################################
# function to do posterior predictive forecasting of future epidemic curve 
# starting from last observed time point
################################################################################

# can only run if modelType %in% 'SIHRD_full', 'SIHRD_inc', 'SIR_inc', 'SIHRD_noAlarm', 'SIR_noAlarm')

postPredForecast <- function(incData, modelType, peak,
                             smoothC, smoothD, hospData, deathData, 
                             paramsSamples, N, S0, I0, H0, D0, R0,
                             Istar0, Dstar0) {
    
    # model-specific constants, data, and inits
    modelInputs <- getModelInput(incData, modelType, peak,
                                 smoothC, smoothD,
                                 hospData, deathData,
                                 N, S0, I0, H0, D0, R0)
    
    # set up timings
    obsTime <- modelInputs$constantsList$tau
    nDaysSim <- 50
    predTime <- obsTime + 1:nDaysSim
    
    modelInputs$constantsList$bw <- 30
    
    modelInputs$constantsList$smoothC0 <- smoothC[1]
    modelInputs$constantsList$Istar0 <- Istar0
    modelInputs$constantsList$Istar0Length <- length(Istar0)
    
    if (modelType %in% c('SIHRD_full')) {
        modelInputs$constantsList$smoothD0 <- smoothD[1]
        modelInputs$constantsList$Dstar0 <- Dstar0
        modelInputs$constantsList$Dstar0Length <- length(Dstar0)
    }
    
    # change tau
    modelInputs$constantsList$tau <- max(predTime)
    
    # get model code
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIR_inc')) {
        modelCode <- get(paste0(modelType, '_sim'))
        
    } else if (modelType %in% c('SIHRD_noAlarm', 'SIR_noAlarm')){
        # don't need separate code to simulate from these models, as they 
        #   have an alarm function which depends on epidemic trajectory
        modelCode <- get(modelType)
    }
    
    # extend data for modeling
    modelInputs$dataList$Istar <- c(modelInputs$dataList$Istar,
                                    rep(NA, nDaysSim))
    
    if (grepl('SIHRD', modelType)) {
        modelInputs$dataList$Hstar <- c(modelInputs$dataList$Hstar,
                                        rep(NA, nDaysSim))
        modelInputs$dataList$Dstar <- c(modelInputs$dataList$Dstar,
                                        rep(NA, nDaysSim))
        
    }
    
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIR_inc')) {
        modelInputs$dataList$smoothC <- NULL
        
        if (modelType == 'SIHRD_full') {
            modelInputs$dataList$smoothD <- NULL
        }
    }
    
    # compile model and simulator
    myModelPred <- nimbleModel(modelCode, 
                               constants = modelInputs$constantsList,
                               data = modelInputs$dataList )
    
    compiledPred  <- compileNimble(myModelPred) 
    
    # add hospData, and deathData if SIHRD to paramsSamples as these are observed at first
    
    if (grepl('SIHRD', modelType)) {
        dataNodes <- c(paste0('Istar[', predTime, ']'), 
                       paste0('Hstar[', predTime, ']'),
                       paste0('Dstar[', predTime, ']'), 
                       paste0('RstarH[', predTime, ']'),
                       paste0('RstarI[', predTime, ']'))
        
        HstarObs <- matrix(hospData, nrow = nrow(paramsSamples), ncol = length(hospData), byrow = T)
        colnames(HstarObs) <- paste0('Hstar[', 1:length(hospData), ']')
        
        DstarObs <- matrix(deathData, nrow = nrow(paramsSamples), ncol = length(deathData), byrow = T)
        colnames(DstarObs) <- paste0('Dstar[', 1:length(deathData), ']')
        
        
        paramsSamples <- cbind(paramsSamples, HstarObs, DstarObs)
        
    } else {
        dataNodes <- paste0('Istar[', predTime, ']')
    }
    
    IstarObs <- matrix(modelInputs$dataList$Istar[1:obsTime], nrow = nrow(paramsSamples), 
                       ncol = length(incData), byrow = T)
    colnames(IstarObs) <- paste0('Istar[', 1:length(incData), ']')
    
    
    paramsSamples <- cbind(paramsSamples, IstarObs)
    
    sim_R <- simulator(myModelPred, dataNodes)
    sim_C <- compileNimble(sim_R)
    
    # get order of parameters
    parentNodes <- myModelPred$getParents(dataNodes, stochOnly = TRUE)
    parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
    parentNodes <- myModelPred$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
    
    nPost <- 10000
    postPredInc <- matrix(NA, nrow = nDaysSim, ncol = nPost)
    postPredHosp <- matrix(NA, nrow = nDaysSim, ncol = nPost)
    postPredDeath <- matrix(NA, nrow = nDaysSim, ncol = nPost)
    set.seed(1)
    for (j in 1:nPost) {
        
        postIdx <- sample(1:nrow(paramsSamples), 1)
        
        trueVals <- paramsSamples[postIdx, parentNodes]
        
        postPredAll <- sim_C$run(trueVals, 1)
        
        postPredInc[,j] <- postPredAll[grep('^Istar', dataNodes)]
        
        if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIHRD_noAlarm')) {
            postPredHosp[,j] <- postPredAll[grep('Hstar', dataNodes)]
            postPredDeath[,j] <- postPredAll[grep('Dstar', dataNodes)]
        }
        
    }
    
    rbind(postPredInc, postPredHosp, postPredDeath)
    
}
