################################################################################
# summarize posterior samples from three chains
# gelman-rubin to assess convergence
# summarize posterior means and credible intervals:
#   model parameters   
#   alarm function
# posterior prediction mean and credible intervals
################################################################################

library(coda)
library(nimble)

# source relevant scripts
source('./scripts/get_model_inputs.R')
source('./scripts/get_WAIC.R')
source('./scripts/post_pred_fit.R')

summarizePost <- function(resThree, incData, modelType,
                          smoothC, smoothH, smoothD, hospData, deathData) {
    
    paramSamples1 <- resThree[[1]][,-grep('alarm|R0|yAlarm|Rstar', colnames(resThree[[1]]))]
    paramSamples2 <- resThree[[2]][,-grep('alarm|R0|yAlarm|Rstar', colnames(resThree[[2]]))]
    paramSamples3 <- resThree[[3]][,-grep('alarm|R0|yAlarm|Rstar', colnames(resThree[[3]]))]
    
    ##############################################################################
    ### gelman-rubin
    
    res_mcmc <- mcmc.list(mcmc(paramSamples1),
                          mcmc(paramSamples2),
                          mcmc(paramSamples3))
    gdiag <- data.frame(gelman.diag(res_mcmc, multivariate = F)$psrf)
    colnames(gdiag) <- c('gr', 'grUpper')
    gdiag$param <- rownames(gdiag)
    rownames(gdiag) <- NULL
    
    ##############################################################################
    ### posterior mean and 95% CI for parameters
    
    paramsPost <- rbind(paramSamples1, paramSamples2, paramSamples3)
    postMeans <- colMeans(paramsPost)
    postCI <- apply(paramsPost, 2, quantile, probs = c(0.025, 0.975))
    postParams <- data.frame(param = names(postMeans),
                             mean = postMeans,
                             lower = postCI[1,],
                             upper = postCI[2,])
    rownames(postParams) <- NULL
    
    ##############################################################################
    ### posterior distribution of alarm function over range observed
    
    # get xC, xH, xD (range of observed incidence/hosp/deaths)
    modelInputs <- getModelInput(incData, modelType, smoothC, smoothH, smoothD,
                                 hospData, deathData)
    n <- length(modelInputs$xC)
    
    alarmSamples1 <- resThree[[1]][,grep('yAlarm', colnames(resThree[[1]]))]
    alarmSamples2 <- resThree[[2]][,grep('yAlarm', colnames(resThree[[2]]))]
    alarmSamples3 <- resThree[[3]][,grep('yAlarm', colnames(resThree[[3]]))]
    alarmSamples <- rbind(alarmSamples1, alarmSamples2, alarmSamples3)
    alarmSamples <- alarmSamples[,c(paste0('yAlarmC[', 1:n, ']'),
                                    paste0('yAlarmH[', 1:n, ']'),
                                    paste0('yAlarmD[', 1:n, ']'))]
    
    postMeans <- colMeans(alarmSamples)
    postCI <- apply(alarmSamples, 2, quantile, probs = c(0.025, 0.975))
    
    postAlarm <- data.frame(xAlarm = c(modelInputs$xC, modelInputs$xH, modelInputs$xD),
                            marg = rep(c('inc', 'hosp', 'death'), each = n),
                            mean = postMeans,
                            lower = postCI[1,],
                            upper = postCI[2,])
    
    rownames(postAlarm) <- NULL
    
    ##############################################################################
    ### posterior distribution of individual alarm functions
    
    tau <- length(incData)
    
    alarmSamples1 <- resThree[[1]][,grep('alarm[[:upper:]]', colnames(resThree[[1]]))]
    alarmSamples2 <- resThree[[2]][,grep('alarm[[:upper:]]', colnames(resThree[[2]]))]
    alarmSamples3 <- resThree[[3]][,grep('alarm[[:upper:]]', colnames(resThree[[3]]))]
    alarmSamples <- rbind(alarmSamples1, alarmSamples2, alarmSamples3)
    alarmSamples <- alarmSamples[,c(paste0('alarmC[', 1:tau, ']'),
                                    paste0('alarmH[', 1:tau, ']'),
                                    paste0('alarmD[', 1:tau, ']'))]
    
    postMeans <- colMeans(alarmSamples)
    postCI <- apply(alarmSamples, 2, quantile, probs = c(0.025, 0.975))
    
    postAlarmTime <- data.frame(time = rep(1:tau, 3),
                                marg = rep(c('inc', 'hosp', 'death'), each = tau),
                                mean = postMeans,
                                lower = postCI[1,],
                                upper = postCI[2,])
    
    rownames(postAlarmTime) <- NULL
    
    ##############################################################################
    ### Posterior distribution of reproductive number over epidemic time
    
    R0Samples1 <-  resThree[[1]][,grep('R0', colnames(resThree[[1]]))]
    R0Samples2 <-  resThree[[2]][,grep('R0', colnames(resThree[[2]]))]
    R0Samples3 <-  resThree[[3]][,grep('R0', colnames(resThree[[3]]))]
    
    # combine posterior parameters with posterior R0
    R0Samples <- rbind(R0Samples1, R0Samples2, R0Samples3)
    
    postMeans <- colMeans(R0Samples)
    postCI <- apply(R0Samples, 2, quantile, probs = c(0.025, 0.975))
    postR0 <- data.frame(time = 1:length(postMeans),
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postR0) <- NULL
    
    ##############################################################################
    ### WAIC values
    
    # samples to use for WAIC calculation differ by model
    samples <- rbind(resThree[[1]], resThree[[2]], resThree[[3]])
    
    waic <- getWAIC(samples = samples, modelType = modelType, incData = incData, 
                    smoothC = smoothC, smoothH = smoothH, smoothD = smoothD,
                    hospData = hospData, deathData = deathData)
    
    ##############################################################################
    ### Posterior predictive fit for full model
    
    if (modelType == 'full') {
        postPred <- postPredFit(incData = incData, smoothC = smoothC,  
                                smoothH = smoothH, smoothD = smoothD, 
                                hospData = hospData, deathData = deathData, 
                                paramsSamples = samples)
        
        postMeans <- rowMeans(postPred)
        postCI <- apply(postPred, 1, quantile, probs = c(0.025, 0.975))
        
        postPredictFit <- data.frame(time = rep(1:tau, 3),
                                     marg = rep(c('inc', 'hosp', 'death'), each = tau),
                                     mean = postMeans,
                                     lower = postCI[1,],
                                     upper = postCI[2,])
    } else {
        postPredictFit <- data.frame(time = NA, 
                                     marg = NA,
                                     mean = NA,
                                     lower = NA,
                                     upper = NA)
    }
    
    ##############################################################################
    ### output
    list(gdiag = gdiag,
         postParams = postParams,
         postAlarm = postAlarm,
         postAlarmTime = postAlarmTime,
         postR0 = postR0,
         waic= waic,
         postPredictFit = postPredictFit)
    
}





