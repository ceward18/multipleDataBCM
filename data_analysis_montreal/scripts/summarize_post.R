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
source('./scripts/model_code.R')
source('./scripts/get_model_inputs.R')
source('./scripts/get_WAIC.R')
source('./scripts/post_pred_fit.R')
source('./scripts/post_pred_forecast.R')

summarizePost <- function(resThree, incData, modelType, peak,
                          smoothC, smoothD, hospData, deathData,
                          N, S0, I0, H0, D0, R0, Istar0, Dstar0) {
    
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIHRD_noAlarm')) {
        paramSamples1 <- resThree[[1]][,-grep('alarm|R0|Rstar|
                                              |comp_init\\[3\\]|comp_init\\[4\\]|comp_init\\[5\\]', colnames(resThree[[1]]))]
        paramSamples2 <- resThree[[2]][,-grep('alarm|R0|Rstar|
                                              |comp_init\\[3\\]|comp_init\\[4\\]|comp_init\\[5\\]', colnames(resThree[[2]]))]
        paramSamples3 <- resThree[[3]][,-grep('alarm|R0|Rstar|
                                              |comp_init\\[3\\]|comp_init\\[4\\]|comp_init\\[5\\]', colnames(resThree[[3]]))]
    } else if (modelType %in% c('SIR_full', 'SIR_inc', 'SIR_noAlarm')) {
        paramSamples1 <- resThree[[1]][,-grep('alarm|R0|comp_init\\[16\\]', colnames(resThree[[1]]))]
        paramSamples2 <- resThree[[2]][,-grep('alarm|R0|comp_init\\[16\\]', colnames(resThree[[2]]))]
        paramSamples3 <- resThree[[3]][,-grep('alarm|R0|comp_init\\[16\\]', colnames(resThree[[3]]))]
    }
    
    ##############################################################################
    ### gelman-rubin
    print('Calculating GR diagnostic...')
    
    res_mcmc <- mcmc.list(mcmc(paramSamples1),
                          mcmc(paramSamples2),
                          mcmc(paramSamples3))
    gdiag <- data.frame(gelman.diag(res_mcmc, multivariate = F)$psrf)
    colnames(gdiag) <- c('gr', 'grUpper')
    gdiag$param <- rownames(gdiag)
    rownames(gdiag) <- NULL
    
    ##############################################################################
    ### posterior mean and 95% CI for parameters
    print('Summarizing posterior distributions...')
    
    paramsPost <- rbind(paramSamples1, paramSamples2, paramSamples3)
    postMeans <- colMeans(paramsPost)
    postCI <- apply(paramsPost, 2, quantile, probs = c(0.025, 0.975))
    postParams <- data.frame(param = names(postMeans),
                             mean = postMeans,
                             lower = postCI[1,],
                             upper = postCI[2,])
    rownames(postParams) <- NULL
    
    ##############################################################################
    ### posterior distribution of alarm function
    
    tau <- length(incData)
    
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIR_full', 'SIR_inc')) {
        
        alarmSamples1 <- resThree[[1]][,grep('alarm', colnames(resThree[[1]]))]
        alarmSamples2 <- resThree[[2]][,grep('alarm', colnames(resThree[[2]]))]
        alarmSamples3 <- resThree[[3]][,grep('alarm', colnames(resThree[[3]]))]
        alarmSamples <- rbind(alarmSamples1, alarmSamples2, alarmSamples3)
        postMeans <- colMeans(alarmSamples)
        postCI <- apply(alarmSamples, 2, quantile, probs = c(0.025, 0.975))
        
        postAlarmTime <- data.frame(time = 1:tau,
                                    mean = postMeans,
                                    lower = postCI[1,],
                                    upper = postCI[2,])
        
        
        rownames(postAlarmTime) <- NULL
        
    } else if (modelType %in% c('SIHRD_noAlarm', 'SIR_noAlarm')) {
        
        postAlarmTime <- data.frame(time = NA,
                                    mean = NA,
                                    lower = NA,
                                    upper = NA)
    }
    rownames(postAlarmTime) <- NULL
    
    ##############################################################################
    ### posterior distribution of R0
    
    R0Samples1 <- resThree[[1]][,grep('R0', colnames(resThree[[1]]))]
    R0Samples2 <- resThree[[2]][,grep('R0', colnames(resThree[[2]]))]
    R0Samples3 <- resThree[[3]][,grep('R0', colnames(resThree[[3]]))]
    R0Samples <- rbind(R0Samples1, R0Samples2, R0Samples3)
    
    postMeans <- colMeans(R0Samples)
    postCI <- apply(R0Samples, 2, quantile, probs = c(0.025, 0.975))
    
    postR0 <- data.frame(time = 1:length(postMeans),
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    
    
    ##############################################################################
    ### WAIC values
    print('Calculating WAIC...')
    
    # samples to use for WAIC calculation differ by model
    samples <- rbind(resThree[[1]], resThree[[2]], resThree[[3]])
    
    waic <- getWAIC(samples = samples, modelType = modelType, peak = peak,
                    smoothC = smoothC, smoothD = smoothD,
                    hospData = hospData, deathData = deathData,
                    N = N, S0 = S0, I0 = I0, H0 = H0, D0 = D0, R0 = R0)
    
    
    ##############################################################################
    ### Posterior predictive fit for models where it makes sense
    print('Calculating posterior predictive fit...')
    
    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIR_inc', 'SIHRD_noAlarm', 'SIR_noAlarm')) {
        postPred <- postPredFit(incData = incData, 
                                modelType = modelType, peak = peak,
                                smoothC = smoothC, smoothD = smoothD, 
                                hospData = hospData, deathData = deathData, 
                                paramsSamples = samples,
                                N = N, S0 = S0, I0 = I0, H0 = H0, D0 = D0, R0 = R0,
                                Istar0 = Istar0, Dstar0 = Dstar0)
        
        # remove NA rows (inc model only)
        postPred <- postPred[!is.na(postPred[,1]),]
        
        postMeans <- rowMeans(postPred)
        postCI <- apply(postPred, 1, quantile, probs = c(0.025, 0.975))
        
        if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIHRD_noAlarm')) {
            
            postPredictFit <- data.frame(time = rep(1:tau, 3),
                                         marg = rep(c('inc', 'hosp', 'death'), each = tau),
                                         mean = postMeans,
                                         lower = postCI[1,],
                                         upper = postCI[2,])
            
        } else {
            
            postPredictFit <- data.frame(time = rep(1:tau, 1),
                                         marg = rep(c('inc'), each = tau),
                                         mean = postMeans,
                                         lower = postCI[1,],
                                         upper = postCI[2,])
        }
        
        
    } else {
        postPredictFit <- data.frame(time = NA, 
                                     marg = NA,
                                     mean = NA,
                                     lower = NA,
                                     upper = NA)
    }
    
    
    ##############################################################################
    ### Posterior predictive forecasting for models where it makes sense
    print('Posterior predictive forecasting...')

    if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIR_inc', 'SIHRD_noAlarm', 'SIR_noAlarm')) {
        postPred <- postPredForecast(incData = incData, 
                                     modelType = modelType, peak = peak,
                                     smoothC = smoothC, smoothD = smoothD, 
                                     hospData = hospData, deathData = deathData, 
                                     paramsSamples = samples,
                                     N = N, S0 = S0, I0 = I0, H0 = H0, D0 = D0, R0 = R0,
                                     Istar0 = Istar0, Dstar0 = Dstar0)
        
        # remove NA rows (inc model only)
        postPred <- postPred[!is.na(postPred[,1]),]
        
        postMeans <- rowMeans(postPred)
        postCI <- apply(postPred, 1, quantile, probs = c(0.025, 0.975))
        
        nDaysSim <- 50
        timeRange <- tau + 1:nDaysSim
        
        if (modelType %in% c('SIHRD_full', 'SIHRD_inc', 'SIHRD_noAlarm')) {
            
            postPredictForecast <- data.frame(time = rep(timeRange, 3),
                                              marg = rep(c('inc', 'hosp', 'death'), each = nDaysSim),
                                              mean = postMeans,
                                              lower = postCI[1,],
                                              upper = postCI[2,])
            
        } else {
            
            postPredictForecast <- data.frame(time = rep(timeRange, 1),
                                              marg = rep(c('inc'), each = nDaysSim),
                                              mean = postMeans,
                                              lower = postCI[1,],
                                              upper = postCI[2,])
        }
        
        
    } else {
        postPredictForecast <- data.frame(time = NA, 
                                     marg = NA,
                                     mean = NA,
                                     lower = NA,
                                     upper = NA)
    }
    
    ##############################################################################
    ### output
    list(gdiag = gdiag,
         postParams = postParams,
         postAlarmTime = postAlarmTime,
         postR0 = postR0,
         waic = waic,
         postPredictFit = postPredictFit,
         postPredictForecast = postPredictForecast)
    
}





