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

summarizePost <- function(resThree, incData, modelType, alarmFit, alarmBase,
                          smoothC, smoothH, smoothD, N, I0, R0) {
    
    paramSamples1 <- resThree[[1]][,-grep('alarm|R0|yAlarm', colnames(resThree[[1]]))]
    paramSamples2 <- resThree[[2]][,-grep('alarm|R0|yAlarm', colnames(resThree[[2]]))]
    paramSamples3 <- resThree[[3]][,-grep('alarm|R0|yAlarm', colnames(resThree[[3]]))]
    
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
    
    alarmSamples1 <- t(resThree[[1]][,grep('yAlarm', colnames(resThree[[1]]))])
    alarmSamples2 <- t(resThree[[2]][,grep('yAlarm', colnames(resThree[[2]]))])
    alarmSamples3 <- t(resThree[[3]][,grep('yAlarm', colnames(resThree[[3]]))])
    alarmSamples <- cbind(alarmSamples1, alarmSamples2, alarmSamples3)
    
    postMeans <- rowMeans(alarmSamples)
    postCI <- apply(alarmSamples, 1, quantile, probs = c(0.025, 0.975))
    
    # get xAlarm (range of observed incidence/hosp/deaths or grid)
    modelInputs <- getModelInput(incData, modelType, alarmFit, alarmBase, 
                                 smoothC, smoothH, smoothD, N, I0, R0)
    
    if (modelType == 'uni') {
        
        postAlarm <- data.frame(xAlarm = modelInputs$xAlarm, 
                                mean = postMeans,
                                lower = postCI[1,],
                                upper = postCI[2,])
        
    } else if (modelType == 'multi') {
        
        postAlarm <- data.frame(xC = modelInputs$grid$xC, 
                                xH = modelInputs$grid$xH, 
                                xD = modelInputs$grid$xD, 
                                mean = postMeans,
                                lower = postCI[1,],
                                upper = postCI[2,])
    }
    
    rownames(postAlarm) <- NULL
    
    ##############################################################################
    ### posterior distribution of individual alarm functions
    
    if (modelType == 'multi' & alarmFit == 'hill') {
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
        
        # get xAlarm (range of observed incidence/hosp/deaths or grid)
        modelInputs <- getModelInput(incData, modelType, alarmFit, alarmBase, 
                                     smoothC, smoothH, smoothD, N, I0, R0)
        
        postAlarmInd <- data.frame(time = rep(1:tau, 3),
                                   marg = c(rep('inc', tau),
                                            rep('hosp', tau),
                                            rep('death', tau)),
                                   mean = postMeans,
                                   lower = postCI[1,],
                                   upper = postCI[2,])
        
        rownames(postAlarmInd) <- NULL
        
    } else {
        
        postAlarmInd <- data.frame(xAlarm = NA, 
                                   marg = NA,
                                   mean = NA,
                                   lower = NA,
                                   upper = NA)
    }
    
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
    
    if (alarmFit == 'gp') {
        # need yAlarm on the logit scale
        yAlarmCols <- grep('yAlarm', colnames(samples))
        samples[,yAlarmCols] <- logit(samples[,yAlarmCols])
        colnames(samples)[yAlarmCols] <- paste0('logit_yAlarm[', 1:length(yAlarmCols), ']')
        
    }
    
    waic <- getWAIC(samples = samples, incData = incData, 
                    modelType = modelType, alarmFit = alarmFit, 
                    alarmBase = alarmBase, 
                    smoothC = smoothC, smoothH = smoothH, smoothD = smoothD, 
                    N = N, I0 = I0, R0 = R0)
    
    
    
    ### output
    list(gdiag = gdiag,
         postParams = postParams,
         postAlarm = postAlarm,
         postAlarmInd = postAlarmInd,
         postR0 = postR0,
         waic= waic)
    
    
    
}





