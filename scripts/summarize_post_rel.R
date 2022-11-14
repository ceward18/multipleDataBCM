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
source('./scripts/get_model_inputs_rel.R')
source('./scripts/get_WAIC_rel.R')

summarizePost <- function(resThree, incData, smoothC, smoothH, smoothD, N, I0, R0) {
    
    paramSamples1 <- resThree[[1]][,-grep('alarm|R0|yC|yH|yD', colnames(resThree[[1]]))]
    paramSamples2 <- resThree[[2]][,-grep('alarm|R0|yC|yH|yD', colnames(resThree[[2]]))]
    paramSamples3 <- resThree[[3]][,-grep('alarm|R0|yC|yH|yD', colnames(resThree[[3]]))]
    
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
    ### posterior distribution of alarm function over its basis
    
    alarmSamples1 <- resThree[[1]][,grep('y', colnames(resThree[[1]]))]
    alarmSamples2 <- resThree[[2]][,grep('y', colnames(resThree[[2]]))]
    alarmSamples3 <- resThree[[3]][,grep('y', colnames(resThree[[3]]))]
    alarmSamples <- rbind(alarmSamples1, alarmSamples2, alarmSamples3)
    
    postMeans <- colMeans(alarmSamples)
    postCI <- apply(alarmSamples, 2, quantile, probs = c(0.025, 0.975))
    
    # get xAlarm
    modelInputs <- getModelInput(incData, smoothC, smoothH, smoothD, N, I0, R0)
    
    postAlarm <- data.frame(xAlarm = c(modelInputs$xC,
                                       modelInputs$xC,
                                       modelInputs$xC), 
                            marg = c(rep('inc', 10),
                                     rep('death', 10),
                                     rep('hosp', 10)),
                            mean = postMeans,
                            lower = postCI[1,],
                            upper = postCI[2,])
    
    rownames(postAlarm) <- NULL
    
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
    
    # need yAlarms on the logit scale
    yAlarmCols <- grep('y', colnames(samples))
    samples[,yAlarmCols] <- logit(samples[,yAlarmCols])
    colnames(samples)[yAlarmCols] <- paste0('logit_',  colnames(samples)[yAlarmCols])
    
    waic <- getWAIC(samples = samples, incData = incData, 
                    smoothC = smoothC, smoothH = smoothH, smoothD = smoothD, 
                    N = N, I0 = I0, R0 = R0)
    
    
    
    ### output
    list(gdiag = gdiag,
         postParams = postParams,
         postAlarm = postAlarm,
         postR0 = postR0,
         waic = waic)
    
    
    
}





