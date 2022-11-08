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
source('./scripts/getModelInputs.R')
source('./scripts/getWAIC.R')

summarizePost <- function(resThree, incData, smoothC, smoothH, smoothD,
                          N, I0, R0, Rstar0, lengthI, 
                          alarmBase, alarmFit, infPeriod) {
    
    paramSamples1 <- resThree[[1]][,-grep('alarm|yC|yH|yD|fCases|fHosp|fDeath|R0', colnames(resThree[[1]]))]
    paramSamples2 <- resThree[[2]][,-grep('alarm|yC|yH|yD|fCases|fHosp|fDeath|R0', colnames(resThree[[2]]))]
    paramSamples3 <- resThree[[3]][,-grep('alarm|yC|yH|yD|fCases|fHosp|fDeath|R0', colnames(resThree[[3]]))]
    
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
    ### posterior distribution of alarm function over epidemic time 
    
    alarmSamples1 <- t(resThree[[1]][,grep('alarm', colnames(resThree[[1]]))])
    alarmSamples2 <- t(resThree[[2]][,grep('alarm', colnames(resThree[[2]]))])
    alarmSamples3 <- t(resThree[[3]][,grep('alarm', colnames(resThree[[3]]))])
    alarmSamples <- cbind(alarmSamples1, alarmSamples2, alarmSamples3)
    
    postMeans <- rowMeans(alarmSamples)
    postCI <- apply(alarmSamples, 1, quantile, probs = c(0.025, 0.975))
    
    postAlarm <- data.frame(time = 1:length(postMeans), 
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
    postR0 <- data.frame(time = 1:length(incData),
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postR0) <- NULL
    
    ##############################################################################
    ### Posterior distribution of marginal alarm contributions over observed range
    
    getXRanges <- getModelInput(alarmFit = alarmFit, incData = incData,
                                smoothC = smoothC, smoothH = smoothH,
                                smoothD = smoothD, infPeriod = infPeriod, 
                                N = N, I0 = I0, R0 = R0, Rstar0 = Rstar0, 
                                lengthI = lengthI)
    
    # cases (alpha[1])
    yCSamples1 <- resThree[[1]][,grep('yC', colnames(resThree[[1]]))]
    yCSamples1 <- yCSamples1 * paramSamples1[,'alpha[1]']
    yCSamples2 <- resThree[[2]][,grep('yC', colnames(resThree[[2]]))]
    yCSamples2 <- yCSamples2 * paramSamples2[,'alpha[1]']
    yCSamples3 <- resThree[[3]][,grep('yC', colnames(resThree[[3]]))]
    yCSamples3 <- yCSamples3 * paramSamples3[,'alpha[1]']
    yCSamples <- rbind(yCSamples1, yCSamples2, yCSamples3)
    
    postMeans <- colMeans(yCSamples)
    postCI <- apply(yCSamples, 2, quantile, probs = c(0.025, 0.975))
    
    postyC <- data.frame(dataSource = 'case',
                         xAlarm = getXRanges$xC,
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postyC) <- NULL
    
    # hospitalizations (alpha[2])
    yHSamples1 <- resThree[[1]][,grep('yH', colnames(resThree[[1]]))]
    yHSamples1 <- yHSamples1 * paramSamples1[,'alpha[2]']
    yHSamples2 <- resThree[[2]][,grep('yH', colnames(resThree[[2]]))]
    yHSamples2 <- yHSamples2 * paramSamples2[,'alpha[2]']
    yHSamples3 <- resThree[[3]][,grep('yH', colnames(resThree[[3]]))]
    yHSamples3 <- yHSamples3 * paramSamples3[,'alpha[2]']
    yHSamples <- rbind(yHSamples1, yHSamples2, yHSamples3)
    
    postMeans <- colMeans(yHSamples)
    postCI <- apply(yHSamples, 2, quantile, probs = c(0.025, 0.975))
    
    postyH <- data.frame(dataSource = 'hosp',
                         xAlarm = getXRanges$xH,
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postyH) <- NULL
    
    # deaths (alpha[3])
    yDSamples1 <- resThree[[1]][,grep('yD', colnames(resThree[[1]]))]
    yDSamples1 <- yDSamples1 * paramSamples1[,'alpha[3]']
    yDSamples2 <- resThree[[2]][,grep('yD', colnames(resThree[[2]]))]
    yDSamples2 <- yDSamples2 * paramSamples2[,'alpha[3]']
    yDSamples3 <- resThree[[3]][,grep('yD', colnames(resThree[[3]]))]
    yDSamples3 <- yDSamples3 * paramSamples3[,'alpha[3]']
    yDSamples <- rbind(yDSamples1, yDSamples2, yDSamples3)
    
    postMeans <- colMeans(yDSamples)
    postCI <- apply(yDSamples, 2, quantile, probs = c(0.025, 0.975))
    
    postyD <- data.frame(dataSource = 'death',
                         xAlarm = getXRanges$xD,
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postyD) <- NULL
    
    postYAlarms <- rbind.data.frame(postyC, postyH, postyD)
    
    ##############################################################################
    ### Posterior distribution of marginal alarm contributions over epidemic time
    
    # cases (alpha[1])
    yCSamples1 <- resThree[[1]][,grep('fCases', colnames(resThree[[1]]))]
    yCSamples1 <- yCSamples1 * paramSamples1[,'alpha[1]']
    yCSamples2 <- resThree[[2]][,grep('fCases', colnames(resThree[[2]]))]
    yCSamples2 <- yCSamples2 * paramSamples2[,'alpha[1]']
    yCSamples3 <- resThree[[3]][,grep('fCases', colnames(resThree[[3]]))]
    yCSamples3 <- yCSamples3 * paramSamples3[,'alpha[1]']
    yCSamples <- rbind(yCSamples1, yCSamples2, yCSamples3)
    
    postMeans <- colMeans(yCSamples)
    postCI <- apply(yCSamples, 2, quantile, probs = c(0.025, 0.975))
    
    postyC <- data.frame(dataSource = 'case',
                         time = 1:length(incData),
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postyC) <- NULL
    
    # hospitalizations (alpha[2])
    yHSamples1 <- resThree[[1]][,grep('fHosp', colnames(resThree[[1]]))]
    yHSamples1 <- yHSamples1 * paramSamples1[,'alpha[2]']
    yHSamples2 <- resThree[[2]][,grep('fHosp', colnames(resThree[[2]]))]
    yHSamples2 <- yHSamples2 * paramSamples2[,'alpha[2]']
    yHSamples3 <- resThree[[3]][,grep('fHosp', colnames(resThree[[3]]))]
    yHSamples3 <- yHSamples3 * paramSamples3[,'alpha[2]']
    yHSamples <- rbind(yHSamples1, yHSamples2, yHSamples3)
    
    postMeans <- colMeans(yHSamples)
    postCI <- apply(yHSamples, 2, quantile, probs = c(0.025, 0.975))
    
    postyH <- data.frame(dataSource = 'hosp',
                         time = 1:length(incData),
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postyH) <- NULL
    
    # deaths (alpha[3])
    yDSamples1 <- resThree[[1]][,grep('fDeath', colnames(resThree[[1]]))]
    yDSamples1 <- yDSamples1 * paramSamples1[,'alpha[3]']
    yDSamples2 <- resThree[[2]][,grep('fDeath', colnames(resThree[[2]]))]
    yDSamples2 <- yDSamples2 * paramSamples2[,'alpha[3]']
    yDSamples3 <- resThree[[3]][,grep('fDeath', colnames(resThree[[3]]))]
    yDSamples3 <- yDSamples3 * paramSamples3[,'alpha[3]']
    yDSamples <- rbind(yDSamples1, yDSamples2, yDSamples3)
    
    postMeans <- colMeans(yDSamples)
    postCI <- apply(yDSamples, 2, quantile, probs = c(0.025, 0.975))
    
    postyD <- data.frame(dataSource = 'death',
                         time = 1:length(incData),
                         mean = postMeans,
                         lower = postCI[1,],
                         upper = postCI[2,])
    rownames(postyD) <- NULL
    
    postFAlarms <- rbind.data.frame(postyC, postyH, postyD)
    
    ##############################################################################
    ### WAIC values
    
    # samples to use for WAIC calculation differ by model
    samples <- rbind(resThree[[1]], resThree[[2]], resThree[[3]])
    
    if (alarmFit == 'gp') {
        # need yAlarm on the logit scale
        yAlarmCols <- grep('y', colnames(samples))
        samples[,yAlarmCols] <- logit(samples[,yAlarmCols])
        colnames(samples)[yAlarmCols] <- paste0('logit_',  colnames(samples)[yAlarmCols])
    }
    
    waic <- getWAIC(samples = samples, incData = incData, 
                    smoothC = smoothC, smoothH = smoothH, smoothD = smoothD, 
                    N = N, I0 = I0, R0 = R0, Rstar0 = Rstar0, lengthI = lengthI,
                    infPeriod = 'fixed',  alarmFit = alarmFit)
    
  
    
    ### output
    list(gdiag = gdiag,
         postParams = postParams,
         postAlarm = postAlarm,
         postR0 = postR0,
         postYAlarms = postYAlarms,
         postFAlarms = postFAlarms,
         waic= waic)
    
    
    
}





