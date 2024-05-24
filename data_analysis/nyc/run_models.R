################################################################################
# Run multivariate alarm models 
# NYC
# Peaks 1, 2, 4
# SIR and SIHRD models based on incidence and deaths
# SIR model based on incidence only
# assuming all cases reported and only undetected cases
################################################################################


task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

### load libraries
library(parallel)
library(nimble)

### source scripts (for movingAverage function)
source('./scripts/model_code.R')

peak <- c('1', '2', '4')
smoothWindow <- 30
modelType <- c('SIHRD_full', 'SIHRD_inc',
               'SIR_full', 'SIR_inc', 'SIHRD_noAlarm', 'SIR_noAlarm')
assumeType <- c('undetected', 'casesOnly')

allModels <- expand.grid(modelType = modelType,
                         peak = peak,
                         assumeType = assumeType,
                         stringsAsFactors = FALSE)



# 36 total - 12 per peak
allModels <- allModels[order(allModels$modelType,
                             allModels$assumeType, 
                             allModels$peak),]
rownames(allModels) <- NULL

tmp <- allModels[seq(1,nrow(allModels), 3),]

batchSize <- 3
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in batchIdx) {
    
    modelType_i <- allModels$modelType[i]
    assumeType_i <- allModels$assumeType[i]
    peak_i <- allModels$peak[i]
    
    print(paste0('Running: ', modelType_i, 
                 ', assumption: ', assumeType_i, 
                 ', peak: ', peak_i))
    
    ### read data
    dat <- read.csv(paste0('./data/nycClean.csv'))
    dat$dailyCases <- round(dat$dailyCases)
    dat$dailyHosp <- round(dat$dailyHosp)
    dat$dailyDeaths <- round(dat$dailyDeaths)
    dat$cumulativeCases <- cumsum(dat$dailyCases)
    
    # population size
    N <- dat$Population[1]
    
    # smoothed incidence/hosp/deaths to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    # (only used by alarm function models)
    dat$smoothC <- head(movingAverage(c(0, dat$dailyCases), smoothWindow), -1)
    dat$smoothD <- head(movingAverage(c(0, dat$dailyDeaths), smoothWindow), -1)
    
    # get data for the specified peak
    incData <- dat$dailyCases[which(dat$peak == peak_i)]
    hospData <- dat$dailyHosp[which(dat$peak == peak_i)]
    deathData <- dat$dailyDeaths[which(dat$peak == peak_i)]
    smoothC <- dat$smoothC[which(dat$peak == peak_i)]
    smoothD <- dat$smoothD[which(dat$peak == peak_i)]
    
    
    # when does peak start
    idxStart <- min(which(dat$peak == peak_i))
    
    # define priors for initial conditions
    # for undetected model, getting number of previously infected and removed 
    # infected before past 7 days and not currently hospitalized or already dead
    # window of time before peak starts (time X-lengthI to time X)
    lengthI <- 7  # average infectious for 7 days
    lengthH <- 15 # average hospitalized for 15 days
    inWindowI <- max(1, (idxStart - lengthI)):(idxStart - 1)
    
    # currently hospitalized
    H0 <- sum(dat$dailyHosp[max(1, (idxStart - lengthH)):(idxStart - 1)])
    
    # already dead
    D0 <- cumsum(dat$dailyDeaths)[max(1, (idxStart - 1))]
    
    # currently infectious if infected in the past 7 days
    I0 <- sum(dat$dailyCases[max(1, (idxStart - lengthI)):(idxStart - 1)])
    
    # should multiply cases, multiplier depends on peak
    #   (https://www.healthdata.org/sites/default/files/covid_briefs/101_briefing_Canada.pdf)
    # wave 1: Feb 25 - 11 July 2020          25% detection
    # wave 2: Aug 23, 2020 - March 20, 2021  40% detection
    # wave 3: July 18 - Dec 4, 2021          25% detection
    # wave 4: Dec 5, 2021 - Mar 12, 2021     20% detection
    if (assumeType_i == 'undetected') {
        probDetectTime <- c(rep(0.25, min(which(dat$peak == 2))),
                            rep(0.4, min(which(dat$peak == 3)) - min(which(dat$peak == 2))),
                            rep(0.25, min(which(dat$peak == 4)) - min(which(dat$peak == 3))),
                            rep(0.2, nrow(dat) - min(which(dat$peak == 4))))
        prev_inf <- cumsum(dat$dailyCases / probDetectTime)[idxStart - lengthI - 1]
        
        I0 <- round(I0 / switch(peak_i, 
                          '1' = 0.25,
                          '2' = 0.4,
                          '3' = 0.25, 
                          '4' = 0.2))
    } else {
        prev_inf <- cumsum(dat$dailyCases)[idxStart - lengthI - 1]
        
        if (peak_i == 1) {
            # for model which assumes all cases observed, H0 > I0 doesn't work
            # move H0 to I0
            I0 <- I0 + H0 
            H0 <- 0
            # if we assume undetected cases, this is fine
        }
    }
    if (peak_i == 1) {
        prev_inf <- 0
    } 
    
    # remove = infected before past 7 days and not currently hospitalized or already dead
    R0 <- prev_inf - H0 - D0
    
    # currently susceptible
    S0 <- N - I0 - H0 - R0 - D0
    c(S0, I0, H0, R0, D0)
    
    # used for posterior predictive fit 
    # previously observed incidence for correct smoothing at beginning of prediction
    Istar0 <- dat$smoothedCases[max(1, (idxStart - smoothWindow + 1)):(idxStart - 1)]
    Dstar0 <- round(dat$dailyDeaths[max(1, (idxStart - smoothWindow + 1)):(idxStart - 1)])
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'modelType_i', 'assumeType_i', 'peak_i',
                           'smoothC', 'smoothD', 
                           'deathData', 'hospData',
                           'N', 'S0', 'I0', 'H0', 'D0', 'R0'))
    
    resThree <- parLapply(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/fit_models.R')
        
        # debugonce(fitAlarmModel)
        fitAlarmModel(incData = incData, modelType = modelType_i,
                      assumeType = assumeType_i, peak = peak_i, 
                      smoothC = smoothC, smoothD = smoothD, 
                      deathData = deathData, hospData = hospData,
                      N = N, S0 = S0, I0 = I0, H0 = H0, D0 = D0, R0 = R0,
                      seed = x)
        
    })
    stopCluster(cl)
    
    source('./scripts/summarize_post.R')
    # debugonce(summarizePost)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   modelType = modelType_i, assumeType = assumeType_i,
                                   peak = peak_i, 
                                   smoothC = smoothC, smoothD = smoothD,
                                   hospData = hospData, deathData = deathData,
                                   N = N, S0 = S0, I0 = I0, H0 = H0, D0 = D0, R0 = R0,
                                   Istar0 = Istar0, Dstar0 = Dstar0)
    
    # if the model did not converge save the chains so these can be examined later
    if (!all(postSummaries$gdiag$gr < 1.1)) {
        
        # create thinned version to save memory
        resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 100),]
        resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 100),]
        resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 100),]
        
        saveRDS(resThree, 
                paste0('./output/chains_', modelType_i, '_', assumeType_i, 
                       '_peak', peak_i, '.rds'))
    }
    
    # save results in separate files
    modelInfo <- data.frame(modelType = modelType_i, 
                            peak = peak_i,
                            assumeType = assumeType_i)
    
    if (i == batchIdx[1]) {
        
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmTimePost <- cbind.data.frame(postSummaries$postAlarmTime, modelInfo)
        R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
        waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)
        predFitPost <- cbind.data.frame(postSummaries$postPredictFit, modelInfo)
        
    } else {
        gr <- rbind.data.frame(gr, 
                               cbind.data.frame(postSummaries$gdiag, modelInfo))
        paramsPost <- rbind.data.frame(paramsPost, 
                                       cbind.data.frame(postSummaries$postParams, modelInfo))
        alarmTimePost <- rbind.data.frame(alarmTimePost, 
                                          cbind.data.frame(postSummaries$postAlarmTime, modelInfo))
        R0Post <- rbind.data.frame(R0Post, 
                                   cbind.data.frame(postSummaries$postR0, modelInfo))
        waicPost <- rbind.data.frame(waicPost, 
                                     cbind.data.frame(postSummaries$waic, modelInfo))
        predFitPost <- rbind.data.frame(predFitPost, 
                                        cbind.data.frame(postSummaries$postPredictFit, modelInfo))
    }
    
} # end loop

idxPrint <- sprintf("%02d",idx)

# save output in RDS form
saveRDS(gr, paste0('./output/gr_Batch', idxPrint, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPost_Batch', idxPrint, '.rds'))
saveRDS(alarmTimePost, paste0('./output/alarmTimePost_Batch', idxPrint, '.rds'))
saveRDS(R0Post, paste0('./output/R0Post_Batch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPost_Batch', idxPrint, '.rds'))
saveRDS(predFitPost, paste0('./output/predFitPostBatch', idxPrint, '.rds'))

