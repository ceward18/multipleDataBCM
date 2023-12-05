################################################################################
# Run multivariate alarm models 
# Montreal
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

peak <- c(1,2,4)
smoothWindow <- 30
modelType <- c('SIHRD_full', 'SIR_full', 'SIR_inc', 'SIHRD_noAlarm', 'SIR_noAlarm')
assumeType <- c('undetected', 'casesOnly')

allModels <- expand.grid(modelType = modelType,
                         peak = peak,
                         assumeType = assumeType,
                         stringsAsFactors = FALSE)


# 30 total - 10 per peak
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
    dat <- read.csv('./data/montrealClean.csv')
    dat$smoothedCases <- round(dat$dailyCases)
    dat$cumulativeCases <- cumsum(dat$smoothedCases)
    
    ### filter out data before peak 1
    
    # constants for all models
    N <- dat$Population[1]
    
    # smoothed incidence/hosp/deaths to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    dat$smoothC <- head(movingAverage(c(0, dat$dailyCases), smoothWindow), -1)
    dat$smoothD <- head(movingAverage(c(0, dat$dailyDeaths), smoothWindow), -1)
    
    # get data for the specified peak
    incData <- dat$smoothedCases[which(dat$peak == peak_i)]
    hospData <- round(dat$dailyHosp[which(dat$peak == peak_i)])
    deathData <- round(dat$dailyDeaths[which(dat$peak == peak_i)])
    
    smoothC <- dat$smoothC[which(dat$peak == peak_i)]
    smoothD <- dat$smoothD[which(dat$peak == peak_i)]
    
    if (peak_i == 1) {
        idxStart <- 2
        incData <- incData[-c(1:idxStart)]
        hospData <- hospData[-c(1:idxStart)]
        deathData <- deathData[-c(1:idxStart)]
        
        smoothC <- smoothC[-c(1:idxStart)]
        smoothD <- smoothD[-c(1:idxStart)]
        idxStart <- min(which(dat$peak == peak_i)) + idxStart
    } else {
        idxStart <- min(which(dat$peak == peak_i))
    }
    
    # define priors for initial conditions
    lengthI <- 5
    cumInf <- dat$cumulativeCases[max(1, (idxStart - 1))]
    cumHosp <- round(cumsum(dat$dailyHosp)[max(1, (idxStart - 1))])
    
    if (peak_i == 1) {
        S0 <- N - cumInf - cumHosp 
    } else {
        S0 <- N - cumInf 
    }
    
    I0 <- sum(dat$smoothedCases[max(1, (idxStart - lengthI)):(idxStart - 1)])
    H0 <- cumHosp
    if (peak_i == 1) {
        # for model which assumes all cases observed, H0 > I0 doesn't work
        if (assumeType_i == 'casesOnly') {
            # move H0 to I0
            I0 <- I0 + H0 
            H0 <- 0
        }
        # if we assume undetected cases, this is fine
    }
    D0 <- round(cumsum(dat$dailyDeaths)[max(1, (idxStart - 1))])
    R0 <- N - S0 - I0 - H0 - D0 # should be > 0
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
        alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
        alarmTimePost <- cbind.data.frame(postSummaries$postAlarmTime, modelInfo)
        R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
        IstarPost <- cbind.data.frame(postSummaries$postIstar, modelInfo)
        waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)
        predFitPost <- cbind.data.frame(postSummaries$postPredictFit, modelInfo)
        
    } else {
        gr <- rbind.data.frame(gr, 
                               cbind.data.frame(postSummaries$gdiag, modelInfo))
        paramsPost <- rbind.data.frame(paramsPost, 
                                       cbind.data.frame(postSummaries$postParams, modelInfo))
        alarmPost <- rbind.data.frame(alarmPost, 
                                      cbind.data.frame(postSummaries$postAlarm, modelInfo))
        alarmTimePost <- rbind.data.frame(alarmTimePost, 
                                          cbind.data.frame(postSummaries$postAlarmTime, modelInfo))
        R0Post <- rbind.data.frame(R0Post, 
                                   cbind.data.frame(postSummaries$postR0, modelInfo))
        IstarPost <- rbind.data.frame(IstarPost, 
                                      cbind.data.frame(postSummaries$postIstar, modelInfo))
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
saveRDS(alarmPost, paste0('./output/alarmPost_Batch', idxPrint, '.rds'))
saveRDS(alarmTimePost, paste0('./output/alarmTimePost_Batch', idxPrint, '.rds'))
saveRDS(R0Post, paste0('./output/R0Post_Batch', idxPrint, '.rds'))
saveRDS(IstarPost, paste0('./output/IstarPost_Batch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPost_Batch', idxPrint, '.rds'))
saveRDS(predFitPost, paste0('./output/predFitPostBatch', idxPrint, '.rds'))


