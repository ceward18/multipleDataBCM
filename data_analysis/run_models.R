################################################################################
# Run multivariate alarm models 
# NYC and Montreal
# Peaks 1, 2, 4
# assuming undetected infections exist
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

### load libraries
library(parallel)
library(nimble)

### source scripts (for movingAverage function)
source('./scripts/model_code.R')

city <- c('nyc', 'montreal')
peak <- c('1', '2', '3')
smoothWindow <- 14
modelType <- c('SIHRD_full', 'SIHRD_inc',
               'SIR_full', 'SIR_inc', 'SIHRD_noAlarm', 'SIR_noAlarm')
timePeriod <- c(4, 6, 8)

allModels <- expand.grid(city = city,
                         modelType = modelType,
                         peak = peak,
                         timePeriod = timePeriod,
                         stringsAsFactors = FALSE)

# 108 total - 18 for each of 3 peaks and 2 cities
allModels <- allModels[order(allModels$city,
                             allModels$modelType,
                             allModels$peak,
                             allModels$timePeriod),]
rownames(allModels) <- NULL

tmp <- allModels[seq(1,nrow(allModels), 3),]

batchSize <- 3
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in batchIdx) {
    
    city_i <- allModels$city[i]
    modelType_i <- allModels$modelType[i]
    peak_i <- allModels$peak[i]
    timePeriod_i <- allModels$timePeriod[i]
    
    print(paste0('Running: ', city_i, 
                 ', model: ', modelType_i, 
                 ', peak: ', peak_i,
                 ', weeks: 1-', timePeriod_i))
    
    ### read data
    dat <- read.csv(paste0('./data/', city_i, 'Clean.csv'))
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
    
    # when does peak start
    idxStart <- min(which(dat$peak == peak_i))
    
    # when does data end (based on time period assumed)
    idxEnd <- idxStart + 7*timePeriod_i - 1
    
    # get data for the specified peak
    incData <- dat$dailyCases[idxStart:idxEnd]
    hospData <- dat$dailyHosp[idxStart:idxEnd]
    deathData <- dat$dailyDeaths[idxStart:idxEnd]
    smoothC <- dat$smoothC[idxStart:idxEnd]
    smoothD <- dat$smoothD[idxStart:idxEnd]
    
    # define priors for initial conditions
    # for undetected model, getting number of previously infected and removed 
    # infected before past 7 days and not currently hospitalized or already dead
    # window of time before peak starts (time X-lengthI to time X)
    lengthI <- 7  # average infectious for 7 days
    lengthH <- 15 # average hospitalized for 15 days
    inWindowI <- max(1, (idxStart - lengthI)):(idxStart - 1)
    inWindowH <- max(1, (idxStart - lengthH)):(idxStart - 1)
    
    
    # currently hospitalized
    H0 <- sum(dat$dailyHosp[inWindowH])
    
    # already dead
    D0 <- cumsum(dat$dailyDeaths)[max(1, (idxStart - 1))]
    
    # currently infectious if infected in the past 7 days
    I0 <- sum(dat$dailyCases[inWindowI])
    
    # prior for probDetect depends on peak + city
    # For NYC: wave 1 17% detection, wave 2 55% detection, wave 3 20% detection
    # For Montreal: wave 1 25% detection, wave 2 40% detection, wave 3 20% detection
    if (city_i == 'nyc') {
        probDetectTime <- c(rep(0.17, min(which(dat$peak == 2))),
                            rep(0.55, min(which(dat$peak == 3)) - min(which(dat$peak == 2))),
                            rep(0.20, nrow(dat) - min(which(dat$peak == 3))))
        
        I0 <- round(I0 / switch(peak_i, 
                                '1' = 0.17,
                                '2' = 0.55,
                                '3' = 0.2))
    } else if (city_i == 'montreal') {
        probDetectTime <- c(rep(0.25, min(which(dat$peak == 2))),
                            rep(0.40, min(which(dat$peak == 3)) - min(which(dat$peak == 2))),
                            rep(0.20, nrow(dat) - min(which(dat$peak == 3))))
        
        I0 <- round(I0 / switch(peak_i, 
                                '1' = 0.25,
                                '2' = 0.4,
                                '3' = 0.2))
    }
    
    if (peak_i == 1) {
        # need to add a few more to I0 so it makes sense with observed hospitalization data
        I0 <- I0 + H0 + 10
        # peak 1 = no previous infectious
        R0 <- D0
    } else {
        prev_inf <- round(cumsum(dat$dailyCases / probDetectTime))[idxStart - lengthI - 1]
        # remove = infected before past 7 days and not currently hospitalized or already dead
        R0 <- prev_inf - H0 - D0
    }
    
    # currently susceptible
    S0 <- N - I0 - H0 - R0 - D0
    c(S0, I0, H0, R0, D0)
    
    # used for posterior predictive fit 
    # previously observed incidence for correct smoothing at beginning of prediction
    Istar0 <- dat$dailyCases[max(1, (idxStart - smoothWindow + 1)):(idxStart - 1)]
    Dstar0 <- round(dat$dailyDeaths[max(1, (idxStart - smoothWindow + 1)):(idxStart - 1)])
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'city_i', 'modelType_i', 'peak_i',
                           'smoothC', 'smoothD', 
                           'deathData', 'hospData',
                           'N', 'S0', 'I0', 'H0', 'D0', 'R0'))
    
    resThree <- parLapply(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/fit_models.R')
        
        # debugonce(fitAlarmModel)
        fitAlarmModel(incData = incData, city = city_i,
                      modelType = modelType_i, peak = peak_i, 
                      smoothC = smoothC, smoothD = smoothD, 
                      deathData = deathData, hospData = hospData,
                      N = N, S0 = S0, I0 = I0, H0 = H0, D0 = D0, R0 = R0,
                      seed = 3)
        
    })
    stopCluster(cl)
    
    source('./scripts/summarize_post.R')
    # debugonce(summarizePost)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   city = city_i,
                                   modelType = modelType_i, peak = peak_i, 
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
                paste0('./output/chains_', modelType_i, '_city', city_i,
                       '_peak', peak_i, '_weeks', timePeriod_i, '.rds'))
    }
    
    # save results in separate files
    modelInfo <- data.frame(city = city_i,
                            modelType = modelType_i, 
                            peak = peak_i,
                            timePeriod = timePeriod_i)
    
    if (i == batchIdx[1]) {
        
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmTimePost <- cbind.data.frame(postSummaries$postAlarmTime, modelInfo)
        R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
        waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)
        predFitPost <- cbind.data.frame(postSummaries$postPredictFit, modelInfo)
        predForecastPost <- cbind.data.frame(postSummaries$postPredictForecast, modelInfo)
        
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
        predForecastPost <- rbind.data.frame(predForecastPost, 
                                        cbind.data.frame(postSummaries$postPredictForecast, modelInfo))
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
saveRDS(predForecastPost, paste0('./output/predForecastPostBatch', idxPrint, '.rds'))
