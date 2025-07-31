################################################################################
# Run multivariate alarm models 
# Montreal
# Peaks 1, 2, 4
# ignoring undetected infections 
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

### load libraries
library(parallel)
library(nimble)

### source scripts (for movingAverage function)
source('./scripts/model_code.R')

city <- c('montreal')
peak <- c('1', '2')
smoothWindow <- 30
modelType <- c('SIHRD_full', 'SIHRD_inc',
               'SIR_full', 'SIR_inc', 
               'SIHRD_noAlarm', 'SIR_noAlarm')
timePeriod <- c(6, 8, 10)

allModels <- expand.grid(city = city,
                         modelType = modelType,
                         peak = peak,
                         timePeriod = timePeriod,
                         stringsAsFactors = FALSE)

# 36 total - 18 for each of 2 peaks
allModels <- allModels[order(allModels$city,
                             allModels$modelType,
                             allModels$peak,
                             allModels$timePeriod),]
rownames(allModels) <- NULL

# prior for probDetect depends on wave 
#   (https://www.healthdata.org/sites/default/files/covid_briefs/101_briefing_Canada.pdf)
#   (Figure 8.1)
# wave 1: Feb 25 - 11 July 2020          25% detection
# wave 2: Aug 23, 2020 - March 20, 2021  40% detection

allModels$probDetectMean <- ifelse(allModels$peak == 1, 0.25, 0.40)

tmp <- allModels[seq(1,nrow(allModels), 3),]

batchSize <- 3
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in batchIdx) {
    
    city_i <- allModels$city[i]
    modelType_i <- allModels$modelType[i]
    peak_i <- allModels$peak[i]
    timePeriod_i <- allModels$timePeriod[i]
    probDetectMean_i <- allModels$probDetectMean[i]
    
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
    #  inflate based on reporting
    I0 <- round(I0 / probDetectMean_i)
    
    
    if (peak_i == 1) {
        # need to add a few more to I0 so it makes sense with observed hospitalization data
        I0 <- I0 + 10
        
        # peak 1 = no previous infectious
        R0 <- D0
    } else {
        
        # account for waning immunity by initializing previous infections to have occurred only in the past 6 months
        lengthR <- 30*6
        inWindowR <- max(1, (idxStart - lengthI - lengthR)):(idxStart - lengthI - 1)
        
        prev_inf <- sum(dat$dailyCases[inWindowR])
        R0 <- prev_inf - H0 - D0
    }
    
    # currently susceptible
    S0 <- N - I0 - H0 - R0 - D0
    c(S0, I0, H0, R0, D0)
    
    # used for posterior predictive fit 
    # previously observed cases for correct smoothing at beginning of prediction
    Istar0 <- dat$dailyCases[max(1, (idxStart - smoothWindow + 1)):(idxStart - 1)]
    Dstar0 <- round(dat$dailyDeaths[max(1, (idxStart - smoothWindow + 1)):(idxStart - 1)])
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'city_i', 'modelType_i', 
                           'peak_i', 'probDetectMean_i',
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
                      probDetectMean = probDetectMean_i,
                      smoothC = smoothC, smoothD = smoothD, 
                      deathData = deathData, hospData = hospData,
                      N = N, S0 = S0, I0 = I0, H0 = H0, D0 = D0, R0 = R0,
                      seed = x)
        
    })
    stopCluster(cl)
    
    source('./scripts/summarize_post.R')
    # debugonce(summarizePost)
    # debugonce(postPredForecast)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   modelType = modelType_i, peak = peak_i, 
                                   probDetectMean = probDetectMean_i,
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
                paste0('./output/chains_', modelType_i, '_', city_i,
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
        IstarPost <- cbind.data.frame(postSummaries$postIstar, modelInfo)
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
        IstarPost <- rbind.data.frame(IstarPost, 
                                      cbind.data.frame(postSummaries$postIstar, modelInfo))
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
saveRDS(IstarPost, paste0('./output/IstarPost_Batch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPost_Batch', idxPrint, '.rds'))
saveRDS(predFitPost, paste0('./output/predFitPostBatch', idxPrint, '.rds'))
saveRDS(predForecastPost, paste0('./output/predForecastPostBatch', idxPrint, '.rds'))
