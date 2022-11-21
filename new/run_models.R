################################################################################
# Run multivariate alarm models 
# London and NYC
# Peaks 1-5
# Cases, hosps, deaths
# Univariate based on GP alarm
# Multivariate based on hill alarms
################################################################################


task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

### load libraries
library(parallel)
library(nimble)

### source scripts (for movingAverage function)
source('./scripts/model_codes.R')

cities <- c('nyc', 'london', 'montreal')
peak <- c('1', '2', '3', '4', '5')
smoothWindow <- c(1, 14, 30, 60)

# univariate models (20 per city)
allModelsUni <- expand.grid(city = cities,
                            modelType = 'uni',
                            alarmFit = c('gp'),
                            alarmBase = c('inc', 'hosp', 'death'),
                            smoothWindow = smoothWindow,
                            peak = peak,
                            stringsAsFactors = FALSE)

# multivaraite models (40 per city)
allModelsMulti <- expand.grid(city = cities,
                              modelType = 'multi',
                              alarmFit = c('hill', 'spline'),
                              alarmBase = 'all',
                              smoothWindow = smoothWindow,
                              peak = peak,
                              stringsAsFactors = FALSE)

allModels <- rbind.data.frame(allModelsUni,
                              allModelsMulti)

# 300 total
allModels <- allModels[order(allModels$alarmFit,
                             allModels$city, 
                             allModels$modelType, 
                             allModels$alarmBase, 
                             allModels$smoothWindow, 
                             allModels$peak),]
rownames(allModels) <- NULL


# used to initialize I0 and R0
lengthI <- 3

# batches by model (5 in each batch)
# spline models should run in their own batch
if (idx <= 48) {
    
    batchSize <- 5
    batchIdx <- batchSize * (idx - 1) + 1:batchSize
    
} else {
    # 49-108 (20 per peak = 60 total)
    
    batchSize <- 1
    batchIdx <- (length(cities) * 80) - (length(cities) * 16) + 
        batchSize * (idx - 1) + 1:batchSize
    
}


for (i in batchIdx) {
    
    city_i <- allModels$city[i]
    modelType_i <- allModels$modelType[i]
    alarmFit_i <- allModels$alarmFit[i]
    alarmBase_i <- allModels$alarmBase[i]
    smoothWindow_i <- allModels$smoothWindow[i]
    peak_i <- allModels$peak[i]
    
    print(paste0('Running: ', city_i, ', ',  modelType_i, ', ', alarmFit_i, 
                 ', ', alarmBase_i, ', ', smoothWindow_i, ', peak: ', peak_i))
    
    ### read data
    dat <- read.csv(paste0('./data/', city_i, 'Clean.csv'))
    dat$smoothedCases <- round(dat$dailyCases)
    dat$cumulativeCases <- cumsum(dat$smoothedCases)
    
    # constants for all models
    N <- dat$Population[1]
    
    # smoothed incidence/hosp/deaths to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    dat$smoothC <- head(movingAverage(c(0, dat$dailyCases), smoothWindow_i), -1)
    dat$smoothH <- head(movingAverage(c(0, dat$dailyHosp), smoothWindow_i), -1)
    dat$smoothD <- head(movingAverage(c(0, dat$dailyDeaths), smoothWindow_i), -1)
    
    # get data for the specified peak
    incData <- dat$smoothedCases[which(dat$peak == peak_i)]
    smoothC <- dat$smoothC[which(dat$peak == peak_i)]
    smoothH <- dat$smoothH[which(dat$peak == peak_i)]
    smoothD <- dat$smoothD[which(dat$peak == peak_i)]
    
    if (peak_i == 1) {
        idxStart <- min(which(dat$peak == peak_i)) + 5
        incData <- incData[-c(1:idxStart)]
        smoothC <- smoothC[-c(1:idxStart)]
        smoothH <- smoothH[-c(1:idxStart)]
        smoothD <- smoothD[-c(1:idxStart)]
    } else {
        idxStart <- min(which(dat$peak == peak_i))
        incData <- incData[-1]
        smoothC <- smoothC[-1]
        smoothH <- smoothH[-1]
        smoothD <- smoothD[-1]
    }
    
    # currently infectious
    I0 <- sum(dat$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)])
    R0 <- dat$cumulativeCases[idxStart] - I0 
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'modelType_i', 'alarmFit_i', 'alarmBase_i',
                           'smoothC', 'smoothH', 'smoothD', 'N', 'I0', 'R0'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/model_fit.R')
        
        # debugonce(fitAlarmModel)
        fitAlarmModel(incData = incData, modelType = modelType_i,
                      alarmFit = alarmFit_i, alarmBase = alarmBase_i,
                      smoothC = smoothC, smoothH = smoothH,
                      smoothD = smoothD, N = N, I0 = I0, R0 = R0, seed = x)
        
        
        
    })
    stopCluster(cl)
    
    
    source('./scripts/summarize_post.R')
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   modelType = modelType_i, alarmFit = alarmFit_i, 
                                   alarmBase = alarmBase_i,
                                   smoothC = smoothC, smoothH = smoothH,
                                   smoothD = smoothD, N = N, I0 = I0, R0 = R0)
    
    # if the model did not converge save the chains so these can be examined later
    if (!all(postSummaries$gdiag$gr < 1.1)) {
        
        # create thinned version
        resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 100),]
        resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 100),]
        resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 100),]
        
        saveRDS(resThree, 
                paste0('./output/chains_', city_i, '_', modelType_i, 
                       '_', alarmFit_i, '_', alarmBase_i,
                       '_peak', peak_i, '_', smoothWindow_i, '.rds'))
    }
    
    # save results in separate files
    modelInfo <- data.frame(city = city_i,
                            modelType = modelType_i, 
                            alarmFit = alarmFit_i, 
                            alarmBase = alarmBase_i,
                            peak = peak_i,
                            smoothWindow = smoothWindow_i)
    
    if (i == batchIdx[1]) {
        
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
        alarmIndPost <- cbind.data.frame(postSummaries$postAlarmInd, modelInfo)
        R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
        waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)
        
    } else {
        gr <- rbind.data.frame(gr, 
                               cbind.data.frame(postSummaries$gdiag, modelInfo))
        paramsPost <- rbind.data.frame(paramsPost, 
                                       cbind.data.frame(postSummaries$postParams, modelInfo))
        alarmPost <- rbind.data.frame(alarmPost, 
                                      cbind.data.frame(postSummaries$postAlarm, modelInfo))
        alarmIndPost <- rbind.data.frame(alarmIndPost, 
                                         cbind.data.frame(postSummaries$postAlarmInd, modelInfo))
        R0Post <- rbind.data.frame(R0Post, 
                                   cbind.data.frame(postSummaries$postR0, modelInfo))
        waicPost <- rbind.data.frame(waicPost, 
                                     cbind.data.frame(postSummaries$waic, modelInfo))
    }
    
} # end loop



idxPrint <- sprintf("%02d",idx)

# save output in RDS form
saveRDS(gr, paste0('./output/gr_Batch', idxPrint, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPost_Batch', idxPrint, '.rds'))
saveRDS(alarmPost, paste0('./output/alarmPost_Batch', idxPrint, '.rds'))
saveRDS(alarmIndPost, paste0('./output/alarmIndPost_Batch', idxPrint, '.rds'))
saveRDS(R0Post, paste0('./output/R0Post_Batch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPost_Batch', idxPrint, '.rds'))


