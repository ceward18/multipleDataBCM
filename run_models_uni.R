################################################################################
# Run univariate alarm models
# Peaks 1-5
# Cases, hosps, deaths
# all based on GP alarm
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

### load libraries
library(parallel)
library(nimble)

### source scripts (for movingAverage function)
source('./scripts/modelCodes.R')

### read data
dat <- read.csv('./data/nycClean.csv')
dat$smoothedCases <- round(dat$dailyCases)
dat$cumulativeCases <- cumsum(dat$smoothedCases)

peak <- c('1', '2', '3', '4', '5')
alarmBase <- c('inc', 'hosp', 'death')
smoothWindow <- 60

# 15
allModels <- expand.grid(peak = peak,
                         alarmBase = alarmBase,
                         smoothWindow = smoothWindow)


allModels <- allModels[order(allModels$alarmBase,
                             allModels$smoothWindow, 
                             allModels$peak),]
rownames(allModels) <- NULL

# list of batches
tmp <- allModels[seq(1, nrow(allModels), 4),]
rownames(tmp) <- NULL

# constants for all models
N <- dat$Population[1]

# used to initialize I0 and R0
lengthI <- 3

# batches by alarmFit (3 batches total)
batchSize <- 5
batchIdx <- batchSize * (idx - 1) + 1:batchSize


for (i in batchIdx) {
    
    peak_i <- allModels$peak[i]
    alarmBase_i <- allModels$alarmBase[i]
    smoothWindow_i <- allModels$smoothWindow[i]
    
    print(paste('Running alarm basis:', alarmBase_i, ', peak:', peak_i))
    
    # what will inform alarm function?
    # (shifted so alarm is informed only by data up to time t-1)
    dat$alarmBase <- switch(as.character(alarmBase_i),
                            'inc' = head(movingAverage(c(0, dat$dailyCases), 
                                                       smoothWindow_i), -1),
                            'hosp' = head(movingAverage(c(0, dat$dailyHosp), 
                                                       smoothWindow_i), -1),
                            'death' = head(movingAverage(c(0, dat$dailyDeaths), 
                                                       smoothWindow_i), -1))
    
    
    # get data for the specified peak
    incData <- dat$smoothedCases[which(dat$peak == peak_i)]
    alarmBase <- dat$alarmBase[which(dat$peak == peak_i)]
    
    # initialize current number of infectious and removed individuals
    if (peak_i == 1) {
        # peak 1 shift 5 days for start
        idxStart <- 5
        incData <- incData[-c(1:idxStart)]
        alarmBase <- alarmBase[-c(1:idxStart)]
    } else {
        idxStart <- min(which(dat$peak == peak_i))
        incData <- incData[-1]
        alarmBase <- alarmBase[-1]
    }
    
    # currently infectious
    I0 <- sum(dat$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)])
    R0 <- dat$cumulativeCases[idxStart] - I0 
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'alarmBase', 'N', 'I0', 'R0'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/model_fit_uni.R')
        
        # debugonce(fitAlarmModel)
        fitAlarmModel(incData = incData, alarmBase = alarmBase,
                      N = N, I0 = I0, R0 = R0, seed = x)
        
    })
    stopCluster(cl)
    
    
    source('./scripts/summarize_post_uni.R')
    # debugonce(summarizePost)
    # debugonce(postPredFit)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   alarmBase = alarmBase, N = N, I0 = I0, R0 = R0)
    
    # if the model did not converge save the chains so these can be examined later
    if (!all(postSummaries$gdiag$gr < 1.1)) {
        
        # create thinned version
        resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 10),]
        resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 10),]
        resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 10),]
        
        saveRDS(resThree, 
                paste0('./output/chains_', alarmBase_i, '_peak', 
                       peak_i, '_', smoothWindow_i, '.rds'))
    }
    
    
    # save results in separate files
    modelInfo <- data.frame(alarmBase = alarmBase_i,
                            peak = peak_i,
                            smoothWindow = smoothWindow_i)
    
    if (i == batchIdx[1]) {
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
        R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
        waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)
        
    } else {
        gr <- rbind.data.frame(gr, 
                               cbind.data.frame(postSummaries$gdiag, modelInfo))
        paramsPost <- rbind.data.frame(paramsPost, 
                                       cbind.data.frame(postSummaries$postParams, modelInfo))
        alarmPost <- rbind.data.frame(alarmPost, 
                                      cbind.data.frame(postSummaries$postAlarm, modelInfo))
        R0Post <- rbind.data.frame(R0Post, 
                                   cbind.data.frame(postSummaries$postR0, modelInfo))
        waicPost <- rbind.data.frame(waicPost, 
                                     cbind.data.frame(postSummaries$waic, modelInfo))
    }
    
} # end loop

idxPrint <- sprintf("%02d",idx)

# save output in RDS form
saveRDS(gr, paste0('./output/grBatch', idxPrint, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPostBatch', idxPrint, '.rds'))
saveRDS(alarmPost, paste0('./output/alarmPostBatch', idxPrint, '.rds'))
saveRDS(R0Post, paste0('./output/R0PostBatch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPostBatch', idxPrint, '.rds'))


