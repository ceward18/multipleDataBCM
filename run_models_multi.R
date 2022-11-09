################################################################################
# Run multivariate alarm models (BMA)
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
source('./scripts/model_codes.R')

### read data
dat <- read.csv('./data/nycClean.csv')
dat$smoothedCases <- round(dat$dailyCases)
dat$cumulativeCases <- cumsum(dat$smoothedCases)

peak <- c('1', '2', '3', '4', '5')
smoothWindow <- 60

# 5
allModels <- expand.grid(peak = peak,
                         smoothWindow = smoothWindow)


allModels <- allModels[order(allModels$smoothWindow, 
                             allModels$peak),]
rownames(allModels) <- NULL

# list of batches
tmp <- allModels[seq(1, nrow(allModels), 5),]
rownames(tmp) <- NULL

# constants for all models
N <- dat$Population[1]

# used to initialize I0 and R0
lengthI <- 3

peak_i <- allModels$peak[idx]
smoothWindow_i <- allModels$smoothWindow[idx]

print(paste('Running multi alarm, peak:', peak_i))

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
    idxStart <- 5
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
clusterExport(cl, list('incData', 'smoothC', 'smoothH', 'smoothD', 'N', 'I0', 'R0'))

resThree <- parLapplyLB(cl, 1:3, function(x) {
    
    library(nimble)
    
    # source relevant scripts
    source('./scripts/model_fit_multi.R')
    
    debugonce(fitAlarmModel)
    samples <- fitAlarmModel(incData = incData, smoothC = smoothC, smoothH = smoothH,
                             smoothD = smoothD, N = N, I0 = I0, R0 = R0, seed = x)
    
})
stopCluster(cl)


source('./scripts/summarize_post_multi.R')
# debugonce(summarizePost)
# debugonce(postPredFit)
postSummaries <- summarizePost(resThree = resThree, incData = incData,
                               smoothC = smoothC, smoothH = smoothH,
                               smoothD = smoothD, N = N, I0 = I0, R0 = R0)

# if the model did not converge save the chains so these can be examined later
if (!all(postSummaries$gdiag$gr < 1.1)) {
    
    # create thinned version
    resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 10),]
    resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 10),]
    resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 10),]
    
    saveRDS(resThree, 
            paste0('./output/chains_multi_peak', 
                   peak_i, '_', smoothWindow_i, '.rds'))
}


# save results in separate files
modelInfo <- data.frame(alarmBase = 'multi',
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



idxPrint <- sprintf("%02d",idx)

# save output in RDS form
saveRDS(gr, paste0('./output/gr_multiBatch', idxPrint, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPost_multiBatch', idxPrint, '.rds'))
saveRDS(alarmPost, paste0('./output/alarmPost_multiBatch', idxPrint, '.rds'))
saveRDS(R0Post, paste0('./output/R0Post_multiBatch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPost_multiBatch', idxPrint, '.rds'))


