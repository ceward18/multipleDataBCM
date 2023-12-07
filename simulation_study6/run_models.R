################################################################################
# Model fitting script
################################################################################

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
idx <- as.numeric(task_id)

# load libraries
library(nimble)
library(parallel)

# source model codes for moving average function
source('./scripts/model_code.R')

# set up grid of models to fit
nSim <- 50
dataType <- c('inc', 'death', 'equal')
modelType <- c('simple', 'full', 'inc', 'simpleNoAlarm', 'fullNoAlarm')
# SIR inc+deaths, SIHRD inc+deaths, SIR inc, SIR no alarm, SIHRD no alarm
assumeType <- c('undetected', 'casesOnly')

# 1500
allFits <- expand.grid(simNumber = 1:nSim,
                       dataType = dataType,
                       modelType = modelType,
                       assumeType = assumeType,
                       stringsAsFactors = FALSE)

# 300
tmp <- allFits[seq(1,nrow(allFits), 5),]
# rownames(tmp) <- NULL

################################################################################

# fit models in batches of 5 (90 batches total)
batchSize <- 5
batchIdx <- batchSize * (idx - 1) + 1:batchSize

for (i in batchIdx) {
    
    dataType_i <- allFits$dataType[i]
    modelType_i <- allFits$modelType[i]
    assumeType_i <- allFits$assumeType[i]
    simNumber_i <- allFits$simNumber[i]
    
    print(paste0('Data Gen: ', dataType_i,
                 ', model fit: ', modelType_i, 
                 ', assumption: ', assumeType_i, 
                 ', simulation: ', simNumber_i))
    
    # load data
    simData <- readRDS(paste0('./data/sim_', dataType_i, '.rds'))
    
    # gather incidence, hospitalizations, and deaths
    incData <- simData[simNumber_i, grep('detectIstar', colnames(simData))]
    hospData <- simData[simNumber_i, grep('fromI.*1\\]', colnames(simData))]
    deathData <- simData[simNumber_i, grep('fromH.*2\\]', colnames(simData))]
    
    # smoothed incidence/hosp/deaths to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    smoothC <- head(movingAverage(c(0, incData), 30), -1)
    smoothD <- head(movingAverage(c(0, deathData), 30), -1)
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'modelType_i', 'assumeType_i', 
                           'smoothC', 'smoothD',
                           'hospData', 'deathData'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/fit_models.R')
        
        # debugonce(fitAlarmModel), 
        
        fitAlarmModel(incData = incData, modelType = modelType_i, 
                      assumeType = assumeType_i,
                      smoothC = smoothC,  smoothD = smoothD,
                      hospData = hospData, deathData = deathData, seed = x)
    })
    stopCluster(cl)
    
    
    # true incidence (to compare to estimated)
    trueInc <- simData[simNumber_i, grep('^Istar', colnames(simData))]
    
    source('./scripts/summarize_post.R')
    # debugonce(summarizePost)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   modelType = modelType_i, assumeType = assumeType_i,
                                   smoothC = smoothC, smoothD = smoothD,
                                   hospData = hospData, deathData = deathData,
                                   trueInc = trueInc)
    
    # if the model did not converge save the chains so these can be examined later
    if (!all(postSummaries$gdiag$gr < 1.1)) {
        
        # # create thinned version
        resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 10),]
        resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 10),]
        resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 10),]
        
        saveRDS(resThree,
                paste0('./output/chains_', dataType_i, '_', modelType_i,
                       '_', assumeType_i, '_',  sprintf("%03d", simNumber_i), '.rds'))
    }
    
    
    # save results in separate files
    modelInfo <- data.frame(dataType = dataType_i,
                            modelType = modelType_i,
                            assumeType = assumeType_i,
                            simNumber = simNumber_i)
    
    # posterior summaries
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

idxPrint <- sprintf("%03d",idx)


# save output in RDS form
saveRDS(gr, paste0('./output/gr_Batch', idxPrint, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPost_Batch', idxPrint, '.rds'))
saveRDS(alarmPost, paste0('./output/alarmPost_Batch', idxPrint, '.rds'))
saveRDS(alarmTimePost, paste0('./output/alarmTimePost_Batch', idxPrint, '.rds'))
saveRDS(R0Post, paste0('./output/R0Post_Batch', idxPrint, '.rds'))
saveRDS(IstarPost, paste0('./output/IstarPost_Batch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPost_Batch', idxPrint, '.rds'))
saveRDS(predFitPost, paste0('./output/predFitPostBatch', idxPrint, '.rds'))


