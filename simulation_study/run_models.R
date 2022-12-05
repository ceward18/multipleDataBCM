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
dataType <- c('inc', 'hosp', 'death', 'equal')
modelType <- c('simple', 'full')

# 400
allFits <- expand.grid(simNumber = 1:nSim,
                       dataType = dataType,
                       modelType = modelType,
                       stringsAsFactors = FALSE)


tmp <- allFits[seq(1,nrow(allFits), 10),]
rownames(tmp) <- NULL

################################################################################

# fit models in batches of 25 (60 batches total)
batchSize <- 10
batchIdx <- batchSize * (idx - 1) + 1:batchSize

for (i in batchIdx) {
    
    dataType_i <- allFits$dataType[i]
    modelType_i <- allFits$modelType[i]
    simNumber_i <- allFits$simNumber[i]
    
    print(paste0('Data Gen: ', dataType_i,
                 ', model fit: ', modelType_i, 
                 ', simulation: ', simNumber_i))
    
    # load data
    simData <- readRDS(paste0('./data/sim_', dataType_i, '.rds'))
    
    # gather incidence, hospitalizations, and deaths
    incData <- simData[simNumber_i, grep('Istar', colnames(simData))]
    hospData <- simData[simNumber_i, grep('fromI.*1\\]', colnames(simData))]
    deathData <- simData[simNumber_i, grep('fromH.*2\\]', colnames(simData))]
    
    # smoothed incidence/hosp/deaths to inform alarm function 
    # (shifted so alarm is informed only by data up to time t-1)
    smoothC <- head(movingAverage(c(0, incData), 30), -1)
    smoothH <- head(movingAverage(c(0, hospData), 30), -1)
    smoothD <- head(movingAverage(c(0, deathData), 30), -1)
    
    # run three chains in parallel
    cl <- makeCluster(3)
    clusterExport(cl, list('incData', 'modelType_i', 'smoothC', 'smoothH', 'smoothD',
                           'hospData', 'deathData'))
    
    resThree <- parLapplyLB(cl, 1:3, function(x) {
        
        library(nimble)
        
        # source relevant scripts
        source('./scripts/fit_models.R')
        
        fitAlarmModel(incData = incData, modelType = modelType_i, 
                      smoothC = smoothC,  smoothH = smoothH, smoothD = smoothD,
                      hospData = hospData, deathData = deathData, seed = x)
    })
    stopCluster(cl)

    source('./scripts/summarize_post.R')
    
    # debugonce(summarizePost)
    postSummaries <- summarizePost(resThree = resThree, incData = incData,
                                   modelType = modelType_i, 
                                   smoothC = smoothC, smoothH = smoothH, smoothD = smoothD,
                                   hospData = hospData, deathData = deathData)
    
    # if the model did not converge save the chains so these can be examined later
    if (!all(postSummaries$gdiag$gr < 1.1)) {
        
        # # create thinned version
        resThree[[1]] <- resThree[[1]][seq(1,nrow(resThree[[1]]), 10),]
        resThree[[2]] <- resThree[[2]][seq(1,nrow(resThree[[2]]), 10),]
        resThree[[3]] <- resThree[[3]][seq(1,nrow(resThree[[3]]), 10),]
        # 
        # saveRDS(resThree, 
        #         paste0('./output/chains_', dataType_i, '_', modelType_i,
        #                '_', simNumber_i, '.rds'))
        
        print(paste0('./output/chains_', dataType_i, '_', modelType_i,
                                    '_', simNumber_i, '.rds'))
    }
    
    
    # save results in separate files
    modelInfo <- data.frame(dataType = dataType_i,
                            modelType = modelType_i,
                            simNumber = simNumber_i)
    
    # posterior summaries
    if (i == batchIdx[1]) {
        
        gr <- cbind.data.frame(postSummaries$gdiag, modelInfo)
        paramsPost <- cbind.data.frame(postSummaries$postParams, modelInfo)
        alarmPost <- cbind.data.frame(postSummaries$postAlarm, modelInfo)
        alarmTimePost <- cbind.data.frame(postSummaries$postAlarmTime, modelInfo)
        R0Post <- cbind.data.frame(postSummaries$postR0, modelInfo)
        waicPost <- cbind.data.frame(postSummaries$waic, modelInfo)
        
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
        waicPost <- rbind.data.frame(waicPost, 
                                     cbind.data.frame(postSummaries$waic, modelInfo))
    }
    
} # end loop

idxPrint <- sprintf("%02d",idx)


# save output in RDS form
saveRDS(gr, paste0('./output/gr_Batch', idxPrint, '.rds'))
saveRDS(paramsPost, paste0('./output/paramsPost_Batch', idxPrint, '.rds'))
saveRDS(alarmPost, paste0('./output/alarmPost_Batch', idxPrint, '.rds'))
saveRDS(alarmTimePost, paste0('./output/alarmTimePost_Batch', idxPrint, '.rds'))
saveRDS(R0Post, paste0('./output/R0Post_Batch', idxPrint, '.rds'))
saveRDS(waicPost, paste0('./output/waicPost_Batch', idxPrint, '.rds'))


