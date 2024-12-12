################################################################################
# combine batch runs
################################################################################

outputFolder <- 'output'
resultsFolder <- 'results'

outputFiles <- list.files(paste0('./', outputFolder))

################################################################################

################################################################################
# Gelman rubin

grFiles <- outputFiles[grep('gr', outputFiles)]

grAll <- readRDS(paste0('./', outputFolder, '/', grFiles[1]))

if (grepl('longer', grFiles[1])) {
    grAll$longer <- 1
} else {
    grAll$longer <- 0
}


for (i in 2:length(grFiles)) {
    gr_i <- readRDS(paste0('./', outputFolder, '/', grFiles[i]))
    
    if (grepl('longer', grFiles[i])) {
        gr_i$longer <- 1
    } else {
        gr_i$longer <- 0
    }
    
    grAll <-rbind.data.frame(grAll, gr_i)
}

grAll <- grAll[order(grAll$city, grAll$modelType, grAll$peak, grAll$timePeriod,
                     grAll$longer, decreasing = T),]

grAll <- grAll[-which(duplicated(grAll[,c('city', 'modelType', 'peak', 
                                          'timePeriod', 'param')])),]

grAll <- grAll[order(grAll$city, grAll$modelType, grAll$peak,
                     grAll$timePeriod, grAll$param),
               -which(colnames(grAll) == 'longer')]

rownames(grAll) <- NULL

saveRDS(grAll, paste0('./', resultsFolder, '/grAll.rds'))

################################################################################
# posterior alarms - for models that estimate all alarm functions

alarmFiles <- outputFiles[grep('alarmTimePost', outputFiles)]

alarmAll <- readRDS(paste0('./', outputFolder, '/', alarmFiles[1]))

if (grepl('longer', grFiles[1])) {
    alarmAll$longer <- 1
} else {
    alarmAll$longer <- 0
}
for (i in 2:length(alarmFiles)) {
    alarm_i <- readRDS(paste0('./', outputFolder, '/', alarmFiles[i]))
    
    if (grepl('longer', grFiles[i])) {
        alarm_i$longer <- 1
    } else {
        alarm_i$longer <- 0
    }
    alarmAll <-rbind.data.frame(alarmAll, alarm_i)
}


# remove NA (noAlarm models)
alarmAll <- alarmAll[!is.na(alarmAll$time),]


alarmAll <- alarmAll[order(alarmAll$city, alarmAll$modelType, alarmAll$peak,
                           alarmAll$timePeriod, alarmAll$longer, decreasing = T),]

alarmAll <- alarmAll[-which(duplicated(alarmAll[,c('city', 'modelType', 'peak', 
                                                   'timePeriod', 'time')])),]

alarmAll <- alarmAll[order(alarmAll$city, alarmAll$modelType, alarmAll$peak,
                           alarmAll$timePeriod, alarmAll$time),
                     -which(colnames(alarmAll) == 'longer')]

rownames(alarmAll) <- NULL

saveRDS(alarmAll, paste0('./', resultsFolder, '/alarmTimePostAll.rds'))

################################################################################
# posterior parameters 

paramsPostFiles <- outputFiles[grep('paramsPost', outputFiles)]

paramsPostAll <- readRDS(paste0('./', outputFolder, '/', paramsPostFiles[1]))

if (grepl('longer', grFiles[1])) {
    paramsPostAll$longer <- 1
} else {
    paramsPostAll$longer <- 0
}

for (i in 2:length(paramsPostFiles)) {
    paramsPost_i <- readRDS(paste0('./', outputFolder, '/', paramsPostFiles[i]))
    
    if (grepl('longer', grFiles[i])) {
        paramsPost_i$longer <- 1
    } else {
        paramsPost_i$longer <- 0
    }
    paramsPostAll <-rbind.data.frame(paramsPostAll, paramsPost_i)
}


paramsPostAll <- paramsPostAll[order(paramsPostAll$city, paramsPostAll$modelType, paramsPostAll$peak,
                                     paramsPostAll$timePeriod, paramsPostAll$longer, decreasing = T),]

paramsPostAll <- paramsPostAll[-which(duplicated(paramsPostAll[,c('city', 'modelType', 'peak', 
                                                                  'timePeriod', 'param')])),]

paramsPostAll <- paramsPostAll[order(paramsPostAll$city, paramsPostAll$modelType, paramsPostAll$peak,
                                     paramsPostAll$timePeriod, paramsPostAll$param),
                               -which(colnames(paramsPostAll) == 'longer')]

rownames(paramsPostAll) <- NULL

saveRDS(paramsPostAll,  paste0('./', resultsFolder, '/paramsPostAll.rds'))

################################################################################
# posterior R0 

R0PostFiles <- outputFiles[grep('R0Post', outputFiles)]

R0PostAll <- readRDS(paste0('./', outputFolder, '/', R0PostFiles[1]))

if (grepl('longer', grFiles[1])) {
    R0PostAll$longer <- 1
} else {
    R0PostAll$longer <- 0
}

for (i in 2:length(R0PostFiles)) {
    R0Post_i <- readRDS(paste0('./', outputFolder, '/', R0PostFiles[i]))
    if (grepl('longer', grFiles[i])) {
        R0Post_i$longer <- 1
    } else {
        R0Post_i$longer <- 0
    }
    R0PostAll <-rbind.data.frame(R0PostAll, R0Post_i)
}


R0PostAll <- R0PostAll[order(R0PostAll$city, R0PostAll$modelType, R0PostAll$peak,
                             R0PostAll$timePeriod, R0PostAll$longer, decreasing = T),]

R0PostAll <- R0PostAll[-which(duplicated(R0PostAll[,c('city', 'modelType', 'peak', 
                                                      'timePeriod', 'time')])),]

R0PostAll <- R0PostAll[order(R0PostAll$city, R0PostAll$modelType, R0PostAll$peak,
                             R0PostAll$timePeriod, R0PostAll$time),
                       -which(colnames(R0PostAll) == 'longer')]

rownames(R0PostAll) <- NULL

saveRDS(R0PostAll,  paste0('./', resultsFolder, '/R0PostAll.rds'))


################################################################################
# WAIC

waicFiles <- outputFiles[grep('waicPost', outputFiles)]

waicAll <- readRDS(paste0('./', outputFolder, '/', waicFiles[1]))
if (grepl('longer', grFiles[1])) {
    waicAll$longer <- 1
} else {
    waicAll$longer <- 0
}

for (i in 2:length(waicFiles)) {
    waic_i <- readRDS(paste0('./', outputFolder, '/', waicFiles[i]))
    if (grepl('longer', grFiles[i])) {
        waic_i$longer <- 1
    } else {
        waic_i$longer <- 0
    }
    waicAll <-rbind.data.frame(waicAll, waic_i)
}


waicAll <- waicAll[order(waicAll$city, waicAll$modelType, waicAll$peak,
                         waicAll$timePeriod, waicAll$longer, decreasing = T),]

waicAll <- waicAll[-which(duplicated(waicAll[,c('city', 'modelType', 'peak', 
                                                'timePeriod')])),]

waicAll <- waicAll[order(waicAll$city, waicAll$modelType, waicAll$peak,
                         waicAll$timePeriod),
                   -which(colnames(waicAll) == 'longer')]

rownames(waicAll) <- NULL

saveRDS(waicAll,  paste0('./', resultsFolder, '/waicAll.rds'))

################################################################################
# posterior predictive fit

postPredFitFiles <- outputFiles[grep('predFitPost', outputFiles)]

postPredFitAll <- readRDS(paste0('./', outputFolder, '/', postPredFitFiles[1]))
if (grepl('longer', grFiles[1])) {
    postPredFitAll$longer <- 1
} else {
    postPredFitAll$longer <- 0
}

for (i in 2:length(postPredFitFiles)) {
    postPredFit_i <- readRDS(paste0('./', outputFolder, '/', postPredFitFiles[i]))
    if (grepl('longer', grFiles[i])) {
        postPredFit_i$longer <- 1
    } else {
        postPredFit_i$longer <- 0
    }
    postPredFitAll <-rbind.data.frame(postPredFitAll, postPredFit_i)
}

# remove NA (simple model)
postPredFitAll <- postPredFitAll[!is.na(postPredFitAll$time),]


postPredFitAll <- postPredFitAll[order(postPredFitAll$city, 
                                       postPredFitAll$modelType, 
                                       postPredFitAll$peak,
                                       postPredFitAll$timePeriod, 
                                       postPredFitAll$longer, decreasing = T),]

postPredFitAll <- postPredFitAll[-which(duplicated(postPredFitAll[,c('city', 'modelType', 'peak', 
                                                                     'timePeriod', 'time', 'marg')])),]

postPredFitAll <- postPredFitAll[order(postPredFitAll$city, postPredFitAll$modelType, 
                                       postPredFitAll$peak,
                                       postPredFitAll$timePeriod,
                                       postPredFitAll$time, postPredFitAll$marg),
                                 -which(colnames(postPredFitAll) == 'longer')]

rownames(postPredFitAll) <- NULL


saveRDS(postPredFitAll,  paste0('./', resultsFolder, '/postPredFitAll.rds'))


################################################################################
# posterior predictive forecasting

postPredFiles <- outputFiles[grep('predForecastPost', outputFiles)]

postPredAll <- readRDS(paste0('./', outputFolder, '/', postPredFiles[1]))
if (grepl('longer', grFiles[1])) {
    postPredAll$longer <- 1
} else {
    postPredAll$longer <- 0
}

for (i in 2:length(postPredFitFiles)) {
    postPred_i <- readRDS(paste0('./', outputFolder, '/', postPredFiles[i]))
    if (grepl('longer', grFiles[i])) {
        postPred_i$longer <- 1
    } else {
        postPred_i$longer <- 0
    }
    postPredAll <-rbind.data.frame(postPredAll, postPred_i)
}

# remove NA (simple model)
postPredAll <- postPredAll[!is.na(postPredAll$time),]

postPredAll <- postPredAll[order(postPredAll$city, postPredAll$modelType, postPredAll$peak,
                                 postPredAll$timePeriod, postPredAll$longer, decreasing = T),]

postPredAll <- postPredAll[-which(duplicated(postPredAll[,c('city', 'modelType', 'peak', 
                                                            'timePeriod', 'time', 'marg')])),]

postPredAll <- postPredAll[order(postPredAll$city, postPredAll$modelType, 
                                 postPredAll$peak,
                                 postPredAll$timePeriod,
                                 postPredAll$time, postPredAll$marg),
                           -which(colnames(postPredAll) == 'longer')]

rownames(postPredAll) <- NULL


saveRDS(postPredAll,  paste0('./', resultsFolder, '/postPredForecastAll.rds'))

