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

for (i in 2:length(grFiles)) {
    gr_i <- readRDS(paste0('./', outputFolder, '/', grFiles[i]))
    grAll <-rbind.data.frame(grAll, gr_i)
}

saveRDS(grAll, paste0('./', resultsFolder, '/grAll.rds'))

################################################################################
# posterior alarms - for models that estimate the alarm function separately

alarmFiles <- outputFiles[grep('alarmPost', outputFiles)]

alarmAll <- readRDS(paste0('./', outputFolder, '/', alarmFiles[1]))

for (i in 2:length(alarmFiles)) {
    alarm_i <- readRDS(paste0('./', outputFolder, '/', alarmFiles[i]))
    alarmAll <-rbind.data.frame(alarmAll, alarm_i)
}

rownames(alarmAll) <- NULL

saveRDS(alarmAll, paste0('./', resultsFolder, '/alarmPostAll.rds'))


################################################################################
# posterior alarms - for models that estimate all alarm functions

alarmFiles <- outputFiles[grep('alarmTimePost', outputFiles)]

alarmAll <- readRDS(paste0('./', outputFolder, '/', alarmFiles[1]))

for (i in 2:length(alarmFiles)) {
    alarm_i <- readRDS(paste0('./', outputFolder, '/', alarmFiles[i]))
    alarmAll <-rbind.data.frame(alarmAll, alarm_i)
}

rownames(alarmAll) <- NULL

saveRDS(alarmAll, paste0('./', resultsFolder, '/alarmTimePostAll.rds'))


################################################################################
# posterior R0 

R0Files <- outputFiles[grep('R0Post', outputFiles)]

R0All <- readRDS(paste0('./', outputFolder, '/', R0Files[1]))

for (i in 2:length(R0Files)) {
    R0_i <- readRDS(paste0('./', outputFolder, '/', R0Files[i]))
    R0All <-rbind.data.frame(R0All, R0_i)
}

rownames(R0All) <- NULL

saveRDS(R0All, paste0('./', resultsFolder, '/R0PostAll.rds'))


################################################################################
# posterior parameters 

paramsPostFiles <- outputFiles[grep('paramsPost', outputFiles)]

paramsPostAll <- readRDS(paste0('./', outputFolder, '/', paramsPostFiles[1]))

for (i in 2:length(paramsPostFiles)) {
    paramsPost_i <- readRDS(paste0('./', outputFolder, '/', paramsPostFiles[i]))
    paramsPostAll <-rbind.data.frame(paramsPostAll, paramsPost_i)
}

saveRDS(paramsPostAll,  paste0('./', resultsFolder, '/paramsPostAll.rds'))

################################################################################
# WAIC

waicFiles <- outputFiles[grep('waicPost', outputFiles)]

waicAll <- readRDS(paste0('./', outputFolder, '/', waicFiles[1]))

for (i in 2:length(waicFiles)) {
    waic_i <- readRDS(paste0('./', outputFolder, '/', waicFiles[i]))
    waicAll <-rbind.data.frame(waicAll, waic_i)
}

saveRDS(waicAll,  paste0('./', resultsFolder, '/waicAll.rds'))

################################################################################
# posterior predictive fit

postPredFitFiles <- outputFiles[grep('predFitPost', outputFiles)]

postPredFitAll <- readRDS(paste0('./', outputFolder, '/', postPredFitFiles[1]))

for (i in 2:length(postPredFitFiles)) {
    postPredFit_i <- readRDS(paste0('./', outputFolder, '/', postPredFitFiles[i]))
    postPredFitAll <-rbind.data.frame(postPredFitAll, postPredFit_i)
}

# remove NA (simple model)
postPredFitAll <- postPredFitAll[!is.na(postPredFitAll$time),]

saveRDS(postPredFitAll,  paste0('./', resultsFolder, '/postPredFitAll.rds'))
