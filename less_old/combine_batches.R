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
if('xAlarm' %in% colnames(alarmAll)) {
    alarmAll$xC <- NA
    alarmAll$xH <- NA
    alarmAll$xD <- NA
} else {
    alarmAll$xAlarm <- NA
}

for (i in 2:length(alarmFiles)) {
    alarm_i <- readRDS(paste0('./', outputFolder, '/', alarmFiles[i]))
    if('xAlarm' %in% colnames(alarm_i)) {
        alarm_i$xC <- NA
        alarm_i$xH <- NA
        alarm_i$xD <- NA
    } else {
        alarm_i$xAlarm <- NA
    }
    alarmAll <-rbind.data.frame(alarmAll, alarm_i)
}

# remove models that didn't estimate an alarm (betat and basic)
alarmAll <- alarmAll[!is.na(alarmAll$mean),]
rownames(alarmAll) <- NULL

saveRDS(alarmAll, paste0('./', resultsFolder, '/alarmPostAll.rds'))


################################################################################
# posterior alarms - for models that estimate all alarm functions

alarmFiles <- outputFiles[grep('alarmIndPost', outputFiles)]

alarmAll <- readRDS(paste0('./', outputFolder, '/', alarmFiles[1]))
if('xAlarm' %in% colnames(alarmAll)) {
    colnames(alarmAll)[which(colnames(alarmAll) == 'xAlarm')] <- 'time'
}

for (i in 2:length(alarmFiles)) {
    alarm_i <- readRDS(paste0('./', outputFolder, '/', alarmFiles[i]))
    if('xAlarm' %in% colnames(alarm_i)) {
        colnames(alarm_i)[which(colnames(alarm_i) == 'xAlarm')] <- 'time'
    }
    alarmAll <-rbind.data.frame(alarmAll, alarm_i)
}

# remove models that didn't estimate an alarm (betat and basic)
alarmAll <- alarmAll[!is.na(alarmAll$mean),]
rownames(alarmAll) <- NULL

saveRDS(alarmAll, paste0('./', resultsFolder, '/alarmIndPostAll.rds'))

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
# R0 posterior over time


r0PostFiles <- outputFiles[grep('R0Post', outputFiles)]

r0PostAll <- readRDS(paste0('./', outputFolder, '/', r0PostFiles[1]))

for (i in 2:length(r0PostFiles)) {
    r0Post_i <- readRDS(paste0('./', outputFolder, '/', r0PostFiles[i]))
    r0PostAll <-rbind.data.frame(r0PostAll, r0Post_i)
}

rownames(r0PostAll) <- NULL

saveRDS(r0PostAll,  paste0('./', resultsFolder, '/r0PostAll.rds'))
