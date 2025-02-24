################################################################################
# Supplemental figures from simulation study 
################################################################################


### set up

library(ggplot2)
library(nimble)
library(grid)
library(gridExtra)
library(ggh4x)
library(knitr)
library(kableExtra)
library(openxlsx)
library(plyr)
library(dplyr)

source('./scripts/model_code.R')

################################################################################
# Fig S1 compare true infections to the estimate using conditional inference 


pdf('./figures/S1_undetectInfEst.pdf', height = 5, width = 10)
layoutMatrix <- matrix(c(1,2,3,4,4,4), byrow = T, nrow = 2)
layout(layoutMatrix, heights = c(0.8, 0.2))

par(mar = c(5.1, 4.1, 4.1, 2.1))
for (dataType in c('inc', 'death', 'equal')) {
    simData <- readRDS(paste0('./data/sim_', dataType, '.rds'))
    
    # gather incidence, hospitalizations, and deaths
    incData <- simData[, grep('^Istar', colnames(simData))]
    incDataObs <- simData[, grep('detectIstar', colnames(simData))]
    incDataEst <- round(incDataObs/0.25)
    
    plotTitle <- switch(dataType,
                        'inc' = bquote(atop(bold('High case importance'), 
                                            bold(alpha == 0.85))),
                        'death' = bquote(atop(bold('High deaths importance'), 
                                              alpha == 0.15)),
                        'equal' = bquote(atop(bold('Equal importance'), 
                                              alpha ==  0.5)))
    
    plot(incData[1,], type = 'l', col = 'black', ylim = c(0, 20000),
         main = plotTitle, cex.lab = 1.2, cex.main = 1.5, lty = 2,
         ylab = 'Count', xlab = 'Epidemic Time')
    for (i in 1:nrow(incData)) {
        lines(incData[i,], col = adjustcolor('black', alpha = 0.4), lty = 2)
        lines(incDataEst[i,], col =  adjustcolor('tomato', alpha = 0.4))
    }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend('top', c('True infections', 'Estimate including undetected'), lwd = 3,
       col = c('black', 'tomato'), 
       lty = c(2, 1, 1, 1),
       cex = 1.5, bty = 'n', xpd = T, horiz = T)
dev.off()


################################################################################
# Fig S3 R0 RMSE across time

# read in GR so the few models that didn't converge are excluded
grAll <- readRDS('./results/grAll.rds')
grAll <- grAll[-grep('comp_init', grAll$param),]

# which didn't converge
notConverge <- grAll[which(round(grAll$gr, 1) > 1.1),  ]
notConvergeModels <-  notConverge[
    !duplicated(notConverge[,-which(colnames(notConverge) 
                                    %in% c('gr', 'grUpper', 'param'))]),
    c('dataType', 'modelType', 'assumeType', 'simNumber')]
notConvergeModels$noConverge <- 1


R0PostAll <- readRDS('./results/R0PostAll.rds')

# get true average R0 for each setting
simDataInc <- readRDS('./data/sim_inc.rds')
simDataDeath <- readRDS('./data/sim_death.rds')
simDataEqual <- readRDS('./data/sim_equal.rds')

# r0
r0Data <- data.frame(time = 1:29,
                     dataType = rep(c('inc', 'death', 'equal'), each = 29),
                     truth = c(colMeans(simDataInc[, grep('R0', colnames(simDataInc))]),
                               colMeans(simDataDeath[, grep('R0', colnames(simDataDeath))]),
                               colMeans(simDataEqual[, grep('R0', colnames(simDataEqual))])))

R0PostAll <- merge(R0PostAll, r0Data, by = c('time', 'dataType'), all.x = T)

R0PostAll <- merge(R0PostAll, notConvergeModels, 
                   by = c('dataType', 'modelType', 'assumeType', 'simNumber'),
                   all.x = T)

R0PostAll <- subset(R0PostAll, is.na(noConverge))

R0PostAll$dataType <- factor(R0PostAll$dataType,
                             levels = c('inc', 'death', 'equal'),
                             labels = c('Cases important',
                                        'Deaths important',
                                        'Equal importance'))

R0PostAll$compartmentType <- ifelse(grepl('SIR', R0PostAll$modelType),
                                    'SIR', 'SIHRD')

R0PostAll$alarmType <- 'No alarm'
R0PostAll$alarmType <- ifelse(grepl('inc', R0PostAll$modelType),
                              'Cases only', R0PostAll$alarmType)

R0PostAll$alarmType <- ifelse(grepl('full', R0PostAll$modelType),
                              'Cases + deaths', R0PostAll$alarmType)

R0PostAll$assumeType <- factor(R0PostAll$assumeType,
                               levels = c('undetected', 'casesOnly'),
                               labels = c('w/ undetected', 'w/o undetected'))



R0postMSE <- R0PostAll %>%
    group_by(dataType, assumeType, compartmentType, alarmType, time) %>%
    summarise(rmse = sqrt(mean((mean - truth)^2)),
              bias = mean(mean - truth)) %>%
    data.frame()

R0postMSE$assumeType <- factor(R0postMSE$assumeType,
                               labels = c('Yes', 'No'))

strip_design <- strip_nested(
    text_x = elem_list_text(size = c(14, 12)),
    by_layer_x = TRUE
)


pdf('./figures/S2_R0Post_RMSE_time.pdf', height = 5, width = 10)
ggplot(subset(R0postMSE), 
       aes(x = time, y = rmse, col = assumeType)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    facet_nested(dataType ~ compartmentType + alarmType,
                 strip = strip_design ) +
    theme_bw() + 
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = 'Epidemic Time', y = 'RMSE', color = 'Undetected\ninfections\nmodeled',
         title =expression('RMSE of'~R[0](t))) +
    scale_color_manual(values = c('steelblue1', 'tomato'))
dev.off()

################################################################################
# Figure S3 - posterior parameter distributions for SIHRD models

paramsPostAll <- readRDS('./results/paramsPostAll.rds')

paramsTruth <- read.xlsx('simParamsSummary.xlsx')


### wide to long
paramsTruth <- reshape(paramsTruth, 
                       varying = colnames(paramsTruth)[-1], 
                       v.names = "truth",
                       timevar = "param", 
                       times = colnames(paramsTruth)[-1], 
                       new.row.names = 1:1000,
                       direction = "long")
paramsTruth <- paramsTruth[-which(colnames(paramsTruth) %in% c('id'))]

# merge with truth
paramsPostAll <- merge(paramsPostAll, paramsTruth, 
                       by = c('param', 'dataType'),
                       all.x = T)

paramsPostAll <- merge(paramsPostAll, notConvergeModels, 
                       by = c('dataType', 'modelType', 'assumeType', 'simNumber'),
                       all.x = T)

paramsPostAll <- subset(paramsPostAll, is.na(noConverge))

paramsPostAll <- paramsPostAll[order(paramsPostAll$dataType, 
                                     paramsPostAll$modelType,
                                     paramsPostAll$assumeType,
                                     paramsPostAll$simNumber, 
                                     paramsPostAll$param),]

paramsPostAll$dataType <- factor(paramsPostAll$dataType,
                                 levels = c('inc', 'death', 'equal'),
                                 labels = c('High cases importance',
                                            'High deaths importance',
                                            'Equal importance'))

paramsPostAll$param <- factor(paramsPostAll$param,
                              levels = c('alpha', 'k',
                                         'beta', 'probDetect',
                                         'nu', 'w0', 
                                         'lambda', 'gamma1', 'gamma2', 'phi'))


paramsPostAll$assumeType <- factor(paramsPostAll$assumeType,
                                   levels = c('casesOnly', 'undetected'),
                                   labels = c('w/o undetected', 'w/ undetected'))


pdf('./figures/S3_paramEsts_SIHRD_full.pdf', height = 9, width = 10)
ggplot(data = subset(paramsPostAll,
                     modelType == 'SIHRD_full' &
                         !is.na(param) ),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, linewidth = 0.8) +
    geom_point(size = 0.7) + 
    geom_errorbar(width=0, position = position_dodge(width = 0.9)) +
    facet_nested(param ~ dataType + assumeType, scales = 'free') +
    labs(x = 'Simulation Number', y = '', title = 'Posterior mean and 95% CI') +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 13, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
dev.off()



pdf('./figures/S4_paramEsts_SIR_full.pdf', height = 7, width = 10)
ggplot(data = subset(paramsPostAll,
                     modelType == 'SIR_full' &
                         !is.na(param) ),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, linewidth = 0.8) +
    geom_point(size = 0.5) + 
    geom_errorbar(width=0, position = position_dodge(width = 0.9)) +
    facet_nested(param ~ dataType + assumeType, scales = 'free') +
    labs(x = 'Simulation Number', y = '', title = 'Posterior mean and 95% CI') +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 13, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
dev.off()


# prob IR undetected vs cases only

tmp_undetected <- subset(paramsPostAll,
              modelType == 'SIR_full' & param %in% c( 'w0', 'nu') & 
                  assumeType == 'w/ undetected' & 
                  dataType == 'Equal importance') 
tmp_casesOnly <- subset(paramsPostAll,
                         modelType == 'SIR_full' & param %in% c( 'w0', 'nu')& 
                            assumeType == 'w/o undetected' & 
                            dataType == 'Equal importance') 

idd_curves <- expand.grid(simNumber = 1:50,
                          assumeType = c('w/ undetected', 'w/o undetected'),
                          day_inf = 1:15, 
                          idd_curve = NA)

for (i in 1:50) {
    
    params <- tmp_undetected[tmp_undetected$simNumber == i, ]
    
    
    idd_curves[idd_curves$simNumber == i & idd_curves$assumeType == 'w/ undetected', 
               'idd_curve'] <- logitDecay(1:15, 
                                         w0=params[params$param == 'w0', 'mean'], 
                                         nu=params[params$param == 'nu', 'mean'])
    
    params <- tmp_casesOnly[tmp_casesOnly$simNumber == i, ]
    
    
    idd_curves[idd_curves$simNumber == i & idd_curves$assumeType == 'w/o undetected', 
               'idd_curve'] <- logitDecay(1:15, 
                                          w0=params[params$param == 'w0', 'mean'], 
                                          nu=params[params$param == 'nu', 'mean'])
}


ggplot(subset(idd_curves, simNumber == 1), 
       aes(x = day_inf, y = idd_curve, group = assumeType, col = assumeType)) + 
    geom_line()


pdf('./figures/S5_paramEsts_SIHRD_inc.pdf', height = 8, width = 10)
ggplot(data = subset(paramsPostAll,
                     modelType == 'SIHRD_inc' &
                         !is.na(param) ),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, linewidth = 0.8) +
    geom_point(size = 0.5) + 
    geom_errorbar(width=0, position = position_dodge(width = 0.9)) +
    facet_nested(param ~ dataType + assumeType, scales = 'free') +
    labs(x = 'Simulation Number', y = '', title = 'Posterior mean and 95% CI') +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 13, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
dev.off()


pdf('./figures/S6_paramEsts_SIR_inc.pdf', height = 6, width = 10)
ggplot(data = subset(paramsPostAll,
                     modelType == 'SIR_inc' &
                         !is.na(param) ),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, linewidth = 0.8) +
    geom_point(size = 0.5) + 
    geom_errorbar(width=0, position = position_dodge(width = 0.9)) +
    facet_nested(param ~ dataType + assumeType, scales = 'free') +
    labs(x = 'Simulation Number', y = '', title = 'Posterior mean and 95% CI') +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 13, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
dev.off()

11
pdf('./figures/S8_paramEsts_SIR_noAlarm.pdf', height = 6, width = 10)
ggplot(data = subset(paramsPostAll,
                     modelType == 'SIR_noAlarm' &
                         !is.na(param) ),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, linewidth = 0.8) +
    geom_point(size = 0.5) + 
    geom_errorbar(width=0, position = position_dodge(width = 0.9)) +
    facet_nested(param ~ dataType + assumeType, scales = 'free') +
    labs(x = 'Simulation Number', y = '', title = 'Posterior mean and 95% CI') +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 13, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
dev.off()





################################################################################
# Table S1 - Estimation of alpha undetected vs detected


# alpha estimation


paramsPostAll <- readRDS('./results/paramsPostAll.rds')
paramsTruth <- read.xlsx('simParamsSummary.xlsx')


### wide to long
paramsTruth <- reshape(paramsTruth, 
                       varying = colnames(paramsTruth)[-1], 
                       v.names = "truth",
                       timevar = "param", 
                       times = colnames(paramsTruth)[-1], 
                       new.row.names = 1:1000,
                       direction = "long")
paramsTruth <- paramsTruth[-which(colnames(paramsTruth) %in% c('id'))]

# merge with truth
paramsPostAll <- merge(paramsPostAll, paramsTruth, 
                       by = c('param', 'dataType'),
                       all.x = T)

paramsPostAll <- paramsPostAll[paramsPostAll$param == 'alpha',]


paramsPostMSE <- paramsPostAll %>%
    group_by(dataType, assumeType, modelType) %>%
    summarise(mse = sqrt(mean((mean - truth)^2)),
              bias = mean(abs(mean - truth)),
              coverage =  mean(truth > lower & truth < upper)) %>%
    data.frame()

paramsPostMSE$assumeType <- factor(paramsPostMSE$assumeType,
                                   levels = c('undetected', 'casesOnly'),
                                   labels = c('Yes', 'No'))


paramsPostMSE$compartmentType <- ifelse(grepl('SIR', paramsPostMSE$modelType),
                                        'SIR', 'SIHRD')

paramsPostMSE$dataType <- factor(paramsPostMSE$dataType,
                                 levels = c('inc', 'death', 'equal'),
                                 labels = c('High case importance',
                                            'High deaths importance',
                                            'Equal importance'))

paramsPostMSE <- paramsPostMSE[,c('dataType', 'compartmentType', 
                                  'assumeType', 'mse')]
paramsPostMSE <- paramsPostMSE[order(paramsPostMSE$dataType,
                                     paramsPostMSE$compartmentType, 
                                     paramsPostMSE$assumeType),]


paramsPostMSE_tab <- reshape(paramsPostMSE, 
                             timevar = "assumeType",
                             idvar = c("dataType", "compartmentType"),
                             direction = "wide")

paramsPostMSE_tab <- reshape(paramsPostMSE_tab, 
                             timevar = "compartmentType",
                             idvar = c("dataType"),
                             direction = "wide")

paramsPostMSE_tab[,2:5] <- round(paramsPostMSE_tab[,2:5], 2)
paramsPostMSE_tab



