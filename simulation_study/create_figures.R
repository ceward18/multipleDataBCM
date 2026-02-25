################################################################################
# Figures of results output 
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
# Figure 1 - example epidemics

pdf('./figures/fig1_simCurves.pdf', height = 4, width = 10)
layoutMatrix <- matrix(c(1,2,3,4,4,4), byrow = T, nrow = 2)
layout(layoutMatrix, heights = c(0.8, 0.2))

par(mar = c(5.1, 4.1, 4.1, 2.1))
for (dataType in c('inc', 'equal', 'death')) {
    simData <- readRDS(paste0('./data/sim_', dataType, '.rds'))
    
    # gather incidence, hospitalizations, and deaths
    incData <- simData[, grep('^Istar', colnames(simData))]
    incDataObs <- simData[, grep('detectIstar', colnames(simData))]
    hospData <- simData[, grep('fromI.*1\\]', colnames(simData))]
    deathData <- simData[, grep('fromH.*2\\]', colnames(simData))]
    
    plotTitle <- switch(dataType,
                        'inc' = bquote(atop(bold('High case importance'), 
                                            bold(alpha == 0.85))),
                        'death' = bquote(atop(bold('High deaths importance'), 
                                              alpha == 0.15)),
                        'equal' = bquote(atop(bold('Equal importance'), 
                                              alpha ==  0.5)))
    
    plot(incData[1,], type = 'l', col = 'grey', ylim = c(0, 1000),
         main = plotTitle, cex.lab = 1.3, cex.main = 1.7, lty = 2,
         ylab = 'Count', xlab = 'Epidemic Time')
    for (i in 1:nrow(incData)) {
        lines(incData[i,], col = adjustcolor('grey70', alpha = 0.4), lty = 2)
        lines(incDataObs[i,], col =  adjustcolor('grey30', alpha = 0.6))
        lines(hospData[i,], col =  adjustcolor('goldenrod2', alpha = 0.6))
        lines(deathData[i,], col =  adjustcolor('royalblue', alpha = 0.6))
    }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend('top', c('Infections', 'Cases', 'Hospitalizations', 'Deaths'), lwd = 3,
       col = c('grey70', 'grey30', 'goldenrod2', 'royalblue'), 
       lty = c(2, 1, 1, 1),
       cex = 1.6, bty = 'n', xpd = T, horiz = T)
dev.off()

################################################################################
# Results start here

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

################################################################################
# Fig 2 - RMSE undetected vs detected in estimation of R0 

R0PostAll <- readRDS('./results/R0PostAll.rds')
R0PostAll <- R0PostAll[-which(R0PostAll$time == 25),]

# get true average R0 for each setting
simDataInc <- readRDS('./data/sim_inc.rds')
simDataDeath <- readRDS('./data/sim_death.rds')
simDataEqual <- readRDS('./data/sim_equal.rds')

# r0
r0Data <- data.frame(time = 1:max(R0PostAll$time),
                     dataType = rep(c('inc', 'death', 'equal'), each = max(R0PostAll$time)),
                     truth = c(colMeans(simDataInc[, grep('R0', colnames(simDataInc))]),
                               colMeans(simDataDeath[, grep('R0', colnames(simDataDeath))]),
                               colMeans(simDataEqual[, grep('R0', colnames(simDataEqual))])))

R0PostAll <- merge(R0PostAll, r0Data, by = c('time', 'dataType'), all.x = T)

R0PostAll <- merge(R0PostAll, notConvergeModels, 
                   by = c('dataType', 'modelType', 'assumeType', 'simNumber'),
                   all.x = T)

R0PostAll <- subset(R0PostAll, is.na(noConverge))

R0PostAll$dataType <- factor(R0PostAll$dataType,
                             levels = c('inc',  'equal', 'death'),
                             labels = c('High case importance',
                                        'Equal importance',
                                        'High deaths importance'))

R0PostAll$compartmentType <- ifelse(grepl('SIR', R0PostAll$modelType),
                                    'SIR', 'SIHRD')

R0PostAll$alarmType <- 'No alarm'
R0PostAll$alarmType <- ifelse(grepl('inc', R0PostAll$modelType),
                              'Cases only', R0PostAll$alarmType)

R0PostAll$alarmType <- ifelse(grepl('full', R0PostAll$modelType),
                              'Cases + deaths', R0PostAll$alarmType)


R0postMSE <- R0PostAll %>%
    group_by(dataType, assumeType, compartmentType, alarmType, time) %>%
    summarise(rmse = sqrt(mean((mean - truth)^2)),
              bias = mean(mean - truth)) %>%
    data.frame()

R0postMSE$assumeType <- factor(R0postMSE$assumeType,
                               levels = c('undetected', 'casesOnly'),
                               labels = c('Modeled', 
                                          'Ignored'))

# use only time 1 and do barplots
R0postMSE_1 <- R0postMSE[R0postMSE$time %in% c(1, 24),]


R0postMSE_1$time <- factor(R0postMSE_1$time, 
                           labels = c('Start of epidemic', 
                                      'End of epidemic'))




pdf('./figures/fig2_R0Post_RMSE.pdf', height = 5, width = 10)
ggplot(R0postMSE_1, 
       aes(x = assumeType,  y = rmse, fill = alarmType)) +
    geom_bar(stat = "identity", position=position_dodge(),
             col = 'black', linewidth = 0) +
    facet_nested( time  ~  dataType + compartmentType, 
                  scales = 'free_y') +
    theme_bw() + 
    # coord_cartesian(ylim =c(0, 0.5)) +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 10),
          axis.text.y = element_text(size = 10, color = 'black'),
          axis.text.x = element_text(size = 9, color = 'black'),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          plot.title = element_text(size = 14, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = 'Incorporation of undetected infections', y = 'RMSE', fill = 'Fitted model',
         title =expression('RMSE of'~R[0](t))) +
    scale_fill_manual(values = c('steelblue1', 
                                 'tomato',
                                 'goldenrod1'))
dev.off()



################################################################################
# Figure 3 - alpha parameters measuring relative importance

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

# ## remove those that did not converge
# paramsPostAll <- merge(paramsPostAll, notConvergeModels,
#                        by = c('dataType', 'modelType', 'assumeType', 'simNumber'),
#                        all.x = T)
# paramsPostAll$noConverge[is.na(paramsPostAll$noConverge)] <- 0
# paramsPostAll <- paramsPostAll[paramsPostAll$noConverge == 0,]

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
                                 levels = c('inc', 'equal', 'death'),
                                 labels = c('High case importance',
                                            'Equal importance',
                                            'High deaths importance'))

paramsPostAll$param <- factor(paramsPostAll$param,
                              levels = c('alpha', 'k',
                                         'beta', 'probDetect',
                                         'nu', 'w0', 
                                         'lambda', 'gamma1', 'gamma2', 'phi'))


paramsPostAll <- paramsPostAll[which(paramsPostAll$modelType %in%
                                         c('SIR_full', 'SIHRD_full')),]
paramsPostAll$modelType <- factor(paramsPostAll$modelType,
                                  levels = c('SIR_full', 'SIHRD_full'),
                                  labels = c('SIR', 'SIHRD'))



pdf('./figures/fig3_alphaPost.pdf', height = 5, width = 9)
ggplot(data = subset(paramsPostAll, 
                     param %in% 'alpha' & assumeType == 'undetected'),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper)) +
    geom_point(size = 0.8, position = position_dodge(width = 0.9), alpha = 0.8) + 
    geom_errorbar(width=0, position = position_dodge(width = 0.9), alpha = 0.8,
                  linewidth = 0.2) +
    geom_hline(aes(yintercept = truth), col = 'red',
               linetype = 2, linewidth = 0.5) +
    facet_nested(modelType ~  dataType ) +
    labs(x = 'Simulation Number', y = expression(alpha), 
         title = expression(paste('Posterior mean and 95% Credible Intervals for ',
                                  alpha)),
         col = '') +
    theme_bw() + 
    ylim(0, 1) + 
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
dev.off()

################################################################################
# Figure 4 - posterior predictive fit for one simulation

postPredFitAll <- readRDS('results/postPredFitAll.rds')

# one randomly selected simulation of those where all models converged
notConvergeSims <- unique(notConvergeModels$simNumber)

convergeSims <- which(!1:50 %in% notConvergeSims)

set.seed(1234)
sim_idx <- sample(convergeSims, 1)

for (dataType in c('inc', 'death', 'equal')) {
    simData <- readRDS(paste0('./data/sim_', dataType, '.rds'))
    
    incData <- simData[sim_idx, grep('^Istar', colnames(simData))]
    caseData <- simData[sim_idx, grep('detect', colnames(simData))]
    hospData <- simData[sim_idx, grep('fromI.*1\\]', colnames(simData))]
    deathData <- simData[sim_idx, grep('fromH.*2\\]', colnames(simData))]
    
    tau <- length(incData)
    
    
    trueEpi <- data.frame(time = 1:tau,
                          truth = c(incData, caseData, hospData, deathData),
                          marg = rep(c('inc', 'cases', 'hosp', 'death'), each = tau),
                          dataType = dataType)
    
    assign(paste0('trueEpidemic', tools::toTitleCase(dataType)), trueEpi)
}

trueEpidemic <- rbind.data.frame(trueEpidemicInc, trueEpidemicDeath, trueEpidemicEqual)

postPredFitSimFinal <- merge(postPredFitAll, trueEpidemic, 
                             by = c('time', 'marg', 'dataType'))


postPredFitSimFinal$dataType <- factor(postPredFitSimFinal$dataType,
                                       levels = c('inc', 'equal', 'death'),
                                       labels = c('Cases importance',
                                                  'Equal importance',
                                                  'Deaths importance'))

postPredFitSimFinal$marg <- factor(postPredFitSimFinal$marg, 
                                   levels = c('inc', 'cases', 'hosp', 'death'),
                                   labels = c('True Incidence', 'Observed Cases',
                                              'Hospitalizations',
                                              'Deaths'))

postPredFitSimFinal$compartmentType <- ifelse(grepl('SIR', postPredFitSimFinal$modelType),
                                              'SIR', 'SIHRD')

postPredFitSimFinal$alarmType <- 'No alarm'
postPredFitSimFinal$alarmType <- ifelse(grepl('inc', postPredFitSimFinal$modelType),
                                        'Cases only', postPredFitSimFinal$alarmType)

postPredFitSimFinal$alarmType <- ifelse(grepl('full', postPredFitSimFinal$modelType),
                                        'Cases + deaths', postPredFitSimFinal$alarmType)


postPredFitSimFinal <- postPredFitSimFinal[postPredFitSimFinal$marg != 'True Incidence',]


pdf('./figures/fig4_postPred.pdf', height = 5, width = 8.5)
ggplot(subset(postPredFitSimFinal, simNumber == sim_idx &  
                  assumeType == 'undetected'), 
       aes(x = time, ymin = lower, ymax = upper, fill = marg)) + 
    geom_line(linewidth = 0.5, aes(y = truth, col = marg)) +
    geom_line(aes(y = mean, col = marg), linetype = 2, linewidth = 0.5) + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(dataType~compartmentType + alarmType,  scales = 'free_y') +
    scale_color_manual(values = c('black', 'goldenrod3', 'royalblue')) + 
    scale_fill_manual(values = c('black', 'goldenrod2', 'royalblue')) +
    theme_bw() + 
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 11),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(y = 'Epidemic Time', x = 'Count') + 
    guides(fill = 'none', col = 'none')
dev.off()









