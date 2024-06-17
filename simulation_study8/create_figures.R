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

source('./scripts/model_code.R')

################################################################################
# Figure 1 - example epidemics

pdf('./figures/simCurves.pdf', height = 5, width = 10)
layoutMatrix <- matrix(c(1,2,3,4,4,4), byrow = T, nrow = 2)
layout(layoutMatrix, heights = c(0.8, 0.2))

par(mar = c(5.1, 4.1, 4.1, 2.1))
for (dataType in c('inc', 'death', 'equal')) {
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
    
    plot(incData[1,], type = 'l', col = 'grey', ylim = c(0, 20000),
         main = plotTitle, cex.lab = 1.2, cex.main = 1.5, lty = 2,
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
legend('top', c('True Incidence', 'Cases', 'Hospitalizations', 'Deaths'), lwd = 3,
       col = c('grey70', 'grey30', 'goldenrod2', 'royalblue'), 
       lty = c(2, 1, 1, 1),
       cex = 1.5, bty = 'n', xpd = T, horiz = T)
dev.off()


################################################################################
# Figure 2 - alpha parameters measuring relative importance


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

paramsPostAll <- paramsPostAll[order(paramsPostAll$dataType, 
                                     paramsPostAll$modelType,
                                     paramsPostAll$assumeType,
                                     paramsPostAll$simNumber, 
                                     paramsPostAll$param),]

paramsPostAll$dataType <- factor(paramsPostAll$dataType,
                                 levels = c('inc', 'death', 'equal'),
                                 labels = c('High case importance',
                                            'High deaths importance',
                                            'Equal importance'))

paramsPostAll$param <- factor(paramsPostAll$param,
                              levels = c('alpha', 'k',
                                         'beta', 'probDetect',
                                         'nu', 'w0', 
                                         'lambda', 'gamma1', 'gamma2', 'phi'))

paramsPostAll$assumeType <- factor(paramsPostAll$assumeType,
                                   levels = c('casesOnly', 'undetected'),
                                   labels = c('Cases Only',
                                              'Undetected'))

paramsPostAll <- paramsPostAll[which(paramsPostAll$modelType %in%
                                         c('SIR_full', 'SIHRD_full')),]
paramsPostAll$modelType <- factor(paramsPostAll$modelType,
                                  levels = c('SIR_full', 'SIHRD_full'),
                                  labels = c('SIR', 'SIHRD'))


pdf('./figures/alphaPost.pdf', height = 5, width = 10)
ggplot(data = subset(paramsPostAll, 
                     param %in% 'alpha'),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper, 
           col = assumeType, group = assumeType)) +
    geom_point(size = 1, position = position_dodge(width = 0.9), alpha = 0.8) + 
    geom_errorbar(width=0, position = position_dodge(width = 0.9), alpha = 0.8) +
    geom_hline(aes(yintercept = truth), col = 'black',
               linetype = 2, linewidth = 0.5) +
    facet_nested(modelType ~  dataType) +
    labs(x = 'Simulation Number', y = '', 
         title = expression(paste('Posterior mean and 95% Credible Intervals for ',
                                  alpha)),
         col = '') +
    theme_bw() + 
    ylim(0, 1) + 
    scale_color_manual(values = c('grey20', 'red')) + 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
dev.off()

paramsPostAll$model_assume_type <- paste0(paramsPostAll$modelType, ', ',
                                          paramsPostAll$assumeType)

paramsPostAll$model_assume_type <- factor(paramsPostAll$model_assume_type,
                                          levels = c('SIR, Undetected',
                                                     'SIR, Cases Only', 
                                                     'SIHRD, Undetected', 
                                                     'SIHRD, Cases Only'),
                                          labels = c('SIR, Undetected',
                                                     'SIR, Detected', 
                                                     'SIHRD, Undetected', 
                                                     'SIHRD, Detected'))


pdf('./figures/alphaPost2.pdf', height = 3.5, width = 12)
ggplot(data = subset(paramsPostAll, 
                     param %in% 'alpha'),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper, 
           col = model_assume_type)) +
    geom_point(size = 1, position = position_dodge(0.9), alpha = 0.8) + 
    geom_errorbar(width=0, position = position_dodge(0.9), alpha = 0.8) +
    geom_hline(aes(yintercept = truth), col = 'black', linetype = 2, linewidth = 0.5) +
    facet_nested( ~ dataType, scales = 'free') +
    labs(x = 'Simulation Number', y = '', 
         title = expression(paste('Posterior mean and 95% Credible Intervals for ',
                                  alpha)),
         col = '') +
    theme_bw() + 
    ylim(0, 1) + 
    scale_color_manual(values = c('black', 'grey50', 'firebrick4', 'indianred1')) + 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 13),
          plot.title = element_text(size = 15, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
dev.off()

################################################################################
# Figure 3 - estimating R0


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

R0PostAll$dataType <- factor(R0PostAll$dataType,
                             levels = c('inc', 'death', 'equal'),
                             labels = c('Cases important',
                                        'Deaths important',
                                        'Equal importance'))

R0PostAll$modelType <- factor(R0PostAll$modelType,
                              levels = c('SIHRD_full',  
                                         'SIR_full', 
                                         'SIHRD_inc','SIR_inc', 
                                         'SIHRD_noAlarm', 'SIR_noAlarm'),
                              labels = c('SIHRD cases + deaths', 
                                         'SIR cases + deaths',
                                         'SIHRD cases only', 
                                         'SIR cases only',
                                         'SIHRD no alarm',
                                         'SIR no alarm'))

R0PostAll$assumeType <- factor(R0PostAll$assumeType,
                               levels = c('undetected', 'casesOnly'),
                               labels = c('w/ undetected', 'w/o undetected'))

pdf('./figures/R0Post.pdf', height = 6, width = 12)
ggplot(subset(R0PostAll, !modelType %in% 
                  c('SIHRD no alarm', 'SIR no alarm')), 
       aes(x = time, y = mean, 
           group = simNumber, col = assumeType)) +
    geom_line() +
    geom_line(aes(y = truth), color = 'black') +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_nested(dataType ~ modelType + assumeType) +
    theme_bw() + 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 16, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    guides(color = 'none') + 
    labs(x = 'Epidemic Time', y = expression(R[0](t))) +
    scale_color_manual(values = c('steelblue1', 'darksalmon'))
dev.off()


################################################################################
# Posterior predictive fit

postPredFitAll <- readRDS('./results/postPredFitAll.rds')
postPredFitAll$marg[postPredFitAll$marg == 'inc' &
                        postPredFitAll$assumeType == 'casesOnly'] <- 'cases'

postPredFitSim <- list()

for (i in 1:8) {
    
    for (dataType in c('inc', 'death', 'equal')) {
        simData <- readRDS(paste0('./data/sim_', dataType, '.rds'))
        
        incData <- simData[i, grep('^Istar', colnames(simData))]
        caseData <- simData[i, grep('detect', colnames(simData))]
        hospData <- simData[i, grep('fromI.*1\\]', colnames(simData))]
        deathData <- simData[i, grep('fromH.*2\\]', colnames(simData))]
        
        tau <- length(incData)
        
        
        trueEpi <- data.frame(time = 1:tau,
                              truth = c(incData, caseData, hospData, deathData),
                              marg = rep(c('inc', 'cases', 'hosp', 'death'), each = tau),
                              dataType = dataType)
        
        assign(paste0('trueEpidemic', tools::toTitleCase(dataType)), trueEpi)
    }
    
    trueEpidemic <- rbind.data.frame(trueEpidemicInc, trueEpidemicDeath, trueEpidemicEqual)
    
    postPredFitSim[[i]] <- merge(subset(postPredFitAll, simNumber == i), trueEpidemic, 
                                 by = c('time', 'marg', 'dataType'))
    
}

# collapse list
postPredFitSimFinal <- do.call('rbind.data.frame', postPredFitSim)


postPredFitSimFinal$dataType <- factor(postPredFitSimFinal$dataType,
                                       levels = c('inc', 'death', 'equal'),
                                       labels = c('Cases importance',
                                                  'Deaths importance',
                                                  'Equal importance'))

postPredFitSimFinal$marg <- factor(postPredFitSimFinal$marg, 
                                   levels = c('inc', 'cases', 'hosp', 'death'),
                                   labels = c('True Incidence', 'Observed Cases',
                                              'Hospitalizations',
                                              'Deaths'))

postPredFitSimFinal$modelType <- factor(postPredFitSimFinal$modelType,
                                        levels = c('SIHRD_full', 'SIHRD_inc',
                                                   'SIR_inc', 
                                                   'SIHRD_noAlarm', 
                                                   'SIR_noAlarm'),
                                        labels = c('SIHRD cases + deaths',
                                                   'SIHRD cases only',
                                                   'SIR cases only',
                                                   'SIHRD no alarm',
                                                   'SIR no alarm'))

postPredFitSimFinal <- postPredFitSimFinal[order(postPredFitSimFinal$dataType,
                                                 postPredFitSimFinal$modelType,
                                                 postPredFitSimFinal$assumeType,
                                                 postPredFitSimFinal$simNumber,
                                                 postPredFitSimFinal$marg,
                                                 postPredFitSimFinal$time),]


postPredFitSimFinal$assumeType <- factor(postPredFitSimFinal$assumeType,
                                         levels = c('undetected', 'casesOnly'),
                                         labels = c('w/ undetected', 'w/o undetected'))

postPredFitSimFinal <- postPredFitSimFinal[postPredFitSimFinal$marg != 'True Incidence',]

p1 <- ggplot(subset(postPredFitSimFinal, simNumber == 1 & 
                        modelType %in% c('SIHRD cases + deaths', 
                                         'SIHRD cases only', 
                                         'SIR cases only')), 
             aes(x = time, ymin = lower, ymax = upper, col = marg, fill = marg)) + 
    geom_line(linewidth = 0.5, aes(y = truth)) +
    geom_line(aes(y = mean), linetype = 2, linewidth = 0.5) + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(dataType~modelType + assumeType,  scales = 'free_y') +
    scale_color_manual(values = c('black', 'goldenrod2', 'royalblue')) + 
    scale_fill_manual(values = c('black', 'goldenrod2', 'royalblue')) +
    theme_bw() + 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(y = 'Epidemic Time', x = 'Count') + 
    guides(fill = 'none', col = 'none')

# combine two plots
p2 <- ggplot(subset(postPredFitSimFinal, simNumber == 1 & 
                        modelType %in% c('SIHRD no alarm', 'SIR no alarm')), 
             aes(x = time, ymin = lower, ymax = upper, col = marg, fill = marg)) + 
    geom_line(linewidth = 0.5, aes(y = truth)) +
    geom_line(aes(y = mean), linetype = 2, linewidth = 0.5) + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(dataType~modelType + assumeType,  scales = 'free_y') +
    scale_color_manual(values = c('black', 'goldenrod2', 'royalblue')) + 
    scale_fill_manual(values = c('black', 'goldenrod2', 'royalblue')) +
    theme_bw() + 
    ylim(0, 100) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(y = 'Epidemic Time', x = 'Count',
         col = '', fill = '')

pdf('./figures/postPred.pdf', height = 5, width = 10)
p1
dev.off()

# pdf('./figures/postPred.pdf', height = 5, width = 14)
# grid.arrange(p1, p2, nrow = 1, widths = c(0.45, 0.55),
#              top=textGrob("Posterior predictive distribution", 
#                           gp = gpar(fontsize = 15, font = 1)))
# dev.off()


# all on one plot and only observed cases

ggplot(subset(postPredFitSimFinal, marg == 'Observed Cases' & simNumber == 1), 
       aes(x = time, ymin = lower, ymax = upper, col = modelType, fill = modelType)) + 
    geom_line(linewidth = 1, aes(y = truth), col = 'black') +
    geom_line(aes(y = mean), linetype = 2, linewidth = 1) + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(simNumber~dataType + assumeType, independent = 'y', scales = 'free') +
    # scale_color_manual(values = c('grey30', 'dodgerblue')) + 
    # scale_fill_manual(values = c('grey30', 'dodgerblue')) + 
    theme_bw()
