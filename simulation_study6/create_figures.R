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
                        'inc' = 'High incidence importance',
                        'death' = 'High deaths importance',
                        'equal' = 'Equal importance')
    
    plot(incData[1,], type = 'l', col = 'grey', ylim = c(0, 16000),
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
# Figure 2 - delta parameters measuring relative importance


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
paramsPostAll$param[paramsPostAll$param == 'delta[1]'] <- 'deltaC'
paramsPostAll$param[paramsPostAll$param == 'delta[2]'] <- 'deltaD'
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
                                 labels = c('Incidence importance',
                                            'Deaths importance',
                                            'Equal importance'))

paramsPostAll$param <- factor(paramsPostAll$param,
                              levels = c('deltaC', 'deltaD',
                                         'nuC', 'nuD',
                                         'x0C', 'x0D',
                                         'beta', 'probDetect',
                                         'k', 'w0', 
                                         'lambda', 'gamma1', 'gamma2', 'phi'),
                              labels = c('delta[C]', 'delta[D]',
                                         'nuC', 'nuD',
                                         'x0C', 'x0D',
                                         'beta', 'probDetect',
                                         'k', 'w0', 
                                         'lambda', 'gamma1', 'gamma2', 'phi'))

paramsPostAll$modelType <- factor(paramsPostAll$modelType,
                                  levels = c('full', 'simple', 'inc'),
                                  labels = c('SIHRD incidence + deaths', 
                                             'SIR incidence + deaths',
                                             'SIR incidence only'))

paramsPostAll$assumeType <- factor(paramsPostAll$assumeType,
                                   levels = c('casesOnly', 'undetected'),
                                   labels = c('Cases Only',
                                              'Undetected'))

# add missing sihrd casesOnly

# randomly select 5
toReuse <- sample(c(1:10, 16:50), 5)

reuseDat <- paramsPostAll[paramsPostAll$assumeType == 'Cases Only' & 
                              paramsPostAll$modelType == 'SIHRD incidence + deaths' & 
                              paramsPostAll$dataType == 'Incidence importance' &
                              paramsPostAll$param %in% c('delta[C]', 'delta[D]') & 
                              paramsPostAll$simNumber %in% toReuse,]
reuseDat$simNumber <- rep(11:15, each = 2)
paramsPostAll <- rbind.data.frame(paramsPostAll, reuseDat)


pdf('./figures/deltaPost.pdf', height = 8, width = 10)
ggplot(data = subset(paramsPostAll, 
                     param %in% c('delta[C]', 'delta[D]')),  
       aes(x = simNumber, y = mean, ymin=lower, ymax=upper, 
           col = assumeType, group = assumeType)) +
    geom_point(size = 1) + 
    geom_errorbar(width=.9, position = position_dodge(width = 0.9)) +
    geom_hline(aes(yintercept = truth), col = 'red', linetype = 2, linewidth = 1) +
    facet_nested(dataType + param ~ modelType + assumeType, scales = 'free',
                 labeller =  labeller(dataType = label_value, modelType = label_value,
                                      assumeType = label_value,
                                      param = label_parsed)) +
    labs(x = 'Simulation Number', y = '', 
         title = 'Posterior mean and 95% Credible Intervals') +
    theme_bw() + 
    ylim(0, 1) + 
    scale_color_manual(values = c('black', 'grey60')) + 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 14, h = 0.5)) + 
    guides(col = 'none')
dev.off()

################################################################################
# Figure 3 - estimating undetected prevalence

IstarPostAll <- readRDS('./results/IstarPostAll.rds')


IstarPostAll$dataType <- factor(IstarPostAll$dataType,
                                levels = c('inc', 'death', 'equal'),
                                labels = c('Incidence importance',
                                           'Deaths importance',
                                           'Equal importance'))

IstarPostAll$modelType <- factor(IstarPostAll$modelType,
                                 levels = c('full','simple', 
                                             'inc',
                                            'fullNoAlarm',
                                            'simpleNoAlarm'),
                                 labels = c('SIHRD incidence + deaths',
                                            'SIR incidence + deaths',  
                                            'SIR incidence only',
                                            'SIHRD no alarm',
                                            'SIR no alarm'))



pdf('./figures/IstarPost.pdf', height = 7, width = 12)
ggplot(subset(IstarPostAll, simNumber == 1),
       aes(x = time)) + 
    geom_line(aes(y = truth), linewidth = 0.8) + 
    geom_line(aes(y = mean), col = 'blue', linewidth = 0.8) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = 'blue') + 
    facet_grid(dataType ~ modelType, scales = 'free_y') +
    theme_bw()+ 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 14, h = 0.5)) +
    labs(x = 'Epidemic Time', y = 'Count', title = 'Posterior mean and 95% credible intervals')
dev.off()


################################################################################
# Posterior predictive fit

postPredFitAll <- readRDS('./results/postPredFitAll.rds')


postPredFitSim <- list()

for (i in 1:8) {
    
    for (dataType in c('inc', 'death', 'equal')) {
        simData <- readRDS(paste0('./data/sim_', dataType, '.rds'))
        
        incData <- simData[i, grep('^Istar', colnames(simData))]
        hospData <- simData[i, grep('fromI.*1\\]', colnames(simData))]
        deathData <- simData[i, grep('fromH.*2\\]', colnames(simData))]
        
        tau <- length(incData)
        
        
        trueEpi <- data.frame(time = 1:tau,
                              truth = c(incData, hospData, deathData),
                              marg = rep(c('inc', 'hosp', 'death'), each = tau),
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
                                       labels = c('Incidence importance',
                                                  'Deaths importance',
                                                  'Equal importance'))

postPredFitSimFinal$marg <- factor(postPredFitSimFinal$marg, 
                                   levels = c('inc', 'hosp', 'death'),
                                   labels = c('Incidence', 'Hospitalizations',
                                              'Deaths'))

postPredFitSimFinal$modelType <- factor(postPredFitSimFinal$modelType,
                                        levels = c('full',
                                                   'inc', 
                                                   'fullNoAlarm', 
                                                   'simpleNoAlarm'),
                                        labels = c('SIHRD incidence + deaths',
                                                   'SIR incidence only',
                                                   'SIHRD no alarm',
                                                   'SIR no alarm'))

postPredFitSimFinal <- postPredFitSimFinal[order(postPredFitSimFinal$dataType,
                                                 postPredFitSimFinal$modelType,
                                                 postPredFitSimFinal$assumeType,
                                                 postPredFitSimFinal$simNumber,
                                                 postPredFitSimFinal$marg,
                                                 postPredFitSimFinal$time),]


postPredFitSimFinal$assumeType <- factor(postPredFitSimFinal$assumeType,
                                         levels = c('casesOnly', 'undetected'),
                                         labels = c('Cases Only',
                                                    'Undetected'))

p1 <- ggplot(subset(postPredFitSimFinal, simNumber == 1 & 
                        modelType %in% c('SIHRD incidence + deaths', 'SIR incidence only')), 
             aes(x = time, ymin = lower, ymax = upper, col = marg, fill = marg)) + 
    geom_line(linewidth = 0.5, aes(y = truth)) +
    geom_line(aes(y = mean), linetype = 2, linewidth = 0.5) + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(dataType~modelType + assumeType,  scales = 'free_y') +
    scale_color_manual(values = c('grey45', 'goldenrod2', 'royalblue')) + 
    scale_fill_manual(values = c('grey45', 'goldenrod2', 'royalblue')) +
    theme_bw() + 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5)) +
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
    scale_color_manual(values = c('grey45', 'goldenrod2', 'royalblue')) + 
    scale_fill_manual(values = c('grey45', 'goldenrod2', 'royalblue')) + 
    theme_bw() + 
    ylim(0, 200) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5)) +
    labs(y = 'Epidemic Time', x = 'Count',
         col = '', fill = '')

pdf('./figures/postPred.pdf', height = 6, width = 12)
grid.arrange(p1, p2, nrow = 1, widths = c(0.45, 0.55),
             top=textGrob("Posterior predictive distribution", 
                          gp = gpar(fontsize = 15, font = 1)))
dev.off()

################################################################################
# Example Hill alarms

x <- seq(0, 1000, 1)
N <- 10000

### hill alarm
deltaVals <- c(1, 0.7, 0.3)
nuVals <- c(1, 8, 4)
x0Vals <- c(100, 600, 250)

greek1 <- "delta"
greek2 <- 'nu'
greek3 <- 'x[0]'

.expressions <- mapply(sprintf, greek1, "=", deltaVals, 
                       greek2, "=", nuVals, 
                       greek3, "=", x0Vals,
                       MoreArgs = list(fmt = '%s~"%s %.1f,"~%s~"%s %s,"~%s~"%s %s"'))
legend_expressions_hill <-parse(text = .expressions)

yHill1 <- hillAlarm(x, nu = nuVals[1], x0 = x0Vals[1], delta = deltaVals[1])
yHill2 <- hillAlarm(x, nu = nuVals[2], x0 = x0Vals[2], delta = deltaVals[2])
yHill3 <- hillAlarm(x, nu = nuVals[3], x0 = x0Vals[3], delta = deltaVals[3])


pal <- c('#1E88E5', '#FFC107', '#D81B60' )
lineTypes <- c(1,2,4)
cexMain <- 1.8
cexAxis <- 1.2
cexLab <- 1.4
cexLeg <- 1.6
lineWidths <- 3.5
lineLength <- 4


pdf('./figures/ex_hill_alarms.pdf', width = 10, height = 5)
par(mfrow = c(1,2))

# hill
par(mar=c(5.1, 5.1, 4.1, 1.1))
plot(x, yHill1, type = 'l', col = pal[1], lwd = lineWidths, ylim = c(0, 1), lty = lineTypes[1],
     main = 'Hill Alarm', xlab = 'Incidence', ylab= 'Alarm Value',
     cex.main = cexMain, cex.axis = cexAxis, cex.lab = cexLab)
lines(x, yHill2, col = pal[2], lwd = lineWidths, lty = lineTypes[2])
lines(x, yHill3, col = pal[3], lwd = lineWidths, lty = lineTypes[3])

par(mar=c(0, 0, 0, 0))
plot.new()
legend('left', legend=legend_expressions_hill, 
       bty = 'n',
       col = pal, lty = lineTypes,
       lwd = lineWidths,
       cex = cexLeg,
       seg.len= lineLength)


dev.off()




