################################################################################
# Figures and tables in paper
# Data Analysis of NYC and Montreal
################################################################################


# load libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(nimble)
library(plyr)
library(dplyr)
library(knitr)
library(kableExtra)
library(ggh4x)

################################################################################
# Plot Montreal data


### set up Montreal data
montreal <- read.csv('data/montrealClean.csv')
montreal$date <- as.Date(montreal$date)
montreal$Peak <- factor(montreal$peak)

# just 2020
montreal <- montreal[montreal$date >= min(montreal$date[which(montreal$peak == 1)]),]

montreal <- montreal[montreal$date <= as.Date('2020-12-01'),]

pal <- c('red', 'blue', 'goldenrod1')

peak1Dates <- montreal$date[min(which(montreal$peak == 1))]
peak1Dates <- c(peak1Dates, peak1Dates + 6*7)

peak2Dates <- montreal$date[min(which(montreal$peak == 2))]
peak2Dates <- c(peak2Dates, peak2Dates + 6*7)


pdf('figures/montreal_data.pdf', height = 3.5, width = 7)
ggplot(montreal, aes(x = date, y = dailyCases)) + 
    geom_rect(aes(xmin = peak1Dates[1], 
                  xmax = peak1Dates[2], ymin = -Inf, ymax = Inf), 
              fill = 'grey90', alpha = 0.2) +
    geom_rect(aes(xmin = peak2Dates[1], 
                  xmax = peak2Dates[2], ymin = -Inf, ymax = Inf), 
              fill = 'grey90', alpha = 0.2) +
    geom_line() + 
    geom_line(aes(y = dailyDeaths), col = 'red')+
    scale_y_continuous(labels = scales::comma, limits = c(0, 550)) +
    scale_x_date(date_breaks = "3 month", date_minor_breaks = "1 month",
                 date_labels = "%b '%y") +
    labs(x = 'Date', y = 'Incidence', title = 'Montreal COVID-19 Case and Death Counts') +
    annotate('text', x = as.Date('2020-03-28'), y = 550, label = 'Wave 1') +
    annotate('text', x = as.Date('2020-09-13'), y = 550, label = 'Wave 2') +
    theme_bw() + 
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          plot.title = element_text(size = 11, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
dev.off()

################################################################################
# posterior alpha parameters by wave

paramsPostAll <- readRDS('results/paramsPostAll.rds')

alphaPost <- paramsPostAll[paramsPostAll$param == 'alpha' & 
                               paramsPostAll$city == 'montreal' & 
                               paramsPostAll$peak %in% 1:2 & 
                               paramsPostAll$timePeriod == 6,]

alphaPost$compartmentModel <- ifelse(alphaPost$modelType == 'SIHRD_full', 
                                     'SIHRD', 'SIR')

pdf('figures/montreal_alpha.pdf', height = 3.5, width = 7)
ggplot(data = alphaPost,  
       aes(x = peak, y = mean, ymin=lower, ymax=upper, 
           col = compartmentModel)) +
    geom_point(size = 2, position = position_dodge(width = 0.9), alpha = 0.8) + 
    geom_errorbar(width=0.3, linewidth = 1, position = position_dodge(width = 0.9), alpha = 0.8) +
    labs(x = 'Wave', y = '', 
         title = expression(paste('Posterior mean and 95% Credible Intervals for ',
                                  alpha)),
         col = '') +
    theme_bw() + 
    ylim(0, 1) + 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    scale_color_manual(values = c('dodgerblue', 'tomato'))
dev.off()


################################################################################
# Posterior predictive fit + alarms

postPredFitAll <- readRDS('./results/postPredFitAll.rds')

postPredFitAll <- postPredFitAll[postPredFitAll$city == 'montreal' & 
                                     postPredFitAll$peak %in% 1:2 & 
                                     postPredFitAll$timePeriod == 6,]


montreal_peak_lengths <- c(sum(montreal$peak == 1, na.rm = T),
                           sum(montreal$peak == 2, na.rm = T))
montreal_dat <- data.frame(city = 'montreal',
                           time = c(1:montreal_peak_lengths[1], 
                                    1:montreal_peak_lengths[2]), 
                           peak = rep(c(1,2), montreal_peak_lengths),
                           date = c(montreal$date[which(montreal$peak == 1)],
                                    montreal$date[which(montreal$peak == 2)]),
                           inc = c(montreal$dailyCases[which(montreal$peak == 1)],
                                   montreal$dailyCases[which(montreal$peak == 2)]),
                           hosp = c(montreal$dailyHosp[which(montreal$peak == 1)],
                                    montreal$dailyHosp[which(montreal$peak == 2)]),
                           death = c(montreal$dailyDeaths[which(montreal$peak == 1)],
                                     montreal$dailyDeaths[which(montreal$peak == 2)]))

# need in long format to merge with posterior predictions

montreal_dat_long <- reshape(montreal_dat, 
                             varying = c("inc", "hosp", "death"), 
                             v.names = "truth",
                             timevar = "marg", 
                             times = c("inc", "hosp", "death"), 
                             new.row.names = 1:10000,
                             direction = "long")

montreal_dat_long <- montreal_dat_long[,-which(colnames(montreal_dat_long) == 'id')]

postPredFitAll <- merge(postPredFitAll, montreal_dat_long, 
                        by = c('city', 'peak', 'time', 'marg'),
                        all.x = T)

postPredFitAll <- postPredFitAll[order(postPredFitAll$city,
                                       postPredFitAll$peak,
                                       postPredFitAll$time,
                                       postPredFitAll$marg, 
                                       postPredFitAll$timePeriod),]

postPredFitAll$marg <- factor(postPredFitAll$marg,
                              levels = c('inc', 'hosp', 'death'),
                              labels = c('Cases', 'Hospitalizations' , 'Deaths'))

postPredFitAll$compartmentType <- ifelse(grepl('SIR', postPredFitAll$modelType),
                                           'SIR', 'SIHRD')
postPredFitAll$alarmType <- 'No alarm'
postPredFitAll$alarmType <- ifelse(grepl('inc', postPredFitAll$modelType),
                                     'Cases only', postPredFitAll$alarmType)

postPredFitAll$alarmType <- ifelse(grepl('full', postPredFitAll$modelType),
                                     'Cases + deaths', postPredFitAll$alarmType)

postPredFitAll$wave <- factor(postPredFitAll$peak,
                              labels = paste0('Wave ', c(1,2)))

postPredFitAll$date <- as.Date(postPredFitAll$date, format = '%m/%d/%Y')


pdf('figures/montreal_postPred.pdf', height = 5, width = 13)
ggplot(postPredFitAll,
       aes(x = date, 
           col = marg, fill = marg, group = marg)) +
    geom_line(aes(y = truth),linewidth = 0.8) +
    geom_line(aes(y = mean), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    facet_nested( wave ~  compartmentType + alarmType , scales = "free", independent = "all") +
    theme_bw() +  
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 16, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = 'Date', y = 'Count',
         col = '', fill = '', 
         title = 'Posterior predictive distribution') + 
    scale_color_manual(values = c('black', 'goldenrod4', 'dodgerblue3'))+ 
    scale_fill_manual(values = c('grey30', 'goldenrod2', 'dodgerblue'))
dev.off()

################################################################################


alarmTimePostAll <- readRDS('./results/alarmTimePostAll.rds')
alarmTimePostAll <- alarmTimePostAll[alarmTimePostAll$city == 'montreal' & 
                                         alarmTimePostAll$peak %in% 1:2 & 
                                         alarmTimePostAll$timePeriod == 6,]


# merge with data for dates
alarmTimePostAll <- merge(alarmTimePostAll, montreal_dat, 
                          by = c('peak', 'time'),
                          all.x = T)


alarmTimePostAll$compartmentType <- ifelse(grepl('SIR', alarmTimePostAll$modelType),
                                    'SIR', 'SIHRD')
alarmTimePostAll$alarmType <- 'No alarm'
alarmTimePostAll$alarmType <- ifelse(grepl('inc', alarmTimePostAll$modelType),
                              'Cases only', alarmTimePostAll$alarmType)

alarmTimePostAll$alarmType <- ifelse(grepl('full', alarmTimePostAll$modelType),
                              'Cases + deaths', alarmTimePostAll$alarmType)

alarmTimePostAll$wave <- factor(alarmTimePostAll$peak,
                              labels = paste0('Wave ', c(1,2)))

alarmTimePostAll$date <- as.Date(alarmTimePostAll$date, format = '%m/%d/%Y')


pdf('figures/montreal_alarmPost.pdf', height = 5, width = 8)
ggplot(alarmTimePostAll,
       aes(x = date, y = mean, ymin = lower, ymax = upper)) +
    geom_line() + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(wave ~ compartmentType + alarmType, scales = "free", independent = "x") +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white')) + 
    guides(col = 'none', fill = 'none') +
    ylim(0, 0.75)+  
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = 'Date', y = 'Alarm',
         col = '', fill = '', 
         title = 'Posterior distribution of alarm over time')
dev.off()

################################################################################
# R0

# need to add actual dates to the R0 and posterior predictions as those are indexed
# by epidemic time

subset(dat_longer, type == 'Cases') %>%
    group_by(city, peak) %>%
    summarize(tau = length(date))

nyc_dat <- data.frame(city = 'NYC',
                      time = c(1:98, 1:227, 1:77), 
                      peak = rep(c(1,2,4), c(98, 227, 77)),
                      date = c(nyc$date[which(nyc$peak == 1)],
                               nyc$date[which(nyc$peak == 2)],
                               nyc$date[which(nyc$peak == 4)]),
                      cases = c(nyc$dailyCases[which(nyc$peak == 1)],
                                nyc$dailyCases[which(nyc$peak == 2)],
                                nyc$dailyCases[which(nyc$peak == 4)]),
                      hosp = c(nyc$dailyHosp[which(nyc$peak == 1)],
                               nyc$dailyHosp[which(nyc$peak == 2)],
                               nyc$dailyHosp[which(nyc$peak == 4)]),
                      death = c(nyc$dailyDeaths[which(nyc$peak == 1)],
                                nyc$dailyDeaths[which(nyc$peak == 2)],
                                nyc$dailyDeaths[which(nyc$peak == 4)]))


montreal_dat <- data.frame(city = 'Montreal',
                           time = c(1:123, 1:329, 1:98), 
                           peak = rep(c(1,2,4), c(123, 329, 98)),
                           date = c(montreal$date[which(montreal$peak == 1)],
                                    montreal$date[which(montreal$peak == 2)],
                                    montreal$date[which(montreal$peak == 4)]),
                           cases = c(montreal$dailyCases[which(montreal$peak == 1)],
                                     montreal$dailyCases[which(montreal$peak == 2)],
                                     montreal$dailyCases[which(montreal$peak == 4)]),
                           hosp = c(montreal$dailyHosp[which(montreal$peak == 1)],
                                    montreal$dailyHosp[which(montreal$peak == 2)],
                                    montreal$dailyHosp[which(montreal$peak == 4)]),
                           death = c(montreal$dailyDeaths[which(montreal$peak == 1)],
                                     montreal$dailyDeaths[which(montreal$peak == 2)],
                                     montreal$dailyDeaths[which(montreal$peak == 4)]))

R0PostAll_nyc <- readRDS('nyc/results/R0PostAll.rds')
R0PostAll_nyc <- merge(R0PostAll_nyc, nyc_dat, by = c('peak', 'time'), all.x = T)

R0PostAll_montreal <- readRDS('montreal/results/R0PostAll.rds')
R0PostAll_montreal <- merge(R0PostAll_montreal, montreal_dat, by = c('peak', 'time'), all.x = T)

R0PostAll <- rbind.data.frame(R0PostAll_nyc, R0PostAll_montreal)

# undetected only
R0PostAll <- R0PostAll[R0PostAll$assumeType == 'undetected',]

R0PostAll$modelType <- factor(R0PostAll$modelType,
                              levels = c('SIHRD_full', 'SIHRD_inc', 
                                         'SIR_full', 'SIR_inc', 
                                         'SIHRD_noAlarm', 'SIR_noAlarm'),
                              labels = c('SIHRD cases + deaths', 
                                         'SIHRD cases only', 
                                         'SIR cases + deaths', 
                                         'SIR cases only', 'SIHRD no alarm',
                                         'SIR no alarm'))

R0PostAll$peak <- paste0('Wave ', R0PostAll$peak)


p1 <- ggplot(subset(R0PostAll, city == 'Montreal' & 
                        modelType %in% c('SIHRD cases + deaths', 
                                         'SIR cases only', 
                                         'SIHRD no alarm')),
             aes(x = date, y = mean, ymin = lower, ymax = upper)) +
    geom_line() + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(peak ~ modelType, scales = "free", independent = "x") +
    theme_bw() + 
    geom_hline(yintercept = 1, linetype = 2) + 
    guides(col = 'none', fill = 'none')  +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = 'Date', y = expression(R[0](t)), title = 'Montreal')+
    scale_x_date(date_breaks ="1 month", date_labels = "%b")

p2 <- ggplot(subset(R0PostAll, city == 'NYC' & 
                        modelType %in% c('SIHRD cases + deaths', 
                                         'SIR cases only', 
                                         'SIHRD no alarm')),
             aes(x = date, y = mean, ymin = lower, ymax = upper)) +
    geom_line() + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(peak ~ modelType, scales = "free", independent = "x") +
    theme_bw() + 
    geom_hline(yintercept = 1, linetype = 2) + 
    guides(col = 'none', fill = 'none')  +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = 'Date', y = NULL,
         title = 'NYC') +
    scale_x_date(date_breaks ="1 month", date_labels = "%b")


pdf('figures/nyc_montreal_r0.pdf', height = 4, width = 9)
grid.arrange(p1, p2, nrow = 1)
dev.off()


################################################################################
# Posterior prediction

# need in long format to merge with posterior predictions

nyc_dat_long <- reshape(nyc_dat, 
                        varying = c("cases", "hosp", "death"), 
                        v.names = "truth",
                        timevar = "marg", 
                        times = c("cases", "hosp", "death"), 
                        new.row.names = 1:10000,
                        direction = "long")

montreal_dat_long <- reshape(montreal_dat, 
                             varying = c("cases", "hosp", "death"), 
                             v.names = "truth",
                             timevar = "marg", 
                             times = c("cases", "hosp", "death"), 
                             new.row.names = 1:10000,
                             direction = "long")



# posterior predictions
postPred_nyc <- readRDS('nyc/results/postPredFitAll.rds')
postPred_nyc <- postPred_nyc[postPred_nyc$marg != 'inc',]
postPred_montreal <- readRDS('montreal/results/postPredFitAll.rds')
postPred_montreal <- postPred_montreal[postPred_montreal$marg != 'inc',]

postPred_nyc <- merge(postPred_nyc, nyc_dat_long, by = c('peak', 'time', 'marg'),
                      all.x = T)
postPred_montreal <- merge(postPred_montreal, montreal_dat_long, by = c('peak', 'time', 'marg'),
                           all.x = T)

postPredAll <- rbind.data.frame(postPred_montreal, postPred_nyc)

postPredAll$modelType <- factor(postPredAll$modelType,
                                levels = c('SIHRD_full', 'SIHRD_inc', 
                                           'SIR_full', 'SIR_inc', 
                                           'SIHRD_noAlarm', 'SIR_noAlarm'),
                                labels = c('SIHRD cases + deaths', 'SIHRD cases only', 
                                           'SIR cases + deaths', 'SIR cases only', 
                                           'SIHRD no alarm', 'SIR no alarm'))


postPredAll$peak <- paste0('Wave ', postPredAll$peak)

p3 <- ggplot(subset(postPredAll, assumeType == 'undetected' & 
                        modelType %in% c('SIHRD cases + deaths', 
                                         'SIR cases only',
                                         'SIHRD no alarm') & 
                        marg == 'cases' & city == 'Montreal'), 
             aes(x = date, col = modelType, fill = modelType)) +
    geom_line(aes(y = truth)) + 
    geom_line(aes(y = mean), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_nested(peak ~ modelType, scales = T, independent = 'x') +
    theme_bw() + 
    guides(col = 'none', fill = 'none') +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_continuous(labels = scales::comma)+ 
    labs(x = 'Date', y = 'Case Count', title = 'Montreal')

p4 <- ggplot(subset(postPredAll, assumeType == 'undetected' & 
                        modelType %in% c('SIHRD cases + deaths', 'SIR cases only', 'SIHRD no alarm') & 
                        marg == 'cases' & city == 'NYC'), 
             aes(x = date, col = modelType, fill = modelType)) +
    geom_line(aes(y = truth)) + 
    geom_line(aes(y = mean), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_nested(peak ~ modelType, scales = T, independent = 'x') +
    theme_bw() + 
    guides(col = 'none', fill = 'none') +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_continuous(labels = scales::comma) + 
    labs(x = 'Date', y = NULL, title = 'NYC')



pdf('figures/nyc_montreal_postPred.pdf', height = 4, width = 9)
grid.arrange(p3, p4, nrow = 1)
dev.off()



p3 <- ggplot(subset(postPredAll, assumeType == 'undetected' & 
                        modelType %in% c('SIHRD cases + deaths', 'SIR cases only', 'SIHRD no alarm') & 
                        marg == 'cases' & city == 'Montreal'), 
             aes(x = date)) +
    geom_line(aes(y = truth)) + 
    geom_line(aes(y = mean), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_nested(peak ~ modelType, scales = T, independent = 'x') +
    theme_bw() + 
    guides(col = 'none', fill = 'none') +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text.x = element_blank(),
          strip.text.y = element_blank(),
          axis.title = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_continuous(labels = scales::comma)+ 
    labs(x = 'Date', y = 'Case Count', title = NULL)+
    scale_x_date(date_breaks ="1 month", date_labels = "%b")

p4 <- ggplot(subset(postPredAll, assumeType == 'undetected' & 
                        modelType %in% c('SIHRD cases + deaths', 'SIR cases only', 'SIHRD no alarm') & 
                        marg %in% c('cases') & 
                        city == 'NYC'), 
             aes(x = date)) +
    geom_line(aes(y = truth)) + 
    geom_line(aes(y = mean), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_nested(peak ~ modelType, scales = T, independent = 'x') +
    theme_bw() + 
    guides(col = 'none', fill = 'none') +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text.y = element_text(size = 12),
          strip.text.x = element_blank(),
          axis.title = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_continuous(labels = scales::comma) + 
    labs(x = 'Date', y = NULL, title = NULL)+
    scale_x_date(date_breaks ="1 month", date_labels = "%b")




pdf('figures/nyc_montreal_resultsComb.pdf', height = 10, width = 13)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()




# hospitalizations and deaths
postPredAll$marg <- factor(postPredAll$marg,
                           levels = c('cases', 'hosp', 'death'),
                           labels = c('cases', 'Hospitalizations', 'Deaths'))

pdf('figures/nyc_montreal_hospDeathPred.pdf', height = 4, width = 8)
ggplot(subset(postPredAll, assumeType == 'undetected' & 
                  modelType %in% c('SIHRD cases + deaths') & 
                  marg %in% c('Deaths', 'Hospitalizations') ), 
       aes(x = date, group = marg, fill = marg)) +
    geom_line(aes(y = truth, col = marg)) + 
    geom_line(aes(y = mean, col = marg), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_nested(city ~ peak , scales = T, independent = 'all') +
    theme_bw() + 
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_continuous(labels = scales::comma) + 
    labs(x = 'Date', y = NULL, title = 'Posterior Prediction of Hospitalizations and Deaths\nSIHRD cases + deaths alarm', col = '', fill = '') +
    scale_color_manual(values = c('goldenrod2', 'royalblue')) +
    scale_fill_manual(values = c('goldenrod2', 'royalblue'))
dev.off()


################################################################################
# alarm over time

alarmTimePostAll <- readRDS('nyc/results/alarmTimePostAll.rds')

# undetected only
alarmTimePostAll <- alarmTimePostAll[alarmTimePostAll$assumeType == 'undetected',]

alarmTimePostAll$modelType <- factor(alarmTimePostAll$modelType,
                                     levels = c('SIHRD_full', 'SIHRD_inc', 
                                                'SIR_full', 'SIR_inc', 
                                                'SIHRD_noAlarm', 'SIR_noAlarm'),
                                     labels = c('SIHRD cases + deaths', 
                                                'SIHRD cases only', 
                                                'SIR cases + deaths', 
                                                'SIR cases only', 'SIHRD no alarm', 'SIR no alarm'))

ggplot(alarmTimePostAll,
       aes(x = time, y = mean, ymin = lower, ymax = upper, 
           col = modelType, fill = modelType)) +
    geom_line() + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(peak~modelType, scales = "free_x", independent = "x") +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white')) + 
    guides(col = 'none', fill = 'none') +
    ylim(c(0, 0.5))

################################################################################
# WAIC

# not comparable between SIR and SIHRD, but comparable between alarms within each


waicAll  <- waicAll [waicAll $city == 'montreal' & 
                         waicAll $peak %in% 1:2 & 
                         waicAll $timePeriod == 6,]
