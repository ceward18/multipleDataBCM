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
# Plot data

### set up NYC data
nyc <- read.csv('nyc/data/nycClean.csv')
nyc$date <- as.Date(nyc$date, format = '%m/%d/%Y')
nyc$Peak <- factor(nyc$peak)
nyc$dailyCases <- round(nyc$dailyCases)
nyc$dailyHosp <- round(nyc$dailyHosp)
nyc$dailyDeaths <- round(nyc$dailyDeaths)


### set up Montreal data
montreal <- read.csv('montreal/data/montrealClean.csv')
montreal$date <- as.Date(montreal$date)
montreal$Peak <- factor(montreal$peak)

# start and end both on the same days

# start
nyc <- nyc[nyc$date >= min(nyc$date[which(nyc$peak == 1)]),]
montreal <- montreal[montreal$date >= min(nyc$date[which(nyc$peak == 1)]),]

# end
nyc <- nyc[nyc$date <= (max(montreal$date[which(montreal$peak == 4)]) + 10),]
montreal <- montreal[montreal$date <= (max(montreal$date[which(montreal$peak == 4)]) + 10),]


pal <- c('red', 'blue', 'goldenrod1')

p1 <- ggplot(nyc, aes(x = date, y = dailyCases)) + 
    geom_line(linetype = 2) + 
    geom_line(data = subset(nyc, peak == 1), col = pal[1], lwd = 1.2) + 
    geom_line(data = subset(nyc, peak == 2), col = pal[2], lwd = 1.2) + 
    geom_line(data = subset(nyc, peak == 4), col = pal[3], lwd = 1.2) + 
    scale_y_continuous(labels = scales::comma, limits = c(0, 44000)) +
    scale_x_date(date_breaks = "3 month", date_minor_breaks = "1 month",
                 date_labels = "%b '%y") +
    labs(x = 'Date', y = 'Incidence', title = 'NYC COVID-19 Case Counts') +
    annotate('text', x = as.Date('2020-04-15'), y = 8000, label = 'Wave 1') +
    annotate('text', x = as.Date('2021-01-15'), y = 8000, label = 'Wave 2') +
    annotate('text', x = as.Date('2021-12-30'), y = 42500, label = 'Wave 3') +
    theme_bw() + 
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          plot.title = element_text(size = 11, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())


p2 <- ggplot(montreal, aes(x = date, y = dailyCases)) + 
    geom_line(linetype = 2) + 
    geom_line(data = subset(montreal, peak == 1), col = pal[1], lwd = 1.2) + 
    geom_line(data = subset(montreal, peak == 2), col = pal[2], lwd = 1.2) + 
    geom_line(data = subset(montreal, peak == 4), col = pal[3], lwd = 1.2) + 
    scale_y_continuous(labels = scales::comma, limits = c(0, 5100)) +
    scale_x_date(date_breaks = "3 month", date_minor_breaks = "1 month",
                 date_labels = "%b '%y") +
    labs(x = 'Date', y = 'Incidence', title = 'Montreal COVID-19 Case Counts') +
    annotate('text', x = as.Date('2020-04-15'), y = 1000, label = 'Wave 1') +
    annotate('text', x = as.Date('2021-01-15'), y = 1500, label = 'Wave 2') +
    annotate('text', x = as.Date('2021-12-30'), y = 5000, label = 'Wave 3') +
    theme_bw() + 
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          plot.title = element_text(size = 11, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())

pdf('figures/nyc_montreal_data.pdf', height = 7, width = 7)
grid.arrange(p1, p2, nrow =2 )
dev.off()

# combine into one long dataset
nyc$city <- 'NYC'
montreal$city <- 'Montreal'
dat <- rbind.data.frame(nyc, montreal)

dat_longer <- reshape(dat, 
             varying = c("dailyCases", "dailyHosp", "dailyDeaths"), 
             v.names = "counts",
             timevar = "type", 
             times = c("dailyCases", "dailyHosp", "dailyDeaths"), 
             new.row.names = 1:5000,
             direction = "long")
dat_longer <- dat_longer[order(dat_longer$date, dat_longer$city),]

dat_longer$type <- factor(dat_longer$type,
                          levels = c('dailyCases', 'dailyHosp', 'dailyDeaths'),
                          labels = c('Cases', 'Hospitalizations', 'Deaths'))

dat_longer$Peak <- factor(dat_longer$peak, 
                          levels = c(1,2,4),
                          labels = c('Wave 1', 'Wave 2', 'Wave 3'))

pdf('figures/nyc_montreal_data2.pdf', height = 4, width = 8)
ggplot(subset(dat_longer, peak %in% c(1, 2, 4)), 
       aes(x = date,  y = counts, group = type, col = type)) + 
    geom_line(linewidth = 0.7) +
    facet_nested_wrap(city~Peak, scales = 'free') +
    scale_color_manual(values = c('black', 'goldenrod2', 'royalblue')) +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(col = '', y = 'Count', x = 'Date') +
    scale_y_continuous(labels = scales::comma) 
dev.off()

################################################################################
# table of posterior alpha parameters by wave

paramsPostAll_nyc <- readRDS('nyc/results/paramsPostAll.rds')
paramsPostAll_montreal <- readRDS('montreal/results/paramsPostAll.rds')

alphaPost_nyc <- paramsPostAll_nyc[paramsPostAll_nyc$param == 'alpha' & 
                                       paramsPostAll_nyc$assumeType == 'undetected',]
alphaPost_montreal <- paramsPostAll_montreal[paramsPostAll_montreal$param == 'alpha' & 
                                                 paramsPostAll_montreal$assumeType == 'undetected',]

tabPrint <- cbind.data.frame(Model = factor(alphaPost_montreal$modelType,
                                            labels = c('SIHRD', 'SIR')),
                             Wave = factor(alphaPost_montreal$peak, 
                                           labels = c('Wave 1', 'Wave 2', 'Wave 3')),
                             Montreal = paste0(sprintf("%.3f", alphaPost_montreal$mean),
                                               ' (',
                                               sprintf("%.3f", alphaPost_montreal$lower),
                                               ', ',
                                               sprintf("%.3f", alphaPost_montreal$upper), ')'),
                             NYC = paste0(sprintf("%.3f", alphaPost_nyc$mean),
                                               ' (',
                                          sprintf("%.3f", alphaPost_nyc$lower),
                                               ', ',
                                          sprintf("%.3f", alphaPost_nyc$upper), ')'))


kable(tabPrint, row.names = F, format = 'latex', align = 'lccc', 
      booktabs = T, escape = F) %>% 
    collapse_rows(columns = 1, latex_hline = 'major') 


################################################################################
# R0

R0PostAll <- readRDS('nyc/results/R0PostAll.rds')

# undetected only
R0PostAll <- R0PostAll[R0PostAll$assumeType == 'undetected',]

R0PostAll$modelType <- factor(R0PostAll$modelType,
                              levels = c('SIHRD_full', 'SIHRD_inc', 
                                         'SIR_full', 'SIR_inc', 
                                         'SIHRD_noAlarm', 'SIR_noAlarm'),
                              labels = c('SIHRD inc+deaths', 
                                         'SIHRD inc', 
                                         'SIR inc+deaths', 
                                         'SIR inc', 'SIHRD no alarm', 'SIR no alarm'))

ggplot(R0PostAll,
       aes(x = time, y = mean, ymin = lower, ymax = upper, 
           col = modelType, fill = modelType)) +
    geom_line() + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(peak~modelType, scales = "free_x", independent = "x") +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white')) + 
    geom_hline(yintercept = 1, linetype = 2) + 
    guides(col = 'none', fill = 'none')


################################################################################
# Posterior prediction

# need to add actual dates to the posterior predictions as those are indexed
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

nyc_dat_long <- reshape(nyc_dat, 
                   varying = c("cases", "hosp", "death"), 
                   v.names = "truth",
                   timevar = "marg", 
                   times = c("cases", "hosp", "death"), 
                   new.row.names = 1:10000,
                   direction = "long")

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


ggplot(subset(postPred_montreal, assumeType == 'undetected'), 
       aes(x = date, col = marg, group = marg, fill = marg)) +
    geom_line(aes(y = truth)) + 
    geom_line(aes(y = mean), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_nested(peak ~ modelType, scales = T, independent = 'x') +
    theme_bw()

ggplot(subset(postPred_nyc, assumeType == 'undetected'), 
       aes(x = date, col = marg, group = marg, fill = marg)) +
    geom_line(aes(y = truth)) + 
    geom_line(aes(y = mean), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    facet_nested(peak ~ modelType, scales = T, independent = 'x') +
    theme_bw()


################################################################################
# alarm over time

alarmTimePostAll <- readRDS('montreal/results/alarmTimePostAll.rds')

# undetected only
alarmTimePostAll <- alarmTimePostAll[alarmTimePostAll$assumeType == 'undetected',]

alarmTimePostAll$modelType <- factor(alarmTimePostAll$modelType,
                              levels = c('SIHRD_full', 'SIHRD_inc', 
                                         'SIR_full', 'SIR_inc', 
                                         'SIHRD_noAlarm', 'SIR_noAlarm'),
                              labels = c('SIHRD inc+deaths', 
                                         'SIHRD inc', 
                                         'SIR inc+deaths', 
                                         'SIR inc', 'SIHRD no alarm', 'SIR no alarm'))

ggplot(alarmTimePostAll,
       aes(x = time, y = mean, ymin = lower, ymax = upper, 
           col = modelType, fill = modelType)) +
    geom_line() + 
    geom_ribbon(alpha = 0.3) +
    facet_nested(peak~modelType, scales = "free_x", independent = "x") +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white')) + 
    guides(col = 'none', fill = 'none') +
    ylim(c(0, 1))
