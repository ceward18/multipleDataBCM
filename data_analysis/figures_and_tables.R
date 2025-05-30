################################################################################
# Figures and tables in paper
# Data Analysis of Montreal and Miami
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
# Fig 6 - Montreal and Miami cases and deaths data

### set up Montreal data
montreal <- read.csv('montreal/data/montrealClean.csv')
montreal$date <- as.Date(montreal$date)
montreal$Peak <- factor(montreal$peak)

# just 2020
montreal <- montreal[montreal$date >= min(montreal$date[which(montreal$peak == 1)]),]
montreal <- montreal[montreal$date < as.Date('2021-01-01'),]

# time to peak for each wave (weeks) 
# peak 2 has a second peak not of interest
which.max(montreal[which(montreal$peak == 1), 'dailyCases']) / 7
which.max(montreal[which(montreal$peak == 2), 'dailyCases'][1:50]) / 7

# peak 1 - 8 weeks
peak1Dates_montreal <- montreal$date[min(which(montreal$peak == 1))]
peak1Dates_montreal <- c(peak1Dates_montreal, peak1Dates_montreal + 8*7)

# peak 2 - 10 weeks
peak2Dates_montreal <- montreal$date[min(which(montreal$peak == 2))]
peak2Dates_montreal <- c(peak2Dates_montreal, peak2Dates_montreal + 10*7)


### set up miami data
miami <- read.csv('miami/data/miamiClean.csv')
miami$date <- as.Date(miami$date)
miami$Peak <- factor(miami$peak)


# time to peak for each wave (weeks)
which.max(miami[which(miami$peak == 1), 'dailyCases']) / 7
which.max(miami[which(miami$peak == 2), 'dailyCases']) / 7

# peak 1 = 6 weeks
peak1Dates_miami <- miami$date[min(which(miami$peak == 1))]
peak1Dates_miami <- c(peak1Dates_miami, peak1Dates_miami + 6*7)

# peak 2 = 10 weeks
peak2Dates_miami <- miami$date[min(which(miami$peak == 2))]
peak2Dates_miami <- c(peak2Dates_miami, peak2Dates_miami + 10*7)



montreal$city <- 'Montreal'
miami$city <- 'Miami'

# convert to long
montreal_long <- reshape(montreal, 
                         varying = c("dailyCases", "dailyDeaths"), 
                         v.names = "count",
                         timevar = "marg", 
                         times = c("dailyCases", "dailyDeaths"), 
                         new.row.names = 1:10000,
                         direction = "long")


miami_long <- reshape(miami, 
                      varying = c("dailyCases", "dailyDeaths"), 
                      v.names = "count",
                      timevar = "marg", 
                      times = c("dailyCases", "dailyDeaths"), 
                      new.row.names = 1:10000,
                      direction = "long")

dat_long <- rbind.data.frame(montreal_long, miami_long)


dat_long$marg <- factor(dat_long$marg,
                        labels = c('Cases', 'Deaths'))

dat_long$peak1_min <- ifelse(dat_long$city == 'Montreal', peak1Dates_montreal[1], peak1Dates_miami[1])
dat_long$peak1_max <- ifelse(dat_long$city == 'Montreal', peak1Dates_montreal[2], peak1Dates_miami[2])


dat_long$peak2_min <- ifelse(dat_long$city == 'Montreal', peak2Dates_montreal[1], peak2Dates_miami[1])
dat_long$peak2_max <- ifelse(dat_long$city == 'Montreal', peak2Dates_montreal[2], peak2Dates_miami[2])

dat_long$peak1_min <- as.Date(dat_long$peak1_min)
dat_long$peak1_max <- as.Date(dat_long$peak1_max)
dat_long$peak2_min <- as.Date(dat_long$peak2_min)
dat_long$peak2_max <- as.Date(dat_long$peak2_max)

dat_long$marg <- factor(dat_long$marg, 
                        labels = c('Cases', 'Deaths'))


cases_height <- 125
deaths_height <- 5.2

ann_text1 <- data.frame(city = rep(c('Montreal', 'Miami'), 2),
                        marg = rep(c('Cases', 'Deaths'), each = 2),
                        x = rep(c(as.Date('2020-04-06'), as.Date('2020-04-03')), 2),
                        y = c(rep(cases_height, 2), rep(deaths_height, 2)),
                        lab = c(rep('Wave 1', 2), '', ''))

ann_text2 <- data.frame(city = rep(c('Montreal', 'Miami'), 2),
                        marg = rep(c('Cases', 'Deaths'), each = 2),
                        x = rep(c(as.Date('2020-09-25'), as.Date('2020-07-08')), 2),
                        y = c(rep(cases_height, 2), rep(deaths_height, 2)),
                        lab = c(rep('Wave 2', 2), '', ''))


pdf('figures/fig6_montreal_miami_data.pdf', height = 5.5, width = 9)
ggplot(dat_long, aes(x = date, y = count/Population* 1e5)) + 
    geom_rect(aes(xmin = peak1_min, 
                  xmax = peak1_max, ymin = -Inf, ymax = Inf), 
              fill = 'grey90', alpha = 0.2) +
    geom_rect(aes(xmin = peak2_min,
                  xmax = peak2_max, ymin = -Inf, ymax = Inf),
              fill = 'grey90', alpha = 0.2) +
    geom_line(aes(col = marg), linewidth = 0.6) + 
    facet_grid(marg ~ city, scales = 'free_y', switch = 'y') + 
    scale_y_continuous(labels = scales::comma) +
    scale_x_date(breaks = seq(as.Date('2020-03-01'), by ='month', length.out = 11)[seq(1, 11, by = 2)],
                 date_labels = "%b") +
    labs(x = 'Date', y = 'Count per 100,000',
         title = 'COVID-19 Case and Death Rates in 2020') +
    geom_text(data = ann_text1, aes(x = x, y = y, label = lab), size = 4) +
    geom_text(data = ann_text2, aes(x = x, y = y, label = lab), size = 4) +
    theme_bw() + 
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 12, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 11),
          strip.placement = "outside",
          strip.switch.pad.grid = unit(0.2, "cm"),
          axis.title.y = element_text(vjust = -12, size = 11)) + 
    scale_color_manual(values = c('black', 'red')) + 
    guides(col = 'none')
dev.off()

################################################################################
# Table 1 - posterior alpha parameters by wave/city


paramsPostAll_miami <- readRDS('miami/results/paramsPostAll.rds')

# 6 weeks for wave 1 and 10 weeks for wave 2
alphaPost_miami <- paramsPostAll_miami[paramsPostAll_miami$param == 'alpha' & 
                                           paramsPostAll_miami$peak %in% 1:2 & 
                                           paramsPostAll_miami$timePeriod %in% c(6, 10),]

alphaPost_miami <- alphaPost_miami[which(alphaPost_miami$peak == 1 & alphaPost_miami$timePeriod == 6 | 
                                             alphaPost_miami$peak == 2 & alphaPost_miami$timePeriod == 10  ),]

paramsPostAll_montreal <- readRDS('montreal/results/paramsPostAll.rds')

# 8 weeks for wave 1 and 10 weeks for wave 2
alphaPost_montreal <- paramsPostAll_montreal[paramsPostAll_montreal$param == 'alpha' & 
                                                 paramsPostAll_montreal$peak %in% 1:2 & 
                                                 paramsPostAll_montreal$timePeriod %in% c(8, 10),]

alphaPost_montreal <- alphaPost_montreal[which(alphaPost_montreal$peak == 1 & alphaPost_montreal$timePeriod == 8 | 
                                                   alphaPost_montreal$peak == 2 & alphaPost_montreal$timePeriod == 10  ),]

alphaPost <- rbind.data.frame(alphaPost_miami, alphaPost_montreal)

alphaPost$compartmentModel <- ifelse(alphaPost$modelType == 'SIHRD_full', 
                                     'SIHRD', 'SIR')

alphaPost$city <- factor(alphaPost$city,
                         labels = c("Miami", "Montreal"))


# table
alphaPost$val <- paste0(sprintf("%.3f", round(alphaPost$mean, 3)), 
                        ' (',
                        sprintf("%.3f", round(alphaPost$lower, 3)),
                        ', ',
                        sprintf("%.3f", round(alphaPost$upper, 3)),
                        ')')


alphaPostTab <- alphaPost[,c('city', 'compartmentModel', 'peak', 'val')]

alphaPostTab <- reshape(alphaPostTab, 
                        timevar = "peak",
                        idvar = c("city", "compartmentModel"),
                        direction = "wide")


kable(alphaPostTab, row.names = F, format = 'latex', align = 'llcc', 
      booktabs = T, escape = F, 
      col.names = linebreak(c('\\textbf{City}',
                              '\\textbf{Model fitted}', 
                              '\\textbf{Wave 1}',
                              '\\textbf{Wave2}'), align = 'c')) %>% 
    collapse_rows(columns = 1, latex_hline = 'major') 





################################################################################
# Fig 7 - alarm functions and R0 over time


# get epidemic time matched with dates
miami_peak_lengths <- c(sum(miami$peak == 1, na.rm = T),
                        sum(miami$peak == 2, na.rm = T))

miami_dat <- data.frame(city = 'miami',
                        time = c(1:miami_peak_lengths[1], 
                                 1:miami_peak_lengths[2]), 
                        peak = rep(c(1,2), miami_peak_lengths),
                        date = c(miami$date[which(miami$peak == 1)],
                                 miami$date[which(miami$peak == 2)]),
                        inc = c(miami$dailyCases[which(miami$peak == 1)],
                                miami$dailyCases[which(miami$peak == 2)]),
                        hosp = c(miami$dailyHosp[which(miami$peak == 1)],
                                 miami$dailyHosp[which(miami$peak == 2)]),
                        death = c(miami$dailyDeaths[which(miami$peak == 1)],
                                  miami$dailyDeaths[which(miami$peak == 2)]))


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


### alarms
alarmTimePostAll_miami <- readRDS('miami/results/alarmTimePostAll.rds')

# 6 weeks for wave 1 and 10 weeks for wave 2
alarmTimePostAll_miami <- alarmTimePostAll_miami[alarmTimePostAll_miami$peak %in% 1:2 & 
                                                     alarmTimePostAll_miami$timePeriod %in% c(6, 10),]

alarmTimePostAll_miami <- alarmTimePostAll_miami[which(alarmTimePostAll_miami$peak == 1 & 
                                                           alarmTimePostAll_miami$timePeriod == 6 | 
                                                           alarmTimePostAll_miami$peak == 2 & 
                                                           alarmTimePostAll_miami$timePeriod == 10),]


alarmTimePostAll_montreal <- readRDS('montreal/results/alarmTimePostAll.rds')

# 8 weeks for wave 1 and 10 weeks for wave 2
alarmTimePostAll_montreal <- alarmTimePostAll_montreal[alarmTimePostAll_montreal$peak %in% 1:2 & 
                                                           alarmTimePostAll_montreal$timePeriod %in% c(8, 10),]

alarmTimePostAll_montreal <- alarmTimePostAll_montreal[which(alarmTimePostAll_montreal$peak == 1 & 
                                                                 alarmTimePostAll_montreal$timePeriod == 8 | 
                                                                 alarmTimePostAll_montreal$peak == 2 & 
                                                                 alarmTimePostAll_montreal$timePeriod == 10),]


# merge with data for dates
alarmTimePostAll_miami <- merge(alarmTimePostAll_miami, miami_dat, 
                                by = c('city', 'peak', 'time'),
                                all.x = T)
alarmTimePostAll_montreal <- merge(alarmTimePostAll_montreal, montreal_dat, 
                                   by = c('city','peak', 'time'),
                                   all.x = T)

# merge cities
alarmTimePostAll <- rbind.data.frame(alarmTimePostAll_miami, alarmTimePostAll_montreal)
alarmTimePostAll$city <- factor(alarmTimePostAll$city,
                                labels = c("Miami", "Montreal"))


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



### R0


R0PostAll_miami <- readRDS('miami/results/R0PostAll.rds')

# 6 weeks for wave 1 and 10 weeks for wave 2
R0PostAll_miami <- R0PostAll_miami[R0PostAll_miami$peak %in% 1:2 & 
                                       R0PostAll_miami$timePeriod %in% c(6, 10),]

R0PostAll_miami <- R0PostAll_miami[which(R0PostAll_miami$peak == 1 & 
                                             R0PostAll_miami$timePeriod == 6 | 
                                             R0PostAll_miami$peak == 2 & 
                                             R0PostAll_miami$timePeriod == 10),]


R0PostAll_montreal <- readRDS('montreal/results/R0PostAll.rds')

# 8 weeks for wave 1 and 10 weeks for wave 2
R0PostAll_montreal <- R0PostAll_montreal[R0PostAll_montreal$peak %in% 1:2 & 
                                             R0PostAll_montreal$timePeriod %in% c(8, 10),]

R0PostAll_montreal <- R0PostAll_montreal[which(R0PostAll_montreal$peak == 1 & 
                                                   R0PostAll_montreal$timePeriod == 8 | 
                                                   R0PostAll_montreal$peak == 2 & 
                                                   R0PostAll_montreal$timePeriod == 10),]

# merge with data for dates
R0PostAll_miami <- merge(R0PostAll_miami, miami_dat, 
                         by = c('city', 'peak', 'time'),
                         all.x = T)
R0PostAll_montreal <- merge(R0PostAll_montreal, montreal_dat, 
                            by = c('city','peak', 'time'),
                            all.x = T)

# merge cities
R0PostAll <- rbind.data.frame(R0PostAll_miami, R0PostAll_montreal)
R0PostAll$city <- factor(R0PostAll$city,
                         labels = c("Miami", "Montreal"))


R0PostAll$compartmentType <- ifelse(grepl('SIR', R0PostAll$modelType),
                                    'SIR', 'SIHRD')
R0PostAll$alarmType <- 'No alarm'
R0PostAll$alarmType <- ifelse(grepl('inc', R0PostAll$modelType),
                              'Cases only', R0PostAll$alarmType)

R0PostAll$alarmType <- ifelse(grepl('full', R0PostAll$modelType),
                              'Cases + deaths', R0PostAll$alarmType)

R0PostAll$wave <- factor(R0PostAll$peak,
                         labels = paste0('Wave ', c(1,2)))

R0PostAll$date <- as.Date(R0PostAll$date, format = '%m/%d/%Y')



### create plots


myTheme <- theme(strip.placement = "outside",
                 strip.background = element_blank(),
                 strip.text = element_text(size = 13),
                 axis.title = element_text(size =11),
                 axis.text = element_text(size = 10),
                 plot.title = element_text(size = 14, h = 0.5),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 strip.clip = "off") 


p_alarm <- ggplot(alarmTimePostAll,
                  aes(x = time, y = mean, ymin = lower, ymax = upper,
                      group = wave, fill = wave)) +
    geom_line(aes(col = wave)) + 
    geom_ribbon(alpha = 0.3) +
    facet_nested( city  ~  compartmentType + alarmType, 
                  scales = "free", independent = 'x') +
    theme_bw() + 
    myTheme +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2),
                       labels = seq(0, 1, 0.2)) +
    labs(x = 'Time since wave start', y = 'Alarm',
         col = '', fill = '', 
         title = 'Posterior distribution of alarm over time')

p_r0 <- ggplot(R0PostAll,
               aes(x = time, y = mean, ymin = lower, ymax = upper,
                   group = wave, fill = wave)) +
    geom_hline(yintercept = 1, linetype = 2) + 
    geom_line(aes(col = wave)) + 
    geom_ribbon(alpha = 0.3) +
    facet_nested( city ~  compartmentType + alarmType, 
                  scales = "free", independent = 'x') +
    theme_bw() + 
    theme(strip.background = element_rect(fill = 'white')) + 
    myTheme + 
    scale_y_continuous(limits = c(0.5, 2.5)) +
    labs(x = 'Time since wave start', y = expression(R[0](t)),
         col = '', fill = '', 
         title = expression('Posterior distribution of'~R[0](t)))

pdf('figures/fig7_data_alarmR0Post2.pdf', height = 8, width = 10)
grid.arrange(p_alarm, p_r0, nrow = 2)
dev.off()





################################################################################
# Figure 8 - Posterior predictive fit + alarms


postPredFitAll_miami <- readRDS('miami/results/postPredFitAll.rds')

# 6 weeks for wave 1 and 8 weeks for wave 2
postPredFitAll_miami <- postPredFitAll_miami[postPredFitAll_miami$peak %in% 1:2 & 
                                                 postPredFitAll_miami$timePeriod %in% c(6, 10),]

postPredFitAll_miami <- postPredFitAll_miami[which(postPredFitAll_miami$peak == 1 &
                                                       postPredFitAll_miami$timePeriod == 6 | 
                                                       postPredFitAll_miami$peak == 2 &
                                                       postPredFitAll_miami$timePeriod == 10),]

postPredFitAll_montreal <- readRDS('montreal/results/postPredFitAll.rds')

# 8 weeks for wave 1 and 8 weeks for wave 2
postPredFitAll_montreal <- postPredFitAll_montreal[postPredFitAll_montreal$peak %in% 1:2 & 
                                                       postPredFitAll_montreal$timePeriod %in% c(8, 10),]


postPredFitAll_montreal <- postPredFitAll_montreal[which(postPredFitAll_montreal$peak == 1 &
                                                             postPredFitAll_montreal$timePeriod == 8 | 
                                                             postPredFitAll_montreal$peak == 2 &
                                                             postPredFitAll_montreal$timePeriod == 10),]



# need in long format to merge with posterior predictions
miami_dat_long <- reshape(miami_dat, 
                          varying = c("inc", "hosp", "death"), 
                          v.names = "truth",
                          timevar = "marg", 
                          times = c("inc", "hosp", "death"), 
                          new.row.names = 1:10000,
                          direction = "long")

miami_dat_long <- miami_dat_long[,-which(colnames(miami_dat_long) == 'id')]

postPredFitAll_miami <- merge(postPredFitAll_miami, miami_dat_long, 
                              by = c('city', 'peak', 'time', 'marg'),
                              all.x = T)



# need in long format to merge with posterior predictions

montreal_dat_long <- reshape(montreal_dat, 
                             varying = c("inc", "hosp", "death"), 
                             v.names = "truth",
                             timevar = "marg", 
                             times = c("inc", "hosp", "death"), 
                             new.row.names = 1:10000,
                             direction = "long")

montreal_dat_long <- montreal_dat_long[,-which(colnames(montreal_dat_long) == 'id')]

postPredFitAll_montreal <- merge(postPredFitAll_montreal, montreal_dat_long, 
                                 by = c('city', 'peak', 'time', 'marg'),
                                 all.x = T)

postPredFitAll <- rbind.data.frame(postPredFitAll_miami, postPredFitAll_montreal)

postPredFitAll <- postPredFitAll[order(postPredFitAll$city,
                                       postPredFitAll$peak,
                                       postPredFitAll$time,
                                       postPredFitAll$marg, 
                                       postPredFitAll$timePeriod),]


# convert outcome to rate scale
postPredFitAll$population <- ifelse(postPredFitAll$city == 'miami', 
                                    miami$Population[1],
                                    montreal$Population[1])
postPredFitAll$truth <- postPredFitAll$truth/postPredFitAll$population * 1e5
postPredFitAll$mean <- postPredFitAll$mean/postPredFitAll$population * 1e5
postPredFitAll$lower <- postPredFitAll$lower/postPredFitAll$population * 1e5
postPredFitAll$upper <- postPredFitAll$upper/postPredFitAll$population * 1e5


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

postPredFitAll$city <- factor(postPredFitAll$city,
                              labels = c("Miami", "Montreal"))

scales_set<- list(
    scale_x_date(date_labels = "%b %d ",date_breaks = '2 weeks'),
    scale_x_date(date_labels = "%b %d ",date_breaks = '2 weeks'),
    scale_x_date(date_labels = "%b %d ",date_breaks = '2 weeks'),
    scale_x_date(date_labels = "%b %d ",date_breaks = '2 weeks'),
    scale_x_date(date_labels = "%b %d ",date_breaks = '2 weeks'),
    
    scale_x_date(date_labels = "%b %d ",  date_breaks = '4 weeks'),
    scale_x_date(date_labels = "%b %d ",  date_breaks = '4 weeks'),
    scale_x_date(date_labels = "%b %d ",  date_breaks = '4 weeks'),
    scale_x_date(date_labels = "%b %d ",  date_breaks = '4 weeks'),
    scale_x_date(date_labels = "%b %d ",  date_breaks = '4 weeks'),
    
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks'),
    scale_x_date(date_labels = "%b %d ", date_breaks = '3 weeks')
)



pdf('figures/fig8_data_postPred.pdf', height = 7, width = 12)
ggplot(subset(postPredFitAll, marg %in% c('Cases', 'Deaths')),
       aes(x = date, fill = marg, group = marg)) +
    geom_line(aes(y = truth),linewidth = 0.8) +
    geom_line(aes(y = mean, col = marg), linetype = 2) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    facet_nested( city + wave ~  compartmentType + alarmType, 
                  scales = "free", independent = 'x') +
    facetted_pos_scales(x = scales_set)+
    theme_bw() +  
    scale_x_date(date_labels = "%b %d ",
                 date_breaks = '2 weeks') +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 13),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 9),
          plot.title = element_text(size = 14, h = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = 'Date', y = 'Count per 100,000',
         col = '', fill = '', 
         title = 'Posterior predictive distributions') + 
    scale_color_manual(values = c('black',  'dodgerblue3'))+ 
    scale_fill_manual(values = c('grey30', 'dodgerblue'))
dev.off()

################################################################################
# Table 2 - WAIC values

# not comparable between SIR and SIHRD, but comparable between alarms within each

waicAll_miami <- readRDS('miami/results/waicAll.rds')

# 6 weeks for wave 1 and 10 weeks for wave 2
waicAll_miami <- waicAll_miami[waicAll_miami$peak %in% 1:2 & 
                                   waicAll_miami$timePeriod %in% c(6, 10),]

waicAll_miami <- waicAll_miami[which(waicAll_miami$peak == 1 & 
                                         waicAll_miami$timePeriod == 6 | 
                                         waicAll_miami$peak == 2 & 
                                         waicAll_miami$timePeriod == 10  ),]


waicAll_montreal <- readRDS('montreal/results/waicAll.rds')

# 8 weeks for wave 1 and 10 weeks for wave 2
waicAll_montreal <- waicAll_montreal[waicAll_montreal$peak %in% 1:2 & 
                                         waicAll_montreal$timePeriod %in% c(8, 10),]

waicAll_montreal <- waicAll_montreal[which(waicAll_montreal$peak == 1 & 
                                               waicAll_montreal$timePeriod == 8 | 
                                               waicAll_montreal$peak == 2 & 
                                               waicAll_montreal$timePeriod == 10  ),]


waicAll <- rbind.data.frame(waicAll_miami, waicAll_montreal)


waicAll$compartmentType <- ifelse(grepl('SIR', waicAll$modelType),
                                  'SIR', 'SIHRD')
waicAll$alarmType <- 'No alarm'
waicAll$alarmType <- ifelse(grepl('inc', waicAll$modelType),
                            'Cases only', waicAll$alarmType)

waicAll$alarmType <- ifelse(grepl('full', waicAll$modelType),
                            'Cases + deaths', waicAll$alarmType)

waicAll$wave <- factor(waicAll$peak,
                       labels = paste0('Wave ', c(1,2)))


waicAll <- waicAll[,c('city','wave',  'compartmentType', 'alarmType', 'waic')]


waicTab <- reshape(waicAll, 
                   timevar = "alarmType",
                   idvar = c("city", "wave", 'compartmentType'),
                   direction = "wide")

waicTab <- reshape(waicTab, 
                   timevar = "compartmentType",
                   idvar = c("city", "wave"),
                   direction = "wide")



waicTab <- waicTab[order(waicTab$city, 
                         waicTab$wave),]
waicTab$city <- factor(waicTab$city,
                       labels = c('Miami', 'Montreal'))


kable(waicTab, row.names = F, format = 'latex', align = 'llcccccc', 
      booktabs = T, escape = F,   digits = 1,
      col.names = linebreak(c('\\textbf{City}',
                              '\\textbf{Wave}', 
                              '\\textbf{Cases + deaths}',
                              '\\textbf{Cases only}',
                              '\\textbf{No alarm}', 
                              '\\textbf{Cases + deaths}',
                              '\\textbf{Cases only}',
                              '\\textbf{No alarm}'), align = 'c')) %>% 
    collapse_rows(columns = 1, latex_hline = 'major') %>%
    add_header_above(c(" " = 2, "SIHRD" = 3, "SIR" = 3))
