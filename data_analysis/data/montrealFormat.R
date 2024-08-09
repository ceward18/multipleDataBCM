################################################################################
# montreal data formatting
# data downloaded from :
# https://www.inspq.qc.ca/covid-19/donnees on Aug 21, 2022
# waves defined on their website:
# wave 1: Feb 23 - 11 July 2020
# wave 2: Aug 23, 2020 - March 20, 2021 (delta)
# wave 5: Dec 5, 2021 - March 12, 2022 (omicron)
################################################################################


library(lubridate)

### Read in downloaded data
cases <- read.csv("graph_1-1.csv")[,c(1,3)]
hosps <- read.csv("graph_2-2.csv")
deaths <- read.csv("graph_3-1.csv")

colnames(cases) <- c('date', 'dailyCases')
colnames(hosps) <- c('date', 'dailyHosp')
colnames(deaths) <- c('date', 'dailyDeaths')

# format dates
cases$date <- as.Date(cases$date)
hosps$date <- as.Date(hosps$date)
deaths$date <- as.Date(deaths$date)

# start data at end of Feb instead of january
cases  <- cases[cases$date >= as.Date('2020-02-23'),]
hosps  <- hosps[hosps$date >= as.Date('2020-02-23'),]

# deaths start march 14, fill in data
deathAdd <- data.frame(date = as.Date(as.Date('2020-02-23'):as.Date('2020-03-13')),
                       dailyDeaths = 0)
deaths  <- rbind.data.frame(deathAdd, deaths)

# don't need data after omicron wave
cases  <- cases[cases$date <= as.Date('2022-03-12'),]
hosps  <- hosps[hosps$date <= as.Date('2022-03-12'),]
deaths  <- deaths[deaths$date <= as.Date('2022-03-12'),]

# merge together
montreal <- merge(cases, hosps, by = 'date')
montreal <- merge(montreal, deaths, by = 'date', all.x = T)


# population
montreal$Population <- 1762949

# peak identifier
montreal$peak <- NA

# peak 1 - Feb 23 - 11 July 2020 (but first case observed march 6)
startDate <- as.Date('2020-03-07')
endDate <- as.Date('2020-07-12')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 1

plot(montreal$date[which(montreal$peak == 1)], 
     montreal$dailyCases[which(montreal$peak == 1)], pch = 16, type = 'l')

# peak 2 - Aug 23, 2020 - March 20, 2021
startDate <- as.Date('2020-08-23')
endDate <- as.Date('2021-03-21')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 2

plot(montreal$date[which(montreal$peak == 2)], 
     montreal$dailyCases[which(montreal$peak == 2)], pch = 16, type = 'l')

# peak 3 - Dec 5, 2021 - March 12, 2022
startDate <- as.Date('2021-12-05')
endDate <- as.Date('2022-03-13')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 3

plot(montreal$date[which(montreal$peak == 3)], 
     montreal$dailyCases[which(montreal$peak == 3)], pch = 16, type = 'l')



movingAverage <- function(x, bw) {
    n <- length(x)
    bw <- floor(bw)
    out <- rep(0, n)
    for (i in 1:n) {
        if (i < bw) {
            t1 = 1
            t2 = i
        } else {
            t1 = i - bw + 1
            t2 = i
        }
        out[i] <- mean(x[t1:t2])
    }
    return(out)
}

# smooth due to weekend effects
# keep raw data for first month
montreal$dailyCases <- c(montreal$dailyCases[1:30], 
                         round(movingAverage(montreal$dailyCases, 3))[-c(1:30)])
montreal$dailyHosp <- c(montreal$dailyHosp[1:30], 
                        round(movingAverage(montreal$dailyHosp, 3))[-c(1:30)])
montreal$dailyDeaths <- c(montreal$dailyDeaths[1:30], 
                          round(movingAverage(montreal$dailyDeaths, 3))[-c(1:30)])


write.csv(montreal, 'montrealClean.csv', quote = F, row.names = F)


