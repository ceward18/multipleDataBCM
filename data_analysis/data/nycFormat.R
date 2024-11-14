################################################################################
# NYC data formatting
################################################################################


library(lubridate)

### Read in data from GitHub
nycCases <-read.csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/master/trends/cases-by-day.csv") 
nycHosp <-read.csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/master/trends/hosp-by-day.csv") 
nycDeath <-read.csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/master/trends/deaths-by-day.csv") 

nycCases <- nycCases[,c('date_of_interest', 'CASE_COUNT')]
colnames(nycCases) <- c('date', 'dailyCases')

nycHosp <- nycHosp[,c('date_of_interest', 'HOSPITALIZED_COUNT')]
colnames(nycHosp) <- c('date', 'dailyHosp')

nycDeath <- nycDeath[,c('date_of_interest', 'DEATH_COUNT')]
colnames(nycDeath) <- c('date', 'dailyDeaths')

nyc <- merge(nycCases, nycHosp, by = 'date')
nyc <- merge(nyc, nycDeath, by = 'date')

# format dates
nyc$date <- as.Date(nyc$date, format = '%m/%d/%Y')
nyc <- nyc[order(nyc$date),]

# restrict dates as not interested after omicron
nyc  <- nyc[nyc$date <= as.Date('2022-03-30'),]

# population
nyc$Population <- 8804190


# peak identifier
nyc$peak <- NA

# peak 1 - Mar 2 - Jun 30, 2020
# montreal: Feb 23 - 11 July 2020
startDate <- as.Date('2020-03-02')
endDate <- as.Date('2020-07-01')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 1

plot(nyc$date[which(nyc$peak == 1)], 
     nyc$dailyCases[which(nyc$peak == 1)], pch = 16, type = 'l')


# peak 2 - Sept 01, 2020 - Jun 1, 2021
# montreal: Aug 23, 2020 - March 20, 2021 (delta)
startDate <- as.Date('2020-09-01')
endDate <- as.Date('2021-03-21')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 2

plot(nyc$date[which(nyc$peak == 2)], 
     nyc$dailyCases[which(nyc$peak == 2)], pch = 16, type = 'l')


# peak 3 - Dec 1, 2021 - March 12, 2022
# montreal: Dec 5, 2021 - March 12, 2022
startDate <- as.Date('2021-12-01')
endDate <- as.Date('2022-03-13')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 3


plot(nyc$date[which(nyc$peak == 3)], 
     nyc$dailyCases[which(nyc$peak == 3)], pch = 16, type = 'l')


# smooth data
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

# smooth out weekend reporting effects
# don't smooth deaths as those are not reported dates,they are from death files
nyc$dailyCases <- ceiling(movingAverage(nyc$dailyCases, 7))
nyc$dailyHosp <- ceiling(movingAverage(nyc$dailyHosp, 7))


par(mfrow = c(1,3))
plot(nyc$date[which(nyc$peak == 1)], 
     nyc$dailyCases[which(nyc$peak == 1)], pch = 16, type = 'l')

plot(nyc$date[which(nyc$peak == 2)], 
     nyc$dailyCases[which(nyc$peak == 2)], pch = 16, type = 'l')

plot(nyc$date[which(nyc$peak == 3)], 
     nyc$dailyCases[which(nyc$peak == 3)], pch = 16, type = 'l')



write.csv(nyc, 'nycClean.csv', quote = F, row.names = F)


# plot(nyc$date[which(nyc$peak == 5)], nyc$dailyCases[which(nyc$peak == 5)], type = 'l')

