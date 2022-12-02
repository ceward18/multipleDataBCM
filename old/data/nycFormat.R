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

# population
nyc$Population <- 8804190


# peak identifier
nyc$peak <- NA

# peak 1 - Feb 29 - Jun 1, 2020
startDate <- as.Date('2020-02-29')
endDate <- as.Date('2020-06-16')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 1

# peak 2 - Oct 01, 2020 - Jun 1, 2021
startDate <- as.Date('2020-10-01')
endDate <- as.Date('2021-05-16')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 2

# peak 3 - Jul 1 - Nov 1, 2021
startDate <- as.Date('2021-07-01')
endDate <- as.Date('2021-11-01')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 3

# peak 4 - Dec 1, 2021 - Feb 15, 2022
startDate <- as.Date('2021-12-01')
endDate <- as.Date('2022-02-16')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 4

# peak 5 - Nov 11, 2021 - Mar 1, 2022
startDate <- as.Date('2022-03-01')
endDate <- as.Date('2022-09-01')
nyc$peak[nyc$date >= startDate & nyc$date < endDate]<- 5

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

nyc$dailyCases <- movingAverage(nyc$dailyCases, 7)
nyc$dailyDeaths <- movingAverage(nyc$dailyDeaths, 7)
nyc$dailyHosp <- movingAverage(nyc$dailyHosp, 7)

write.csv(nyc, 'nycClean.csv', quote = F, row.names = F)


# plot(nyc$date[which(nyc$peak == 5)], nyc$dailyCases[which(nyc$peak == 5)], type = 'l')

