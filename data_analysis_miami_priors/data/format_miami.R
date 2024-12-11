################################################################################
# Miami-Data county data formatting
# hospitalizations from https://rwilli5.github.io/MiamiCovidProject/Trajectory/
################################################################################

library(lubridate)

# cases and deaths
miami <-read.csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties-2020.csv") 
miami <- miami[grepl('Miami-Dade', miami$county),]

# cases and deaths are cumulative
miami$dailyCases <- diff(c(0, miami$cases))
miami$dailyDeaths <- diff(c(0, miami$deaths))
miami$dailyDeaths[which(miami$dailyDeaths< 0)] <- 0

# hospitalizations
hosp <- read.csv('miami_hospitalizations.csv')
hosp$date <- as.Date(hosp$Date, format = '%m/%d/%Y')
hosp$dailyHosp <- c(hosp$AdmitPrevDay[-1], NA)
hosp <- hosp[,c('date', 'dailyHosp')]

# merge hospitalizations
miami$date <- as.Date(miami$date)
miami <- miami[,c('date', 'dailyCases', 'dailyDeaths')]
miami <- merge(miami, hosp, all.x = T)

# restrict dates as not interested after 2020
miami  <- miami[miami$date <= as.Date('2021-01-01'),]

# population (2020 census)
miami$Population <- 2701767

# smooth out weekend reporting effects

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


miami$dailyCases <- ceiling(movingAverage(miami$dailyCases, 7))
miami$dailyDeaths <- ceiling(movingAverage(miami$dailyDeaths, 7))

plot(miami$date, 
     miami$dailyCases, pch = 16, type = 'l')

# peak identifier
miami$peak <- NA

# peak 1 - Mar 12 - May 30, 2020
startDate <- as.Date('2020-03-13')
endDate <- as.Date('2020-06-01')
miami$peak[miami$date >= startDate & miami$date < endDate]<- 1

plot(miami$date[which(miami$peak == 1)], 
     miami$dailyCases[which(miami$peak == 1)], pch = 16, type = 'l')


# peak 2 - Jun 01, 2020 - Oct 30, 2020
startDate <- as.Date('2020-06-01')
endDate <- as.Date('2020-10-30')
miami$peak[miami$date >= startDate & miami$date < endDate]<- 2

plot(miami$date[which(miami$peak == 2)], 
     miami$dailyCases[which(miami$peak == 2)], pch = 16, type = 'l')

write.csv(miami, 'miamiClean.csv', quote = F, row.names = F)

