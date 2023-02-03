
################################################################################
# COVID London
# Downloaded 10/23/2021
# https://coronavirus.data.gov.uk/details/download
# newCasesBySpecimenDate
# newDailyNsoDeathsByDeathDate
# newAdmissions
################################################################################


london <- read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&areaCode=E12000007&metric=newCasesBySpecimenDate&metric=newDailyNsoDeathsByDeathDate&format=csv") 
hosp <- read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&areaCode=E40000003&metric=newAdmissions&format=csv") 


london <- merge(london, hosp, by = 'date', all.x = T)

london$Date <- as.Date(london$date, format = '%Y-%m-%d')
london <- london[order(london$Date),]

london <- london[,c('Date', 'newCasesBySpecimenDate', 
                    'newDailyNsoDeathsByDeathDate', 'newAdmissions')]
colnames(london) <- c('date', 'dailyCases', 'dailyDeaths', 'dailyHosp')

# population
london$Population <- 8799800


# replace NA's with 0
london$dailyDeaths[is.na(london$dailyDeaths)] <- 0
london$dailyHosp[is.na(london$dailyHosp)] <- 0


# peak identifier
london$peak <- NA

# peak 1 - Mar 1 - Jun 15, 2020
startDate <- as.Date('2020-03-01')
endDate <- as.Date('2020-06-01')
london$peak[london$date >= startDate & london$date < endDate]<- 1

plot(london$date[which(london$peak == 1)], 
     london$dailyCases[which(london$peak == 1)], pch = 16)

# peak 2 - Sep 15, 2020 - Apr 1, 2021
startDate <- as.Date('2020-09-15')
endDate <- as.Date('2021-03-01')
london$peak[london$date >= startDate & london$date < endDate]<- 2

plot(london$date[which(london$peak == 2)], 
     london$dailyCases[which(london$peak == 2)], pch = 16)

# peak 3 - Jul 1 - Nov 1, 2021
startDate <- as.Date('2021-05-01')
endDate <- as.Date('2021-10-01')
london$peak[london$date >= startDate & london$date < endDate]<- 3


plot(london$date[which(london$peak == 3)], 
     london$dailyCases[which(london$peak == 3)], pch = 16)

# peak 4 - Nov 15, 2021 - Feb 25, 2022
startDate <- as.Date('2021-11-15')
endDate <- as.Date('2022-02-26')
london$peak[london$date >= startDate & london$date < endDate]<- 4

plot(london$date[which(london$peak == 4)], 
     london$dailyCases[which(london$peak == 4)], pch = 16)

# peak 5 - Feb 26, 2022 - May 1, 2022
startDate <- as.Date('2022-02-26')
endDate <- as.Date('2022-05-01')
london$peak[london$date >= startDate & london$date < endDate]<- 5

plot(london$date[which(london$peak == 5)], 
     london$dailyCases[which(london$peak == 5)], pch = 16)


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

london$dailyCases <- movingAverage(london$dailyCases, 7)
london$dailyDeaths <- movingAverage(london$dailyDeaths, 7)
london$dailyHosp <- movingAverage(london$dailyHosp, 7)


write.csv(london, 'londonClean.csv', quote = F, row.names = F)
