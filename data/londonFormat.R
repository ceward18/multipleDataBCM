
################################################################################
# COVID London
# Downloaded 10/23/2021
# https://coronavirus.data.gov.uk/details/download
# newCasesBySpecimenDate


################################################################################


london <- read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=region&areaCode=E12000007&metric=newCasesBySpecimenDate&metric=newDailyNsoDeathsByDeathDate&format=csv") 
hosp <- read.csv("https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&areaCode=E40000003&metric=newAdmissions&format=csv") 


london <- merge(london, hosp, by = 'date', all.x = T)

london$Date <- as.Date(london$date, format = '%Y-%m-%d')
london <- london[order(london$Date),]

london <- london[,c('Date', 'newCasesBySpecimenDate', 
                    'newDailyNsoDeathsByDeathDate', 'newAdmissions')]
colnames(london) <- c('date', 'daily_cases', 'daily_deaths', 'daily_hosp')

# population
london$Population <- 9002400

# peak identifier
london$peak <- NA

# peak 1 - Mar 1 - Jun 15, 2020
startDate <- as.Date('2020-03-01')
endDate <- as.Date('2020-06-15')
london$peak[london$date >= startDate & london$date < endDate]<- 1

# peak 2 - Sep 15, 2020 - Apr 1, 2021
startDate <- as.Date('2020-09-15')
endDate <- as.Date('2021-04-01')
london$peak[london$date >= startDate & london$date < endDate]<- 2


# peak 3 - Jul 1 - Nov 1, 2021
startDate <- as.Date('2021-05-01')
endDate <- as.Date('2021-10-01')
london$peak[london$date >= startDate & london$date < endDate]<- 3


plot(london$date[which(london$peak == 3)], 
     london$daily_cases[which(london$peak == 3)], pch = 16)

# peak 4 - Nov 11, 2021 - Mar 1, 2022
startDate <- as.Date('2021-11-15')
endDate <- as.Date('2022-03-01')
london$peak[london$date >= startDate & london$date < endDate]<- 4

plot(london$date[which(london$peak == 4)], 
     london$daily_cases[which(london$peak == 4)], pch = 16, type = 'l')

plot(london$date, london$daily_cases, pch = 16)

write.csv(london, 'londonClean.csv', quote = F, row.names = F)
