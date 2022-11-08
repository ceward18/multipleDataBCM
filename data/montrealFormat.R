################################################################################
# montreal data formatting
################################################################################


library(lubridate)

### Read in data from Montreal website
montreal <-read.csv2("https://santemontreal.qc.ca/fileadmin/fichiers/Campagnes/coronavirus/situation-montreal/courbe.csv") 
montreal <- montreal[!is.na(montreal$Nouveaux.cas),]

montreal <- montreal[,c('Date', 'Nouveaux.cas')]
colnames(montreal) <- c('date', 'dailyCases')

# format dates
montreal$date <- as.Date(montreal$date)

# population
montreal$Population <- 2065647

# peak identifier
montreal$peak <- NA

# peak 1 - Feb 25 - July, 2020
startDate <- as.Date('2020-02-25')
endDate <- as.Date('2020-07-08')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 1

# peak 2 - Sep 2020 - July 2021
startDate <- as.Date('2020-09-01')
endDate <- as.Date('2021-07-01')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 2

# peak 3 - Jul 1 - Nov 1, 2021
startDate <- as.Date('2021-07-01')
endDate <- as.Date('2021-11-01')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 3

# peak 4 - Nov 1, 2021 - Mar 15, 2022
startDate <- as.Date('2021-11-01')
endDate <- as.Date('2022-03-15')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 4

# peak 5 -Mar 15, 2022 - Now
startDate <- as.Date('2022-03-16')
endDate <- as.Date('2022-05-13')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 5


write.csv(montreal, 'montrealClean.csv', quote = F, row.names = F)


