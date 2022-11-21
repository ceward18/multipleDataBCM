################################################################################
# montreal data formatting
# data downloaded from :
# https://www.inspq.qc.ca/covid-19/donnees on Nov 21, 2022
# waves defined on their website:
# wave 1: Feb 25 - 11 July 2020
# wave 2: Aug 23, 2020 - March 20, 2021
# wave 3: March 21 - July 17, 2021
# wave 4: July 18 - Dec 4, 2021
# wave 5: Dec 5, 2021 - March 12, 2022
# wave 6: March 13 - May 28, 2022
################################################################################


library(lubridate)

### Read in downloaded data
cases <- read.csv("graph_1-1_page_par_region.csv")[,c(1,2)]
hosps <- read.csv("graph_3-2_page_par_region.csv")
deaths <- read.csv("graph_2-2_page_par_region.csv")

colnames(cases) <- c('date', 'dailyCases')
colnames(hosps) <- c('date', 'non-ICU', 'ICU')
colnames(deaths) <- c('date', 'intermediate', 'home', 'seniorLiving', 'LTC')

hosps$dailyHosp <- rowSums(hosps[,-1])
deaths$dailyDeaths <- rowSums(deaths[,-1])

hosps <- hosps[,c('date', 'dailyHosp')]
deaths <- deaths[,c('date', 'dailyDeaths')]

# format dates
cases$date <- as.Date(cases$date)
hosps$date <- as.Date(hosps$date)
deaths$date <- as.Date(deaths$date)


# merge together

montreal <- merge(cases, hosps, by = 'date')
montreal <- merge(montreal, deaths, by = 'date', all.x = T)

montreal$dailyDeaths[is.na(montreal$dailyDeaths)] <- 0

# population
montreal$Population <- 1762949

# peak identifier
montreal$peak <- NA

# peak 1 - Feb 25 - 11 July 2020 (first case observed march 6)
startDate <- as.Date('2020-03-11')
endDate <- as.Date('2020-07-12')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 1


plot(montreal$date[which(montreal$peak == 1)], 
     montreal$dailyCases[which(montreal$peak == 1)], pch = 16)

# peak 2 - Aug 23, 2020 - March 20, 2021
# include "wave 3": March 21 - July 17, 2021
startDate <- as.Date('2020-08-23')
endDate <- as.Date('2021-07-18')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 2

plot(montreal$date[which(montreal$peak == 2)], 
     montreal$dailyCases[which(montreal$peak == 2)], pch = 16)

# peak 3 - July 18 - Dec 4, 2021
startDate <- as.Date('2021-07-18')
endDate <- as.Date('2021-12-05')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 3

plot(montreal$date[which(montreal$peak == 3)], 
     montreal$dailyCases[which(montreal$peak == 3)], pch = 16)

# peak 4 - Dec 5, 2021 - March 12, 2022
startDate <- as.Date('2021-12-05')
endDate <- as.Date('2022-03-13')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 4

plot(montreal$date[which(montreal$peak == 4)], 
     montreal$dailyCases[which(montreal$peak == 4)], pch = 16)

# peak 5 - March 13 - May 28, 2022
startDate <- as.Date('2022-03-13')
endDate <- as.Date('2022-05-29')
montreal$peak[montreal$date >= startDate & montreal$date < endDate]<- 5

plot(montreal$date[which(montreal$peak == 5)], 
     montreal$dailyCases[which(montreal$peak == 5)], pch = 16)

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

montreal$dailyCases <- movingAverage(montreal$dailyCases, 7)
montreal$dailyDeaths <- movingAverage(montreal$dailyDeaths, 7)
montreal$dailyHosp <- movingAverage(montreal$dailyHosp, 7)


write.csv(montreal, 'montrealClean.csv', quote = F, row.names = F)


