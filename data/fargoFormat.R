################################################################################
# Fargo, ND data formatting
################################################################################

cases <- read.csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv')
deaths <- read.csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv')
hosp <- read.csv('https://healthdata.gov/resource/g62h-syeh.csv')


# just ND
cases <- cases[cases$Province_State == 'North Dakota' & cases$Admin2 == 'Cass',]
deaths <- deaths[deaths$Province_State == 'North Dakota' & deaths$Admin2 == 'Cass',]
hosp <- hosp[hosp$state == 'MN',]


# wide to long
dateCols <- regmatches(colnames(cases), regexpr('X.*', colnames(cases))) 
casesLong <- reshape(cases, 
                   varying = dateCols, 
                   v.names = "cumulative_cases",
                   timevar = "date", 
                   times = dateCols, 
                   new.row.names = 1:300000,
                   direction = "long")


dateCols <- regmatches(colnames(deaths), regexpr('X.*', colnames(deaths))) 
deathsLong <- reshape(deaths, 
                     varying = dateCols, 
                     v.names = "cumulative_cases",
                     timevar = "date", 
                     times = dateCols, 
                     new.row.names = 1:300000,
                     direction = "long")

hosp <- hosp[,c('date', 'hospital_onset_covid')]
hosp <- hosp[order(hosp$date),]

hosp$date <- as.Date(substr(hosp$date, 1, 10), format = "%Y-%m-%d")
plot(hosp$date, hosp$hospital_onset_covid)
