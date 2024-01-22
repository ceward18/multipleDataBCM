################################################################################
#### estimate priors for hospitalization and deaths rates
# gamma1 ~ dgamma(2, 10) # IR 
# gamma2 ~ dgamma(2, 10) # HR 
# lambda ~ dgamma(0.5, 10) # IH 
# phi ~ dgamma(1, 10)    # HD 
################################################################################

dat <- read.csv('./data/montrealClean.csv')

### peak 1 (lambda = 0.05, phi = 0.17)
peakDat <- dat[which(dat$peak == 1),]

# hospitalizations/deaths occur X days after (lambda = 0.45, phi = 0.8)
# compare peaks
which.max(peakDat$dailyCases)
which.max(peakDat$dailyHosp)
which.max(peakDat$dailyDeaths)

summary(tail(peakDat$dailyHosp, -3) / head(peakDat$dailyCases, -3))
summary(tail(peakDat$dailyDeaths, -10) / head(peakDat$dailyHosp, -10))

### peak 2 
peakDat <- dat[which(dat$peak == 2),]

# hospitalizations/deaths occur X days after (lambda = 0.05, phi = 0.2)
# compare peaks
which.max(peakDat$dailyCases)
which.max(peakDat$dailyHosp)
which.max(peakDat$dailyDeaths)

summary(tail(peakDat$dailyHosp, -3) / head(peakDat$dailyCases, -3))
summary(tail(peakDat$dailyDeaths, -12) / head(peakDat$dailyHosp, -12))

### peak 4

peakDat <- dat[which(dat$peak == 4),]

# hospitalizations/deaths occur X days after (lambda = 0.05, phi = 0.1)
# compare peaks
which.max(peakDat$dailyCases)
which.max(peakDat$dailyHosp)
which.max(peakDat$dailyDeaths)

summary(tail(peakDat$dailyHosp, -10) / head(peakDat$dailyCases, -10))
summary(tail(peakDat$dailyDeaths, -6) / head(peakDat$dailyHosp, -6))
