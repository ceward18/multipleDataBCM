################################################################################
# NYC data analysis comparing incidence, deaths, hospitalizations
################################################################################


### load libraries
library(parallel)
library(nimble)
library(ggplot2)

# source scripts (for movingAverage function)
source('./scripts/modelCodes.R')

### read data
dat <- read.csv('./Data/nycClean.csv')

dat$smoothedCases <- round(dat$dailyCases)

# constants for all models
N <- dat$Population[1]
lengthI <- 2

peak <- 1
smoothWindow <- 14

# smoothed incidence/hosp/deaths to inform alarm function 
# (shifted so alarm is informed only by data up to time t-1)
dat$smoothC <- head(movingAverage(c(0, dat$dailyCases), smoothWindow), -1)
dat$smoothH <- head(movingAverage(c(0, dat$dailyHosp), smoothWindow), -1)
dat$smoothD <- head(movingAverage(c(0, dat$dailyDeaths), smoothWindow), -1)

incData <- dat$smoothedCases[which(dat$peak == peak)]
smoothC <- dat$smoothC[which(dat$peak == peak)]
smoothH <- dat$smoothH[which(dat$peak == peak)]
smoothD <- dat$smoothD[which(dat$peak == peak)]

if (peak == 1) {
    idxStart <- 5
    incData <- incData[-c(1:idxStart)]
    smoothC <- smoothC[-c(1:idxStart)]
    smoothH <- smoothH[-c(1:idxStart)]
    smoothD <- smoothD[-c(1:idxStart)]
} else {
    idxStart <- min(which(dat$peak == peak))
    incData <- incData[-1]
    smoothC <- smoothC[-1]
    smoothH <- smoothH[-1]
    smoothD <- smoothD[-1]
}

# currently infectious/removed
I0 <- sum(dat$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)])
R0 <- cumsum(dat$smoothedCases)[idxStart] - I0 

Rstar0 <- dat$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)]


# get appropriate model code
modelCode <- get('SIR_hill_multi')

# constants that are the same for all models
S0 <- N - I0 - R0
tau <- length(incData)

### initial conditions probability
initProb <- c(S0, I0, N - S0 - I0)/N
SIR_init <- rmulti(1, N, initProb)

### constants
n <- 50

minC <- floor(min(smoothC))
minH <- floor(min(smoothH))
minD <- floor(min(smoothD))

maxC <- ceiling(max(smoothC))
maxH <- ceiling(max(smoothH))
maxD <- ceiling(max(smoothD))

xC <- seq(0, maxC, length.out = n)
xH <- seq(0, maxH, length.out = n)
xD <- seq(0, maxD, length.out = n)

constantsList <- list(tau = tau,
                      N = N,
                      initProb = initProb,
                      n = n,
                      xC = xC,
                      xH = xH,
                      xD = xD,
                      minC = minC,
                      maxC = maxC,
                      minH = minH,
                      maxH = maxH,
                      minD = minD,
                      maxD = maxD,
                      maxInf = 10,
                      priorWeights = rep(1, 3))

### data
dataList <- list(Istar = incData,
                 smoothC= smoothC,
                 smoothH= smoothH,
                 smoothD= smoothD)

### inits 
initsList <- list(SIR_init = SIR_init,
                  beta = runif(1, 1/7, 1),
                  alpha = rep(1/3, 3),
                  nuC = runif(1, 5, 50),
                  x0C = max(rnorm(1, maxC/2, 10), 1),
                  nuH = runif(1, 5, 50),
                  x0H = max(rnorm(1, maxH/2, 10), 1),
                  nuD = runif(1, 5, 50),
                  x0D = max(rnorm(1, maxD/2, 10), 1),
                  w0 = rnorm(1, 2, 0.1),
                  k = rgamma(1, 100, 100))



### create nimble model
myModel <- nimbleModel(modelCode, 
                       data = dataList, 
                       constants = constantsList,
                       inits = initsList)

myConfig <- configureMCMC(myModel)

# need to ensure all stochastic nodes are monitored for WAIC calculation
myConfig$addMonitors(c('yC', 'yH', 'yD', 
                       'alarmC', 'alarmH', 'alarmD', 'alarm'))


myMCMC <- buildMCMC(myConfig)
compiled <- compileNimble(myModel, myMCMC) 

samples <- runMCMC(compiled$myMCMC, 
                   niter = 50000,
                   nburnin = 10000,
                   thin = 1,
                   setSeed  = 1)

par(mfrow = c(2,2))
plot(samples[,'beta'], type = 'l')
plot(samples[,'alpha[1]'], type = 'l')
plot(samples[,'alpha[2]'], type = 'l')
plot(samples[,'alpha[3]'], type = 'l')

par(mfrow = c(1,2))
plot(samples[,'w0'], type = 'l')
plot(samples[,'k'], type = 'l')


par(mfrow = c(2,3))
plot(samples[,'nuC'], type = 'l')
plot(samples[,'nuH'], type = 'l')
plot(samples[,'nuD'], type = 'l')
plot(samples[,'x0C'], type = 'l')
plot(samples[,'x0H'], type = 'l')
plot(samples[,'x0D'], type = 'l')


yCSamples <- samples[,grep('yC', colnames(samples))]
yHSamples <- samples[,grep('yH', colnames(samples))]
yDSamples <- samples[,grep('yD', colnames(samples))]
ySamples <- cbind.data.frame(yCSamples, yHSamples, yDSamples)

toPlot <- data.frame(x = c(xC, xH, xD),
                     type = rep(c('cases', 'hosp', 'deaths'), each = n),
                     mean = colMeans(ySamples),
                     lower = apply(ySamples, 2, quantile, probs = 0.025),
                     upper = apply(ySamples, 2, quantile, probs = 0.975))

ggplot(toPlot, aes(x = x, y= mean)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    facet_grid(~type, scales = 'free') + 
    theme_bw()




