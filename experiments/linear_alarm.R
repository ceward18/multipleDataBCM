################################################################################
# Linear alarm function with three covariates
################################################################################

library(nimble)
library(ggplot2)
library(lattice)

# nimbleFunction for logistic decay
logitDecay <- nimbleFunction(     
    run = function(x = double(1), w0 = double(0), k = double(0)) {
        returnType(double(1))
        
        result <- 1 / (1 + exp(k * (x - w0)))
        
        return(result)
    })

# calculate moving average for smoothing
movingAverage <- nimbleFunction(     
    run = function(x = double(1), bw = double(0)) {
        returnType(double(1))
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
    })

# get effective reproductive number at time t using a forward sum
get_R0 <- nimbleFunction(     
    run = function(betat = double(1), N = double(0), S = double(1),
                   maxInf = double(0), iddCurve = double(1)) {
        returnType(double(1))
        
        # probability of transition given 1 infectious individual (vector length t)
        pi_SI <- 1 - exp(- betat / N )
        
        nTime <- length(pi_SI)
        bw <- maxInf
        sumSmooth <- rep(NA, nTime - bw)
        for(k in 1:(nTime - bw)){
            t1 <- k
            t2 <- k + bw - 1
            pi_SI <- 1 - exp(-betat[t1:t2] * iddCurve / N)
            sumSmooth[k] <- sum(pi_SI) * S[k]
        }
        
        return(sumSmooth)
    })


# linear alarm
linearAlarm <- nimbleFunction(     
    run = function(x1 = double(1), x2 = double(1), x3 = double(1),
                   alpha = double(1)) {
        returnType(double(1))
        
        logit(result) <- alpha[1] * x1 + alpha[2] * x2 + alpha[3] * x3
        
        return(result)
    })


# Hill-Langmuir alarm function with fixed max at 1
hillAlarm <- nimbleFunction(     
    run = function(x = double(0), nu = double(0), x0 = double(0)) {
        returnType(double(0))
        
        x <- alpha[1] * x1 + alpha[2] * x2 + alpha[3] * x3
        
        result <- 1 / (1 + (x0 / x) ^ nu)
        
        return(result)
    })

SIR_linear_multi <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1,1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1,2:maxInf] <- 0
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # calculate alarm function
        # linear alarm function
        alarm[t] <- alpha[1] * smoothC[t] + 
            alpha[2] * smoothH[t] + 
            alpha[3] * smoothD[t]
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0_update[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                        N = N, S = S[1:tau], maxInf = maxInf,
                                        iddCurve = idd_curve[1:maxInf])
    
    for (i in 1:n) {
        yAlarm[i] <-  alpha[1] * xC[i] + 
            alpha[2] * xH[i] + 
            alpha[3] * xD[i]
    }
    
    # constrain alarm to be between 0 and 1
    minYAlarm <- min(alarm[1:tau])
    maxYAlarm <- max(alarm[1:tau])
    constrain_min ~ dconstraint(minYAlarm >= 0)
    constrain_max ~ dconstraint(maxYAlarm <= 1)
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    alpha[1] ~ dgamma(0.1, 0.1)
    alpha[2] ~ dgamma(0.1, 0.1)
    alpha[3] ~ dgamma(0.1, 0.1)
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
    
})


# read data
dat <- read.csv('./Data/nycClean.csv')

dat$smoothedCases <- round(dat$dailyCases)

# constants for all models
N <- dat$Population[1]
lengthI <- 3

peak <- 4
smoothWindow <- 60

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

Rstar0 <- dat$smoothedCases[max(1, idxStart - lengthI + 1):(idxStart)]

# constants that are the same for all models
S0 <- N - I0 - R0
tau <- length(incData)

### initial conditions probability
initProb <- c(S0, I0, N - S0 - I0)/N
SIR_init <- rmulti(1, N, initProb)

### constants
n <- 10

# rescale alarm components so all on the same scale (0 to 1)
smoothC <- smoothC / max(smoothC)
smoothH <- smoothH / max(smoothH)
smoothD <- smoothD / max(smoothD)

maxC <- ceiling(max(smoothC))
maxH <- ceiling(max(smoothH))
maxD <- ceiling(max(smoothD))

xC <- seq(0, maxC, length.out = n)
xH <- seq(0, maxH, length.out = n)
xD <- seq(0, maxD, length.out = n)

xGrid <- expand.grid(xC, xH, xD)

constantsList <- list(tau = tau,
                      N = N,
                      initProb = initProb,
                      n = nrow(xGrid),
                      xC = xGrid[,1],
                      xH = xGrid[,2],
                      xD = xGrid[,3],
                      maxInf = 10)

### data
dataList <- list(Istar = incData,
                 smoothC = smoothC,
                 smoothH = smoothH,
                 smoothD = smoothD,
                 constrain_min = 1,
                 constrain_max = 1)

### inits 
initsList <- list(SIR_init = SIR_init,
                  beta = runif(1, 1/7, 1),
                  alpha = runif(3, 0, 0.0001),
                  w0 = rnorm(1, 2, 0.1),
                  k = rgamma(1, 100, 100))




### create nimble model
myModel <- nimbleModel(SIR_linear_multi, 
                       data = dataList, 
                       constants = constantsList,
                       inits = initsList)

myConfig <- configureMCMC(myModel)

# need to ensure all stochastic nodes are monitored for WAIC calculation
myConfig$addMonitors(c('yAlarm', 'alarm', 'R0_update'))


myMCMC <- buildMCMC(myConfig)
compiled <- compileNimble(myModel, myMCMC) 

samples <- runMCMC(compiled$myMCMC, 
                   niter = 200000,
                   nburnin = 100000,
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


par(mfrow = c(1,2))
plot(samples[,'SIR_init[1]'], type = 'l')
plot(samples[,'SIR_init[2]'], type = 'l')


##########################################################################

alarmSamples <- samples[,grep('R0', colnames(samples))]

toPlot <- data.frame(x = 1:(tau - 10),
                     mean = colMeans(alarmSamples),
                     lower = apply(alarmSamples, 2, quantile, probs = 0.025),
                     upper = apply(alarmSamples, 2, quantile, probs = 0.975))

ggplot(toPlot, aes(x = x, y= mean)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    theme_bw()

##########################################################################

alarmSamples <- samples[,grep('alarm', colnames(samples))]

toPlot <- data.frame(x = 1:tau,
                     mean = colMeans(alarmSamples),
                     lower = apply(alarmSamples, 2, quantile, probs = 0.025),
                     upper = apply(alarmSamples, 2, quantile, probs = 0.975))

ggplot(toPlot, aes(x = x, y= mean)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    theme_bw()

##########################################################################

yAlarmSamples <- samples[,grep('yAlarm', colnames(samples))]

toPlot <- data.frame(xC = xGrid[,1],
                     xH = xGrid[,2],
                     xD = xGrid[,3],
                     mean = colMeans(yAlarmSamples),
                     lower = apply(yAlarmSamples, 2, quantile, probs = 0.025),
                     upper = apply(yAlarmSamples, 2, quantile, probs = 0.975))


wireframe(mean ~ xC * xH, data=toPlot, shade=TRUE, scales=list(arrows=FALSE))
wireframe(mean ~ xC * xD, data=toPlot, shade=TRUE, scales=list(arrows=FALSE))
wireframe(mean ~ xH * xD, data=toPlot, shade=TRUE, scales=list(arrows=FALSE))
