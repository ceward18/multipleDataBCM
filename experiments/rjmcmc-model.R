################################################################################
# multivariate model using RJMCMC
################################################################################


### load libraries
library(parallel)
library(nimble)
library(ggplot2)

### read data
dat <- read.csv('./Data/nycClean.csv')

dat$smoothedCases <- round(dat$dailyCases)

################################################################################
### Helper functions

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
assign('movingAverage', movingAverage, envir = .GlobalEnv)

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
assign('get_R0', get_R0, envir = .GlobalEnv)

# squared exponential covariance for gaussian process
sqExpCov <- nimbleFunction(     
    run = function(dists = double(2), sigma = double(0), l = double(0)) {
        returnType(double(2))
        n <- dim(dists)[1]
        result <- matrix(nrow = n, ncol = n, init = FALSE)
        sigma2 <- sigma*sigma
        l2 <- 2 * l^2
        for(i in 1:n) {
            for(j in 1:n) {
                
                result[i, j] <- sigma2*exp(- dists[i,j]^2 / l2)
                
                if (i == j) {
                    result[i, j] <- result[i, j] + 1e-6
                }
            }
        }
        
        return(result)
    })

# nimbleFunction for logistic decay
logitDecay <- nimbleFunction(     
    run = function(x = double(1), w0 = double(0), k = double(0)) {
        returnType(double(1))
        
        result <- 1 / (1 + exp(k * (x - w0)))
        
        return(result)
    })


# linear interpolation function to get alarm values for each observed incidence value
nim_approx <- nimbleFunction(     
    run = function(x = double(1), y = double(1), xout = double(0)) {
        returnType(double(0))
        
        # if xout is > max(x), return the closest value
        if (xout >= max(x)) {
            return(y[which(x == max(x))[1]])
        }
        
        # x values on either side of xout
        xPlacement <- 1 * (xout < x)
        # last 0 and first 1
        leftIdx <- max(which(xPlacement == 0))
        rightIdx <- min(which(xPlacement == 1))
        
        x0 <- x[leftIdx]
        y0 <- y[leftIdx]
        
        x1 <- x[rightIdx]
        y1 <- y[rightIdx]
        
        # linear interpolation for x
        out <- y0  + (xout - x0) * ((y1 - y0)/(x1 - x0))
        
        return(out)
    })


# determine mean for prior on lengthscale parameter
getl <- function(maxDist) {
    sqrt(- (maxDist/2) ^2 / (2 * log(0.025)))
}

# determine shape and scale for prior on rho
myF <- function(x, min, mid,  max) {
    a <- x[1]; b <- x[2]
    
    # mean is b/(a - 1)
    distMean <-  b/(a - 1) 
    
    # variance is b^2 / ((a-1)^2 * (a-2))
    distSD <- sqrt(b^2 / ((a-1)^2 * (a-2)))
    sdEst <- 4
    
    summat <- c(distMean - mid,
                distSD - sdEst)
    sum(summat^2)
}


################################################################################
### Model Code

SIR_gp_multi <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1,1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1,2:maxInf] <- 0
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) { 
        
        # compute alarm- each piece is between 0 and 1
        fCases[t] <- nim_approx(xC[1:n], yC[1:n], smoothC[t])
        fHosp[t] <- nim_approx(xH[1:n], yH[1:n], smoothH[t])
        fDeath[t] <- nim_approx(xD[1:n], yD[1:n], smoothD[t])
        
        # indicators for each component
        alarm[t] <- z[1] * fCases[t] + z[2] * fHosp[t] + z[3] * fDeath[t]
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    yC[1] <- 0
    yH[1] <- 0
    yD[1] <- 0
    mu[1:n] <- mu0 * ones[1:n]
    
    # cases - GP
    covC[1:n, 1:n] <- sqExpCov(distC[1:n, 1:n], sigma, lC)
    logit(yC[2:n]) ~ dmnorm(mu[2:n], cov = covC[2:n, 2:n])
    
    # hospitalizations - GP
    covH[1:n, 1:n] <- sqExpCov(distH[1:n, 1:n], sigma, lH)
    logit(yH[2:n]) ~ dmnorm(mu[2:n], cov = covH[2:n, 2:n])
    
    # deaths - GP
    covD[1:n, 1:n] <- sqExpCov(distD[1:n, 1:n], sigma, lD)
    logit(yD[2:n]) ~ dmnorm(mu[2:n], cov = covD[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(100, 50)
    lC ~ dinvgamma(cc, dc)
    lH ~ dinvgamma(ch, dh)
    lD ~ dinvgamma(cd, dd)
    
    # indicator variable for model to be used (cases, hosps, deaths)
    z[1:3] ~ dmulti(size = 1, prob = zProb[1:3])
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
})


################################################################################
### Model Set up


# constants for all models
N <- dat$Population[1]
lengthI <- 2

peak <- 1
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

Rstar0 <- dat$smoothedCases[max(1, (idxStart - lengthI + 1)):(idxStart)]

# constants that are the same for all models
S0 <- N - I0 - R0
tau <- length(incData)

### initial conditions probability
initProb <- c(S0, I0, N - S0 - I0)/N
SIR_init <- rmulti(1, N, initProb)

### constants
n <- 10

maxC <- ceiling(max(smoothC))
maxH <- ceiling(max(smoothH))
maxD <- ceiling(max(smoothD))

xC <- seq(0, maxC, length.out = n)
xH <- seq(0, maxH, length.out = n)
xD <- seq(0, maxD, length.out = n)


distC <- as.matrix(dist(matrix(xC)))
distH <- as.matrix(dist(matrix(xH)))
distD <- as.matrix(dist(matrix(xD)))

midDist <- getl(max(distC[lower.tri(distC)]))
valsC <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                     mid = midDist)$par, 2)

midDist <- getl(max(distH[lower.tri(distH)]))
valsH <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                     mid = midDist)$par, 2)

midDist <- getl(max(distD[lower.tri(distD)]))
valsD <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                     mid = midDist)$par, 2)

constantsList <- list(tau = tau,
                      N = N,
                      initProb = initProb,
                      n = n,
                      distC = distC,
                      distH = distH,
                      distD = distD,
                      mu0 = 1,
                      ones = logit(seq(0.0001, 0.9999, length.out= n)),
                      xC = xC,
                      xH = xH,
                      xD = xD,
                      maxInf = 10,
                      priorWeights = rep(1, 3),
                      cc = valsC[1],
                      dc = valsC[2],
                      ch = valsH[1],
                      dh = valsH[2],
                      cd = valsD[1],
                      dd = valsD[2],
                      zProb = rep(1/3, 3))

### data
dataList <- list(Istar = incData,
                 smoothC= smoothC,
                 smoothH= smoothH,
                 smoothD= smoothD)

### inits 
initsList <- list(SIR_init = SIR_init,
                  beta = runif(1, 1/7, 1),
                  z = c(0, 1, 0),
                  lC = rinvgamma(1, valsC[1], valsC[2]),
                  lH = rinvgamma(1, valsH[1], valsH[2]),
                  lD = rinvgamma(1, valsD[1], valsD[2]),
                  sigma = rgamma(1, 100, 50),
                  w0 = rnorm(1, 2, 0.1),
                  k = rgamma(1, 100, 100))




### create nimble model
myModel <- nimbleModel(SIR_gp_multi, 
                       data = dataList, 
                       constants = constantsList,
                       inits = initsList)

myConfig <- configureMCMC(myModel)

# need to ensure all stochastic nodes are monitored for WAIC calculation
myConfig$addMonitors(c('yC', 'yH', 'yD', 'alarm'))

# RJ MCMC
# 
# configureRJ(myConfig,
#             targetNodes = c('fCases', 'fHosp', 'fDeath'),
#             indicatorNodes = 'z')


myMCMC <- buildMCMC(myConfig)
compiled <- compileNimble(myModel, myMCMC) 

samples <- runMCMC(compiled$myMCMC, 
                   niter = 100,
                   nburnin = 0,
                   thin = 1,
                   setSeed  = 1)

par(mfrow = c(2,2))
plot(samples[,'beta'], type = 'l')
plot(samples[,'z[1]'], type = 'l')
plot(samples[,'z[2]'], type = 'l')
plot(samples[,'z[3]'], type = 'l')

par(mfrow = c(1,2))
plot(samples[,'w0'], type = 'l')
plot(samples[,'k'], type = 'l')


par(mfrow = c(2,2))
plot(samples[,'lC'], type = 'l')
plot(samples[,'lH'], type = 'l')
plot(samples[,'lD'], type = 'l')
plot(samples[,'sigma'], type = 'l')



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



