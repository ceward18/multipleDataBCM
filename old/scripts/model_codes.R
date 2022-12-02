################################################################################
# Script containing model code and associated nimbleFunctions
# to be sourced by model fitting script
################################################################################


################################################################################
### Helper functions

# nimbleFunction for logistic decay IDD curve
logitDecay <- nimbleFunction(     
    run = function(x = double(1), w0 = double(0), k = double(0)) {
        returnType(double(1))
        
        result <- 1 / (1 + exp(k * (x - w0)))
        
        return(result)
    })
assign('logitDecay', logitDecay, envir = .GlobalEnv)


# multivariate threshold alarm function
multiThresholdAlarm <- nimbleFunction(     
    run = function(x1 = double(0), x2 = double(0), x3 = double(0), N = double(0), 
                   delta1 = double(0), delta2 = double(0), delta3 = double(0), 
                   H1 = double(0), H2 = double(0), H3 = double(0)) {
        returnType(double(0))
        
        result <- delta1 * (x1 / N > H1) +
            delta2 * (x2 / N > H2) + 
            delta3 * (x3 / N > H3)
        
        return(result)
    })
assign('multiThresholdAlarm', multiThresholdAlarm, envir = .GlobalEnv)

# Hill-Langmuir alarm function
hillAlarm <- nimbleFunction(     
    run = function(x = double(0), nu = double(0), x0 = double(0), 
                   delta = double(0)) {
        returnType(double(0))
        
        result <- delta / (1 + (x0 / x) ^ nu)
        
        return(result)
    })
assign('hillAlarm', hillAlarm, envir = .GlobalEnv)


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

# squared exponential covariance for Gaussian process
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
assign('sqExpCov', sqExpCov, envir = .GlobalEnv)

# spline alarm define first as an R function
splineAlarmR <- function(xC, xH, xD, b, knotsC, knotsH, knotsD) {
    xBasis <- cbind(splines::ns(xC, knots = knotsC),
                    splines::ns(xH, knots = knotsH),
                    splines::ns(xD, knots = knotsD))
    c(xBasis %*% b)
    
}

# then convert to something that can be compiled in nimble
splineAlarm <- nimbleRcall(function(xC = double(1), xH = double(1), xD = double(1),
                                    b = double(1), 
                                    knotsC = double(1), knotsH = double(1), knotsD = double(1)){}, 
                           Rfun = 'splineAlarmR',
                           returnType = double(1)
)
assign('splineAlarm', splineAlarm, envir = .GlobalEnv)


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
assign('nim_approx', nim_approx, envir = .GlobalEnv)

nim_approx_tri <- nimbleFunction(     
    run = function(x = double(1), y = double(1), z = double(1),
                   grid = double(2), f = double(1), 
                   xout = double(0), yout = double(0), zout = double(0)) {
        returnType(double(0))
        
        # if xout is > max(x), return the closest value
        if (xout >= max(x)) xout <- max(x) - 1e-5
        if (yout >= max(y)) yout <- max(y) - 1e-5
        if (zout >= max(z)) zout <- max(z) - 1e-5
        
        # x values on either side of xout
        xPlacement <- 1 * (xout < x)
        yPlacement <- 1 * (yout < y)
        zPlacement <- 1 * (zout < z)
        
        x0 <- x[max(which(xPlacement == 0))]
        x1 <- x[min(which(xPlacement == 1))]
        
        y0 <- y[max(which(yPlacement == 0))]
        y1 <- y[min(which(yPlacement == 1))]
        
        z0 <- z[max(which(zPlacement == 0))]
        z1 <- z[min(which(zPlacement == 1))]
        
        xd <- (xout - x0) / (x1 - x0)
        yd <- (yout - y0) / (y1 - y0)
        zd <- (zout - z0) / (z1 - z0)
        
        # function values across cube
        c000 <- f[which(grid[,1] == x0 & grid[,2] == y0 & grid[,3] == z0)]
        c100 <- f[which(grid[,1] == x1 & grid[,2] == y0 & grid[,3] == z0)]
        c001 <- f[which(grid[,1] == x0 & grid[,2] == y0 & grid[,3] == z1)]
        c101 <- f[which(grid[,1] == x1 & grid[,2] == y0 & grid[,3] == z1)]
        c010 <- f[which(grid[,1] == x0 & grid[,2] == y1 & grid[,3] == z0)]
        c110 <- f[which(grid[,1] == x1 & grid[,2] == y1 & grid[,3] == z0)]
        c011 <- f[which(grid[,1] == x0 & grid[,2] == y1 & grid[,3] == z1)]
        c111 <- f[which(grid[,1] == x1 & grid[,2] == y1 & grid[,3] == z1)]
        
        # interpolate along x
        c00 <- c000 * (1 - xd) + c100 * xd
        c01 <- c001 * (1 - xd) + c101 * xd
        c10 <- c010 * (1 - xd) + c110 * xd
        c11 <- c011 * (1 - xd) + c111 * xd
        
        # interpolate along y
        c0 <- c00 * (1 - yd) + c10 * yd
        c1 <- c01 * (1 - yd) + c11 * yd
        
        # interpolate along z
        out <- c0 * (1 - zd) + c1 * zd
        
        return(out[1])
    })
assign('nim_approx_tri', nim_approx_tri, envir = .GlobalEnv)


# function to simulate from myModel
simulator <- nimbleFunction(
    setup = function(model, dataNodes) {
        parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
        # exclude data from parent nodes
        parentNodes <- parentNodes[-which(parentNodes %in% dataNodes)]
        parentNodes <- model$expandNodeNames(parentNodes, returnScalarComponents = TRUE)
        cat("Stochastic parents of data are: ", paste(parentNodes, sep = ','), ".\n")
        simNodes <- model$getDependencies(parentNodes, self = FALSE,
                                          downstream = T)
        
        nData <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
    },
    run = function(params = double(1), nSim = double()) {
        simDat <- matrix(nrow = nSim, ncol = nData)   
        for(i in 1:nSim) {
            values(model, parentNodes) <<- params
            model$simulate(simNodes, includeData = TRUE)
            simDat[i, ] <- values(model, dataNodes)
        }
        return(simDat)
        returnType(double(2))
    })
assign('simulator', simulator, envir = .GlobalEnv)

# determine mean for prior on lengthscale parameter
getl <- function(maxDist) {
    sqrt(- (maxDist/2) ^2 / (2 * log(0.025)))
}

# determine invgamma shape and scale for prior on lengthscale parameter
myF <- function(x, mid) {
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
### Univariate GP Model (alarm based only on one thing)

SIR_gp_uni <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1,1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1,2:maxInf] <- 0
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], alarmBasis[t])
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    
    yAlarm[1] <- 0
    mu[1:n] <- mu0 * ones[1:n]
    
    # cases
    cov[2:n, 2:n] <- sqExpCov(dists[2:n, 2:n], sigma, l)
    logit(yAlarm[2:n]) ~ dmnorm(mu[2:n], cov = cov[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(100, 50)
    l ~ dinvgamma(lShape, lScale)
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
})

################################################################################
### Multivariate Threshold Alarm

SIR_thresh_multi <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1,1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1,2:maxInf] <- 0
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        #  sum of each component
        alarm[t] <- multiThresholdAlarm(smoothC[t], smoothH[t], smoothD[t], N,
                                        delta1, delta2, delta3,
                                        H1, H2, H3)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    # compute alarm over input grid
    for (i in 1:nn) {
        yAlarm[i] <- multiThresholdAlarm(xC_grid[i], xH_grid[i], xD_grid[i], N,
                                         delta1, delta2, delta3,
                                         H1, H2, H3)
    }
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    delta1 ~ dbeta(1, 1)
    delta2 ~ dbeta(1, 1)
    delta3 ~ dbeta(1, 1)
    H1 ~ dunif(minC/N, maxC/N)
    H2 ~ dunif(minH/N, maxH/N)
    H3 ~ dunif(minD/N, maxD/N)
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
    
    # constrain deltas to sum to 1
    constrain_deltas ~ dconstraint(delta1 + delta2 + delta3 <= 1)
    
})



################################################################################
### Relative GP Model (alarm based only on convex combo of three things)
# correlation accounted for by model?

SIR_gp_rel <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1,1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1,2:maxInf] <- 0
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm- each piece is between 0 and 1
        yAlarmC[t] <- nim_approx(xC[1:n], yC[1:n], smoothC[t])
        yAlarmH[t] <- nim_approx(xH[1:n], yH[1:n], smoothH[t])
        yAlarmD[t] <- nim_approx(xD[1:n], yD[1:n], smoothD[t])
        
        # weighted sum of each component
        alarm[t] <- alpha[1] * yAlarmC[t] +
            alpha[2] * yAlarmH[t] + 
            alpha[3] * yAlarmD[t]
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    
    yC[1] <- 0
    yH[1] <- 0
    yD[1] <- 0
    mu[1:n] <- mu0 * ones[1:n]
    
    # cases
    covC[2:n, 2:n] <- sqExpCov(distC[2:n, 2:n], sigma, lC)
    logit(yC[2:n]) ~ dmnorm(mu[2:n], cov = covC[2:n, 2:n])
    
    # hospitalizations
    covH[2:n, 2:n] <- sqExpCov(distH[2:n, 2:n], sigma, lH)
    logit(yH[2:n]) ~ dmnorm(mu[2:n], cov = covH[2:n, 2:n])
    
    # deaths
    covD[2:n, 2:n] <- sqExpCov(distD[2:n, 2:n], sigma, lD)
    logit(yD[2:n]) ~ dmnorm(mu[2:n], cov = covD[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(100, 50)
    lC ~ dinvgamma(lCShape, lCScale)
    lH ~ dinvgamma(lHShape, lHScale)
    lD ~ dinvgamma(lDShape, lDScale)
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
    
    # relative importance
    alpha[1:3] ~ ddirch(priorWeights[1:3])
})

################################################################################
### Hierarchical Hill Model 

SIR_gp_hier <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1,1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1,2:maxInf] <- 0
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm- each piece is between 0 and 1
        yAlarmC[t] <- nim_approx(xC[1:n], yC[1:n], smoothC[t])
        yAlarmH[t] <- nim_approx(xH[1:n], yH[1:n], smoothH[t])
        yAlarmD[t] <- nim_approx(xD[1:n], yD[1:n], smoothD[t])
        
        # weighted sum of each component
        alarm[t] <- alpha[1] * yAlarmC[t] +
            alpha[2] * yAlarmH[t] + 
            alpha[3] * yAlarmD[t]
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    
    yC[1] <- 0
    yH[1] <- 0
    yD[1] <- 0
    mu[1:n] <- mu0 * ones[1:n]
    
    # cases
    covC[2:n, 2:n] <- sqExpCov(distC[2:n, 2:n], sigma, lC)
    logit(yC[2:n]) ~ dmnorm(mu[2:n], cov = covC[2:n, 2:n])
    
    # hospitalizations
    covH[2:n, 2:n] <- sqExpCov(distH[2:n, 2:n], sigma, lH)
    logit(yH[2:n]) ~ dmnorm(mu[2:n], cov = covH[2:n, 2:n])
    
    # deaths
    covD[2:n, 2:n] <- sqExpCov(distD[2:n, 2:n], sigma, lD)
    logit(yD[2:n]) ~ dmnorm(mu[2:n], cov = covD[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(100, 50)
    lC ~ dinvgamma(lCShape, lCScale)
    lH ~ dinvgamma(lHShape, lHScale)
    lD ~ dinvgamma(lDShape, lDScale)
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
    
    # relative importance
    alpha[1:3] ~ ddirch(priorWeights[1:3])
})


################################################################################
### Multivariate GP Model (alarm based only on one of three things)

# RJ MCMC proposal for z[1:3]
zUpdate <- nimbleFunction(
    name = 'zUpdate',                              
    contains = sampler_BASE,                     
    setup = function(model, mvSaved, target, control) {                 # REQUIRED setup arguments
        calcNodes <- model$getDependencies(target) 
        
    },                                                                  # setup can't return anything
    run = function() {
        currentValue <- model[[target]]                                   
        currentLogProb <- model$getLogProb(calcNodes)    
        
        # independent proposal
        proposalValue <- rmulti(n = 1, size = 1, prob = rep(1, 3))
        
        # put proposal value in model
        model[[target]] <<- proposalValue                                
        proposalLogProb <- model$calculate(calcNodes)                     
        logAcceptanceRatio <- proposalLogProb - currentLogProb     
        
        accept <- decide(logAcceptanceRatio)       
        
        if(accept)   {
            copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
            copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }                                                            
        
    },
    methods = list(   # required method for sampler_BASE base class
        reset = function() {}
    )
)

assign('zUpdate', zUpdate, envir = .GlobalEnv)

SIR_gp_model_switch <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1,1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1,2:maxInf] <- 0
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) { 
        
        # compute alarm- each piece is between 0 and 1
        alarmC[t] <- nim_approx(xC[1:n], yC[1:n], smoothC[t])
        alarmH[t] <- nim_approx(xH[1:n], yH[1:n], smoothH[t])
        alarmD[t] <- nim_approx(xD[1:n], yD[1:n], smoothD[t])
        
        # indicators for each component
        alarm[t] <- z[1] * alarmC[t] + z[2] * alarmH[t] + z[3] * alarmD[t]
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
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
    lC ~ dinvgamma(lCShape, lCScale)
    lH ~ dinvgamma(lHShape, lHScale)
    lD ~ dinvgamma(lDShape, lDScale)
    
    # indicator variable for model to be used (cases, hosp, deaths)
    z[1:3] ~ dmulti(prob = zProb[1:3], size = 1)
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
})



################################################################################
### Multivariate GP Model (alarm non-linear function of all three)

SIR_spline_multi <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1,1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1,2:maxInf] <- 0
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) { 
        
        # tri-linear interpolation
        alarm[t] <- nim_approx_tri(x = xC[1:n], y = xH[1:n], z = xD[1:n],
                                   grid = grid[1:nn, 1:3], 
                                   f = yAlarm[1:nn],
                                   xout = smoothC[t], yout = smoothH[t], zout = smoothD[t])
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    yAlarm[1:nn] <- splineAlarm(xC_grid[1:nn], 
                                xH_grid[1:nn], 
                                xD_grid[1:nn],
                                b[1:nb], 
                                knotsC[1:2], 
                                knotsH[1:2], 
                                knotsD[1:2])
    
    # constrain yAlarm to be between 0 and 1
    minYAlarm <- min(yAlarm[1:nn])
    maxYAlarm <- max(yAlarm[1:nn])
    constrain_min ~ dconstraint(minYAlarm >= 0)
    constrain_max ~ dconstraint(maxYAlarm <= 1)
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    
    for (i in 1:nb) {
        b[i] ~ dnorm(0, sd = 100)
    }
    for (i in 1:2) {
        knotsC[i] ~ dunif(min = minC, max = maxC)
        knotsH[i] ~ dunif(min = minH, max = maxH)
        knotsD[i] ~ dunif(min = minD, max = maxD)
    }
    
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
})


