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
### Multivariate GP Model (alarm based only on one of three things)

# RJ MCMC proposal for z[1:3]
zUpdate <- nimbleFunction(
    name = 'z',                              
    contains = sampler_BASE,                     
    setup = function(model, mvSaved, target, control) {                 # REQUIRED setup arguments
        calcNodes <- model$getDependencies(target) 
        
        # number of update attempts is some % of the final epidemic size (total # of removals)
        nUpdates <- 200
    },                                                                  # setup can't return anything
    run = function() {
        currentValue <- model[[target]]                                   
        currentLogProb <- model$getLogProb(calcNodes)                    
        
        # repeat proposal many times 
        for (it in 1:nUpdates) {
            
            # three possible moves:
            moveType <- ceiling(runif(1, 0, 3))
            
            proposalValue <- currentValue
            
            nTimePoints <- length(currentValue)
            
            if (moveType == 1) {
                # add a removal time
                addIdx <- runif(1, 1, nTimePoints + 1)
                proposalValue[addIdx] <- proposalValue[addIdx] + 1
                
                # g(old|new) - g(new|old)
                # subtract from new - add to old
                possibleSubtract <- which(proposalValue > 0)
                g <- -log(length(possibleSubtract)) + log(nTimePoints)
                
                
            } else if (moveType == 2) {
                # move a removal time
                possibleSubtract <- which(currentValue > 0)
                subtractIdx <- possibleSubtract[runif(1, 1, length(possibleSubtract) + 1)]
                addIdx <- runif(1, 1, nTimePoints + 1)
                
                proposalValue[subtractIdx] <- proposalValue[subtractIdx] - 1
                proposalValue[addIdx] <- proposalValue[addIdx] + 1
                
                # g(old|new) - g(new|old)
                # possibly have different number of values to subtract from 
                newPossibleSubtract <- which(proposalValue > 0)
                g <- -log(length(newPossibleSubtract)) +log(length(possibleSubtract))
                
            } else if (moveType == 3) {
                # subtract a removal time
                possibleSubtract <- which(currentValue > 0)
                subtractIdx <- possibleSubtract[runif(1, 1, length(possibleSubtract) + 1)]
                proposalValue[subtractIdx] <- proposalValue[subtractIdx] - 1
                
                # g(old|new) - g(new|old)
                # add to new - subtract from old
                g <- -log(nTimePoints) + log(length(possibleSubtract)) 
                
            }
            
            # put proposal value in model
            model[[target]] <<- proposalValue                                
            proposalLogProb <- model$calculate(calcNodes)                     
            logAcceptanceRatio <- proposalLogProb - currentLogProb + g            
            
            accept <- decide(logAcceptanceRatio)                              
            
            if (accept) {
                # no changes to model object needed
                currentLogProb <- proposalLogProb
                currentValue <- proposalValue
                
            } else {
                # reject proposal and revert model to current state
                model[[target]] <<- currentValue
                
                # current full conditional (calculate overwrites the stored value)
                currentLogProb <- model$calculate(calcNodes) 
            }
            
        } # end loop
        
        # synchronize model -> mvSaved after nUpdates
        copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        
    },
    methods = list(   # required method for sampler_BASE base class
        reset = function() {}
    )
)

assign('zUpdate', zUpdate, envir = .GlobalEnv)

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
    lC ~ dinvgamma(cc, dc)
    lH ~ dinvgamma(ch, dh)
    lD ~ dinvgamma(cd, dd)
    
    # indicator variable for model to be used (cases, hosp, deaths)
    z[1:3] ~ dmulti(size = 1, prob = zProb[1:3])
    
    # IDD Curve
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
})

