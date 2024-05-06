################################################################################
# Nimble model codes for simulation study
# alarm based on only incidence and deaths
# add in only a proportion of cases are detected
################################################################################


# nimbleFunction for logistic decay IDD curve
logitDecay <- nimbleFunction(     
    run = function(x = double(1), w0 = double(0), nu = double(0)) {
        returnType(double(1))
        
        result <- 1 / (1 + exp(nu * (x - w0)))
        
        return(result)
    })
assign('logitDecay', logitDecay, envir = .GlobalEnv)

# power alarm function
powerAlarm <- nimbleFunction(     
    run = function(x = double(0), N = double(0), k = double(0)) {
        returnType(double(0))
        
        result <- 1 - (1 - x / N)^(1 / k)
        
        return(result)
    })
assign('powerAlarm', powerAlarm, envir = .GlobalEnv)

# power alarm function with multiple data sources
powerAlarm2 <- nimbleFunction(     
    run = function(x1 = double(0), x2 = double(0), N = double(0), k = double(0)) {
        returnType(double(0))
        
        result <- 1 - (1 - (x1 + x2) / N)^(1 / k)
        
        return(result)
    })
assign('powerAlarm2', powerAlarm2, envir = .GlobalEnv)


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

# get effective reproductive number for full model with transitions to H and R
get_R0_full <- nimbleFunction(     
    run = function(betat = double(1), N = double(0), gamma1 = double(0), 
                   lambda = double(0), S = double(1), maxInf = double(0)) {
        returnType(double(1))
        
        # probability of transition given 1 infectious individual (vector length t)
        pi_SI <- 1 - exp(- betat / N )
        
        # infectious probability of removal
        pi_IR <- 1 - exp(-gamma1)
        pi_IH <- 1 - exp(-lambda)
        
        # for exponential periods
        multVec <- c(1, rep(NA, maxInf))
        for (i in 2:length(multVec)) {
            multVec[i] <- multVec[i-1] * (1 - pi_IR - pi_IH)
        }
        
        nTime <- length(pi_SI)
        bw <- length(multVec)
        sumSmooth <- rep(NA, nTime - bw)
        for(k in 1:(nTime - bw)){
            t1 <- k
            t2 <- k + bw - 1
            sumSmooth[k] <- sum(pi_SI[t1:t2] * multVec) * S[k]
        }
        
        return(sumSmooth)
    })
assign('get_R0_full', get_R0_full, envir = .GlobalEnv)

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


################################################################################
# SIHRD Model

# model used for simulation purposes

SIHRD_sim <-  nimbleCode({
    
    S[1] <- S0
    I[1] <- I0
    H[1] <- 0
    R[1] <- 0
    D[1] <- 0
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    fromIProb[1:3] <- c(probIH, probIR, 1 - probIH - probIR)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    fromHProb[1:3] <- c(probHR, probHD, 1 - probHR - probHD)
    
    ### first time point
    
    # weighted sum of each component
    alarm[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    # SIHRD model
    Istar[1] ~ dbin(probSI[1], S[1])
    # Istar[t] <- detectIstar[t] + undetectIstar[t]
    detectIstar[1] ~ dbin(probDetect, Istar[1]) # e.g. 25% of cases are detected
    
    fromI[1, 1:3] ~ dmulti(prob = fromIProb[1:3], size = I[1])
    Hstar[1] <- fromI[1, 1]
    RstarI[1] <- fromI[1, 2]
    fromH[1, 1:3] ~ dmulti(prob = fromHProb[1:3], size = H[1])
    RstarH[1] <- fromH[1, 1]
    Dstar[1] <- fromH[1, 2]
    
    # update S, I, H, R, D
    S[2] <- S[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Hstar[1] - RstarI[1]
    H[2] <- H[1] + Hstar[1] - RstarH[1] - Dstar[1] 
    R[2] <- R[1] + RstarI[1] + RstarH[1]
    D[2] <- D[1] + Dstar[1] 
    
    ### rest of time points
    for(t in 2:tau) {
        
        # alarm is only based on observed data
        smoothC[t] <- movingAverage(detectIstar[1:(t-1)], bw)[t-1]
        smoothD[t] <- movingAverage(Dstar[1:(t-1)], bw)[t-1]
        
        # compute alarms for each component
        alarm[t] <- powerAlarm2(alpha * smoothC[t], (1 - alpha) * smoothD[t],
                                N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        # I -> H or R
        fromI[t, 1:3] ~ dmulti(prob = fromIProb[1:3], size = I[t])
        Hstar[t] <- fromI[t, 1]
        RstarI[t] <- fromI[t, 2]
        # H -> R or D
        fromH[t, 1:3] ~ dmulti(prob = fromHProb[1:3], size = H[t])
        RstarH[t] <- fromH[t, 1]
        Dstar[t] <- fromH[t, 2]
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    # estimated effective R0
    R0[1:(tau-11)] <- get_R0_full(betat = beta * (1 - alarm[1:tau]), 
                                  N = N, gamma1 = gamma1, lambda = lambda,
                                  S = S[1:tau], maxInf = 10)
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions - priors don't matter as this is used for simulation only
    k ~ dgamma(0.1, 0.1)
    
    # relative importance of cases
    alpha ~ dbeta(1,1)
    
})

################################################################################
# Model scripts:
#   - SIHRD model with alarm based on incidence and deaths (SIHRD_full)
#   - SIR model with alarm based on incidence and deaths (SIR_full)
#   - SIR model with alarm based on incidence only (SIR_inc)
#   - SIHRD model with no alarm (SIHRD_noAlarm)
#   - SIR model with no alarm (SIR_noAlarm)
# Models either allow for undetected cases to be estimated (undetected), or
#   assume all cases were observed (casesOnly)
# Second model codes for simulation are needed for SIHRD_full, SIR_inc because
#   smoothed cases/deaths must be computed for simulation
################################################################################

# SIHRD model with alarm based on incidence and deaths
SIHRD_full_undetected <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### loop over time
    for(t in 1:tau) {
        
        # weighted sum of each component
        alarm[t] <- powerAlarm2(x1 = alpha * smoothC[t], 
                                x2 = (1 - alpha) * smoothD[t],
                                N = N, k = k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf-1)] <- get_R0_full(betat = beta * (1 - alarm[1:tau]), 
                                        N = N, gamma1 = gamma1, lambda = lambda,
                                        S = S[1:tau], maxInf = maxInf)
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    alpha ~ dbeta(1, 1)
    
})

################################################################################

# SIHRD model with alarm based on incidence and deaths
SIHRD_full_undetected_sim <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### first time point
    alarm[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    # SIHRD model
    # S -> I
    Istar[1] ~ dbin(probSI[1], S[1])
    # Istar[t] <- detectIstar[t] + undetectIstar[t]
    detectIstar[1] ~ dbin(probDetect, Istar[1]) # e.g. 25% of cases are detected
    # I -> H or R using sequential binomial
    Hstar[1] ~ dbin(probIH, I[1])
    RstarI[1] ~ dbin(probIR / (1 - probIH), I[1] - Hstar[1])
    # H -> R or D using sequential binomial
    Dstar[1] ~ dbin(probHD, H[1])
    RstarH[1] ~ dbin(probHR / (1 - probHD), H[1] - Dstar[1])
    
    # update S, I, H, R, D
    S[2] <- S[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Hstar[1] - RstarI[1]
    H[2] <- H[1] + Hstar[1] - RstarH[1] - Dstar[1] 
    R[2] <- R[1] + RstarI[1] + RstarH[1]
    D[2] <- D[1] + Dstar[1] 
    
    
    ### loop over rest of time points
    for(t in 2:tau) {
        
        # compute moving average up to t-1
        smoothC[t] <- movingAverage(detectIstar[1:(t-1)], bw)[t-1]
        smoothD[t] <- movingAverage(Dstar[1:(t-1)], bw)[t-1]
        
        # weighted sum of each component
        alarm[t] <- powerAlarm2(alpha * smoothC[t], (1 - alpha) * smoothD[t],
                                N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    alpha ~ dbeta(1, 1)
    
})

################################################################################

# SIHRD model with alarm based on incidence and deaths
SIHRD_full_casesOnly <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### loop over time
    for(t in 1:tau) {
        
        # weighted sum of each component
        alarm[t] <- powerAlarm2(alpha * smoothC[t], (1 - alpha) * smoothD[t],
                                N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf-1)] <- get_R0_full(betat = beta * (1 - alarm[1:tau]), 
                                        N = N, gamma1 = gamma1, lambda = lambda,
                                        S = S[1:tau], maxInf = maxInf)
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    alpha ~ dbeta(1, 1)
    
})

################################################################################

# SIHRD model with alarm based on incidence and deaths
SIHRD_full_casesOnly_sim <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### first time point
    alarm[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    # SIHRD model
    # S -> I
    Istar[1] ~ dbin(probSI[1], S[1])
    # I -> H or R using sequential binomial
    Hstar[1] ~ dbin(probIH, I[1])
    RstarI[1] ~ dbin(probIR / (1 - probIH), I[1] - Hstar[1])
    # H -> R or D using sequential binomial
    Dstar[1] ~ dbin(probHD, H[1])
    RstarH[1] ~ dbin(probHR / (1 - probHD), H[1] - Dstar[1])
    
    # update S, I, H, R, D
    S[2] <- S[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Hstar[1] - RstarI[1]
    H[2] <- H[1] + Hstar[1] - RstarH[1] - Dstar[1] 
    R[2] <- R[1] + RstarI[1] + RstarH[1]
    D[2] <- D[1] + Dstar[1] 
    
    
    ### loop over rest of time points
    for(t in 2:tau) {
        
        # compute moving average up to t-1
        smoothC[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
        smoothD[t] <- movingAverage(Dstar[1:(t-1)], bw)[t-1]
        
        # weighted sum of each component
        alarm[t] <- powerAlarm2(alpha * smoothC[t], (1 - alpha) * smoothD[t],
                                N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    alpha ~ dbeta(1, 1)
    
})

################################################################################

# SIHRD model with alarm based on incidence only
SIHRD_inc_undetected <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### loop over time
    for(t in 1:tau) {
        
        # alarm based on cases only
        alarm[t] <- powerAlarm(smoothC[t], N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf-1)] <- get_R0_full(betat = beta * (1 - alarm[1:tau]), 
                                        N = N, gamma1 = gamma1, lambda = lambda,
                                        S = S[1:tau], maxInf = maxInf)
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
})

################################################################################

# SIHRD model with alarm based on incidence only
SIHRD_inc_undetected_sim <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### first time point
    alarm[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    # SIHRD model
    # S -> I
    Istar[1] ~ dbin(probSI[1], S[1])
    # Istar[t] <- detectIstar[t] + undetectIstar[t]
    detectIstar[1] ~ dbin(probDetect, Istar[1]) # e.g. 25% of cases are detected
    # I -> H or R using sequential binomial
    Hstar[1] ~ dbin(probIH, I[1])
    RstarI[1] ~ dbin(probIR / (1 - probIH), I[1] - Hstar[1])
    # H -> R or D using sequential binomial
    Dstar[1] ~ dbin(probHD, H[1])
    RstarH[1] ~ dbin(probHR / (1 - probHD), H[1] - Dstar[1])
    
    # update S, I, H, R, D
    S[2] <- S[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Hstar[1] - RstarI[1]
    H[2] <- H[1] + Hstar[1] - RstarH[1] - Dstar[1] 
    R[2] <- R[1] + RstarI[1] + RstarH[1]
    D[2] <- D[1] + Dstar[1] 
    
    
    ### loop over rest of time points
    for(t in 2:tau) {
        
        # compute moving average up to t-1
        smoothC[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
        
        # alarm based on cases only
        alarm[t] <- powerAlarm(smoothC[t], N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf-1)] <- get_R0_full(betat = beta * (1 - alarm[1:tau]), 
                                        N = N, gamma1 = gamma1, lambda = lambda,
                                        S = S[1:tau], maxInf = maxInf)
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
})


################################################################################

# SIHRD model with alarm based on incidence only
SIHRD_inc_casesOnly <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### loop over time
    for(t in 1:tau) {
        
        # alarm based on cases only
        alarm[t] <- powerAlarm(smoothC[t], N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf-1)] <- get_R0_full(betat = beta * (1 - alarm[1:tau]), 
                                        N = N, gamma1 = gamma1, lambda = lambda,
                                        S = S[1:tau], maxInf = maxInf)
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
})

################################################################################

# SIHRD model with alarm based on incidence only
SIHRD_inc_casesOnly_sim <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### first time point
    alarm[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * I[1] / N)
    
    # SIHRD model
    # S -> I
    Istar[1] ~ dbin(probSI[1], S[1])
    # I -> H or R using sequential binomial
    Hstar[1] ~ dbin(probIH, I[1])
    RstarI[1] ~ dbin(probIR / (1 - probIH), I[1] - Hstar[1])
    # H -> R or D using sequential binomial
    Dstar[1] ~ dbin(probHD, H[1])
    RstarH[1] ~ dbin(probHR / (1 - probHD), H[1] - Dstar[1])
    
    # update S, I, H, R, D
    S[2] <- S[1] - Istar[1]
    I[2] <- I[1] + Istar[1] - Hstar[1] - RstarI[1]
    H[2] <- H[1] + Hstar[1] - RstarH[1] - Dstar[1] 
    R[2] <- R[1] + RstarI[1] + RstarH[1]
    D[2] <- D[1] + Dstar[1] 
    
    
    ### loop over rest of time points
    for(t in 2:tau) {
        
        # compute moving average up to t-1
        smoothC[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
        
        # alarm based on cases only
        alarm[t] <- powerAlarm(smoothC[t], N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf-1)] <- get_R0_full(betat = beta * (1 - alarm[1:tau]), 
                                        N = N, gamma1 = gamma1, lambda = lambda,
                                        S = S[1:tau], maxInf = maxInf)
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
})

################################################################################

# SIHRD model with no alarm
SIHRD_noAlarm_undetected <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### loop over time
    for(t in 1:tau) {
        
        probSI[t] <- 1 - exp(- beta * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t] 
        
        betat[t] <- beta
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf-1)] <- get_R0_full(betat = betat[1:tau], 
                                        N = N, gamma1 = gamma1, lambda = lambda,
                                        S = S[1:tau], maxInf = maxInf)
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
})

################################################################################

# SIHRD model with no alarm
SIHRD_noAlarm_casesOnly <-  nimbleCode({
    
    # compartment initial values
    comp_init[1:5] ~ dmulti(prob = initProb[1:5], size = N)
    S[1] <- comp_init[1] - 1
    I[1] <- comp_init[2] + 1 # add 1 to ensure I0 > 0
    H[1] <- comp_init[3]
    R[1] <- comp_init[4]
    D[1] <- comp_init[5]
    
    probIH <- 1 - exp(-lambda)
    probIR <- 1 - exp(-gamma1)
    
    probHR <- 1 - exp(-gamma2)
    probHD <- 1 - exp(-phi)
    
    ### loop over time
    for(t in 1:tau) {
        
        probSI[t] <- 1 - exp(- beta * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
        # I -> H or R using sequential binomial
        Hstar[t] ~ dbin(probIH, I[t])
        RstarI[t] ~ dbin(probIR / (1 - probIH), I[t] - Hstar[t])
        # H -> R or D using sequential binomial
        Dstar[t] ~ dbin(probHD, H[t])
        RstarH[t] ~ dbin(probHR / (1 - probHD), H[t] - Dstar[t])
        
        # update S, I, H, R, D
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1] <- I[t] + Istar[t] - Hstar[t] - RstarI[t]
        H[t + 1] <- H[t] + Hstar[t] - RstarH[t] - Dstar[t] 
        R[t + 1] <- R[t] + RstarI[t] + RstarH[t]
        D[t + 1] <- D[t] + Dstar[t]
        
        betat[t] <- beta
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf-1)] <- get_R0_full(betat = betat[1:tau], 
                                        N = N, gamma1 = gamma1, lambda = lambda,
                                        S = S[1:tau], maxInf = maxInf)
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(3, 100) # IH (mean 0.03)
    phi ~ dgamma(10, 100)    # HD (mean 0.1)
    
})

################################################################################

# SIR model with alarm based on incidence and deaths
SIR_full_undetected <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 11 = I0, 12 = R0
    comp_init[1:12] ~ dmulti(prob = initProb[1:12], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # weighted sum of each component
        alarm[t] <- powerAlarm2(alpha * smoothC[t], (1 - alpha) * smoothD[t],
                                N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        # SIR model
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    alpha ~ dbeta(1, 1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.25)
    nu ~ dgamma(100, 100)
    
})

################################################################################

# SIR model with alarm based on incidence and deaths
SIR_full_casesOnly <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 11 = I0, 12 = R0
    comp_init[1:12] ~ dmulti(prob = initProb[1:12], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # weighted sum of each component
        alarm[t] <- powerAlarm2(alpha * smoothC[t], (1 - alpha) * smoothD[t],
                                N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        # SIR model
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
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    alpha ~ dbeta(1, 1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.25)
    nu ~ dgamma(100, 100)
    
})

################################################################################

# SIR model with no alarm
SIR_noAlarm_undetected <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 11 = I0, 12 = R0
    comp_init[1:12] ~ dmulti(prob = initProb[1:12], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### rest of time points
    for(t in 1:tau) {
        
        probSI[t] <- 1 - exp(- beta * sum(idd_curve[1:maxInf] * 
                                              I[t, 1:maxInf]) / N)
        
        # SIR model
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
        betat[t] <- beta
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = betat[1:tau], 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.25)
    nu ~ dgamma(100, 100)
    
    
})

################################################################################

# SIR model with no alarm
SIR_noAlarm_casesOnly <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 11 = I0, 12 = R0
    comp_init[1:12] ~ dmulti(prob = initProb[1:12], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### rest of time points
    for(t in 1:tau) {
        
        probSI[t] <- 1 - exp(- beta * sum(idd_curve[1:maxInf] *
                                              I[t, 1:maxInf]) / N)
        
        # SIR model
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
        betat[t] <- beta
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = betat[1:tau], 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.25)
    nu ~ dgamma(100, 100)
    
    
})

################################################################################

# SIR model with incidence only alarm
SIR_inc_undetected <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 11 = I0, 12 = R0
    comp_init[1:12] ~ dmulti(prob = initProb[1:12], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # alarm is function of incidence only
        alarm[t] <- powerAlarm(smoothC[t], N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        # SIR model
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.25)
    nu ~ dgamma(100, 100)
    
})

################################################################################

# SIR model with incidence only alarm
SIR_inc_undetected_sim <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 11 = I0, 12 = R0
    comp_init[1:12] ~ dmulti(prob = initProb[1:12], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    # time 1
    alarm[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * 
                             sum(idd_curve[1:maxInf] * I[1, 1:maxInf]) / N)
    
    # SIR model
    Istar[1] ~ dbin(probSI[1], S[1])
    # Istar[t] <- detectIstar[t] + undetectIstar[t]
    detectIstar[1] ~ dbin(probDetect, Istar[1]) # e.g. 25% of cases are detected
    
    # update S and I
    S[2] <- S[1] - Istar[1]
    I[2, 2:maxInf] <- I[1, 1:(maxInf - 1)]  # shift current I by one day
    I[2, 1] <- Istar[1]                     # add newly infectious
    
    ### rest of time points
    for(t in 2:tau) {
        
        # moving average incidence
        smoothC[t] <- movingAverage(detectIstar[1:(t-1)], bw)[t-1]
        
        # alarm is function of incidence only
        alarm[t] <- powerAlarm(smoothC[t], N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        # SIR model
        Istar[t] ~ dbin(probSI[t], S[t])
        # Istar[t] <- detectIstar[t] + undetectIstar[t]
        detectIstar[t] ~ dbin(probDetect, Istar[t]) # e.g. 25% of cases are detected
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(250, 750)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.25)
    nu ~ dgamma(100, 100)
    
})

################################################################################

# SIR model with incidence only alarm
SIR_inc_casesOnly <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 11 = I0, 12 = R0
    comp_init[1:12] ~ dmulti(prob = initProb[1:12], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # alarm is function of incidence only
        alarm[t] <- powerAlarm(smoothC[t], N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        # SIR model
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
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.25)
    nu ~ dgamma(100, 100)
    
})

################################################################################

# SIR model with incidence only alarm
SIR_inc_casesOnly_sim <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 11 = I0, 12 = R0
    comp_init[1:12] ~ dmulti(prob = initProb[1:12], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    # time 1
    alarm[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarm[1]) * 
                             sum(idd_curve[1:maxInf] * I[1, 1:maxInf]) / N)
    
    # SIR model
    Istar[1] ~ dbin(probSI[1], S[1])
    
    # update S and I
    S[2] <- S[1] - Istar[1]
    I[2, 2:maxInf] <- I[1, 1:(maxInf - 1)]  # shift current I by one day
    I[2, 1] <- Istar[1]                     # add newly infectious
    
    ### rest of time points
    for(t in 2:tau) {
        
        # moving average incidence
        smoothC[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
        
        # alarm is function of incidence only
        alarm[t] <- powerAlarm(smoothC[t], N, k)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        # SIR model
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.25)
    nu ~ dgamma(100, 100)
    
})

################################################################################
### Special proposal function for removal times in exponential model

RstarUpdate <- nimbleFunction(
    name = 'Rstar',                              
    contains = sampler_BASE,                     
    setup = function(model, mvSaved, target, control) {                 # REQUIRED setup arguments
        calcNodes <- model$getDependencies(target) 
        # percent <- if(!is.null(control$percent)) control$percent else 0.05   
        
        # number of update attempts at each iteration
        nUpdates <- if(!is.null(control$nUpdates)) control$nUpdates else 500
        
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
    methods = list(                              # required method for sampler_BASE base class
        reset = function() {}
    )
)

assign('RstarUpdate', RstarUpdate, envir = .GlobalEnv)
