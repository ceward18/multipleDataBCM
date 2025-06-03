################################################################################
# Nimble model codes for simulation study
# alarm based on only incidence and deaths
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


# get smoothC or smoothD at time t when doing posterior prediction
get_smooth <- nimbleFunction(     
    run = function(Xstar = double(1), t = double(0), Xstar0 = double(1), 
                   Xstar0Length = double(0),  bw = double(0)) {
        returnType(double(0))
        
        if (t < bw) {
            
            # incorporate previously observed incidence
            result <- movingAverage(c(Xstar0, Xstar), bw)[Xstar0Length + t - 1]
            
        } else {
            result <- movingAverage(Xstar, bw)[t - 1]
        }
        
        return(result)
    })
assign('get_smooth', get_smooth, envir = .GlobalEnv)

# get effective reproductive number at time t using a forward sum
get_R0 <- nimbleFunction(     
    run = function(betat = double(1), N = double(0), S = double(1),
                   maxInf = double(0), iddCurve = double(1)) {
        returnType(double(1))
        
        nTime <- length(betat)
        bw <- maxInf
        sumSmooth <- rep(NA, nTime - bw)
        for(k in 1:(nTime - bw)){
            t1 <- k
            t2 <- k + bw - 1
            pi_SI <- 1 - exp(-betat[t1:t2] * iddCurve / N)
            sumSmooth[k] <- sum(pi_SI * S[t1:t2])
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
            sumSmooth[k] <- sum(pi_SI[t1:t2] * multVec * S[t1:t2]) 
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
# Model scripts:
#   - SIHRD model with alarm based on incidence and deaths (SIHRD_full)
#   - SIHRD model with alarm based on incidence only (SIHRD_inc)
#   - SIR model with alarm based on incidence and deaths (SIR_full)
#   - SIR model with alarm based on incidence only (SIR_inc)
#   - SIHRD model with no alarm (SIHRD_noAlarm)
#   - SIR model with no alarm (SIR_noAlarm)
################################################################################

################################################################################

# SIHRD model with alarm based on incidence and deaths
SIHRD_full <-  nimbleCode({
    
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
        alarm[t] <- powerAlarm(x = alpha * smoothC[t] + (1 - alpha) * smoothD[t],
                               N = N, 
                               k = k)
        
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
    gamma1 ~ dgamma(1429, 10000)    # IR 
    gamma2 ~ dgamma(6700, 100000)   # HR 
    lambda ~ dgamma(20, 1000)   # IH 
    phi ~ dgamma(50, 1000)      # HD
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    alpha ~ dbeta(1, 1)
    
})

################################################################################

# SIHRD model with alarm based on incidence and deaths
SIHRD_full_sim <-  nimbleCode({
    
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
    smoothC[1] <- smoothC0
    smoothD[1] <- smoothD0
    alarm[1] <- powerAlarm(x = alpha * smoothC[1] + (1 - alpha) * smoothD[1],
                           N = N, 
                           k = k)
    
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
        smoothC[t] <- get_smooth(Istar[1:(t-1)], t, Istar0[1:Istar0Length], Istar0Length, bw)
        smoothD[t] <- get_smooth(Dstar[1:(t-1)], t, Dstar0[1:Dstar0Length], Dstar0Length, bw)
        
        # weighted sum of each component
        alarm[t] <- powerAlarm(x = alpha * smoothC[t] + (1 - alpha) * smoothD[t],
                               N = N, 
                               k = k)
        
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
    gamma1 ~ dgamma(1429, 10000)    # IR 
    gamma2 ~ dgamma(6700, 100000)   # HR 
    lambda ~ dgamma(20, 1000)   # IH 
    phi ~ dgamma(50, 1000)      # HD
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    alpha ~ dbeta(1, 1)
    
})

################################################################################
# SIHRD model with alarm based on incidence only

SIHRD_inc <-  nimbleCode({
    
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
    gamma1 ~ dgamma(1429, 10000)    # IR 
    gamma2 ~ dgamma(6700, 100000)   # HR 
    lambda ~ dgamma(20, 1000)   # IH 
    phi ~ dgamma(50, 1000)      # HD
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
})

################################################################################
# SIHRD model with alarm based on incidence only

SIHRD_inc_sim <-  nimbleCode({
    
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
    smoothC[1] <- smoothC0
    alarm[1] <- powerAlarm(smoothC[1], N, k)
    
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
        smoothC[t] <- get_smooth(Istar[1:(t-1)], t, Istar0[1:Istar0Length], Istar0Length, bw)
        
        # weighted sum of each component
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
    gamma1 ~ dgamma(1429, 10000)    # IR 
    gamma2 ~ dgamma(6700, 100000)   # HR 
    lambda ~ dgamma(20, 1000)   # IH 
    phi ~ dgamma(50, 1000)      # HD
    
    # alarm functions
    k ~ dgamma(0.1, 0.1)
    
})

################################################################################

# SIR model with alarm based on incidence and deaths
SIR_full <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 15 = I0, 16 = R0
    comp_init[1:(maxInf + 2)] ~ dmulti(prob = initProb[1:(maxInf + 2)], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # weighted sum of each component
        alarm[t] <- powerAlarm(x = alpha * smoothC[t] + (1 - alpha) * smoothD[t],
                               N = N, 
                               k = k)
        
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
    w0 ~ dnorm(7, sd = 0.1)
    nu ~ dgamma(1000, 1000)
    
    
    
})

################################################################################

# SIR model with incidence only alarm
SIR_inc <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 15 = I0, 16 = R0
    comp_init[1:(maxInf + 2)] ~ dmulti(prob = initProb[1:(maxInf + 2)], size = N)
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
    w0 ~ dnorm(7, sd = 0.1)
    nu ~ dgamma(1000, 1000)
    
})

################################################################################

# SIR model with incidence only alarm
SIR_inc_sim <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 15 = I0, 16 = R0
    comp_init[1:(maxInf + 2)] ~ dmulti(prob = initProb[1:(maxInf + 2)], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### first time point
    smoothC[1] <- smoothC0
    alarm[1] <- powerAlarm(smoothC[1], N, k)
    
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
        
        # compute moving average up to t-1
        smoothC[t] <- get_smooth(Istar[1:(t-1)], t, 
                                 Istar0[1:Istar0Length], Istar0Length, bw)
        
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
    w0 ~ dnorm(7, sd = 0.1)
    nu ~ dgamma(1000, 1000)
    
})

################################################################################

# SIR model with no alarm
SIR_noAlarm <-  nimbleCode({
    
    # compartment initial values
    # 1 = S0, 2 - 15 = I0, 16 = R0
    comp_init[1:(maxInf + 2)] ~ dmulti(prob = initProb[1:(maxInf + 2)], size = N)
    S[1] <- comp_init[1] - 1
    I[1, 1] <- comp_init[2] + 1
    I[1, 2:maxInf] <- comp_init[3:(maxInf + 1)]
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, nu)
    
    ### rest of time points
    for(t in 1:tau) {
        
        probSI[t] <- 1 - exp(- beta * sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
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
    w0 ~ dnorm(7, sd = 0.1)
    nu ~ dgamma(1000, 1000)
    
    
})

################################################################################

# SIHRD model with no alarm
SIHRD_noAlarm <-  nimbleCode({
    
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
    gamma1 ~ dgamma(1429, 10000)    # IR 
    gamma2 ~ dgamma(6700, 100000)   # HR
    lambda ~ dgamma(20, 1000)   # IH 
    phi ~ dgamma(50, 1000)      # HD
    
})

################################################################################
### Special proposal function for removal times in exponential model

RstarUpdate <- nimbleFunction(
    name = 'Rstar',                              
    contains = sampler_BASE,                     
    setup = function(model, mvSaved, target, control) {                 # REQUIRED setup arguments
        calcNodes <- model$getDependencies(target) 
        # percent <- if(!is.null(control$percent)) control$percent else 0.05   
        
        # number of update attempts
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
    methods = list(                              # required method for sampler_BASE base class
        reset = function() {}
    )
)

assign('RstarUpdate', RstarUpdate, envir = .GlobalEnv)


