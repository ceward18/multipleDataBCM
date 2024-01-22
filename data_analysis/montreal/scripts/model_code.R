################################################################################
# Nimble model codes for simulation study
# alarm based on only incidence and deaths
################################################################################


# nimbleFunction for logistic decay IDD curve
logitDecay <- nimbleFunction(     
    run = function(x = double(1), w0 = double(0), k = double(0)) {
        returnType(double(1))
        
        result <- 1 / (1 + exp(k * (x - w0)))
        
        return(result)
    })
assign('logitDecay', logitDecay, envir = .GlobalEnv)

# Hill-Langmuir alarm function
hillAlarm <- nimbleFunction(     
    run = function(x = double(0), nu = double(0), 
                   x0 = double(0), delta = double(0)) {
        returnType(double(0))
        
        result <- delta / (1 + (x0 / x) ^ nu)
        
        return(result)
    })
assign('hillAlarm', hillAlarm, envir = .GlobalEnv)

# nimbleFunction for Gaussian copula with beta distribution
multiBeta <- nimbleFunction(     
    run = function(Z = double(1)) {
        returnType(double(1))
        
        fInv <- pnorm(Z)
        result <- qbeta(fInv, 1, 1)
        
        return(result)
    })
assign('multiBeta', multiBeta, envir = .GlobalEnv)

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
# Model scripts:
#   - SIHRD model with alarm based on incidence and deaths (SIHRD_full)
#   - SIR model with alarm based on incidence and deaths (SIR_full)
#   - SIR model with alarm based on incidence only (SIR_inc)
#   - SIHRD model with no alarm (SIHRD_noAlarm)
#   - SIR model with no alarm (SIR_noAlarm)
# Models either allow for undetected cases to be estimated (undetected), or
#   assume all cases were observed (casesOnly)
################################################################################

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
        
        # compute alarms for each component
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, delta[1])
        alarmD[t] <- hillAlarm(smoothD[t], nuD, x0D, delta[2])
        
        # weighted sum of each component
        alarm[t] <- alarmC[t] + alarmD[t]
        
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
    
    ### compute alarms for summaries
    for (i in 1:n) {
        
        # compute alarms for each component
        yAlarmC[i] <- hillAlarm(xC[i], nuC, x0C, delta[1])
        yAlarmD[i] <- hillAlarm(xD[i], nuD, x0D, delta[2])
        
    }
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(detectA, detectB)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(lambdaShape, lambdaRate) # IH (mean 0.05)
    phi ~ dgamma(phiShape, phiRate)    # HD (mean 0.1)
    
    # alarm functions
    
    # delta prior using Gaussian copula
    Z[1:2] ~ dmnorm(zeros[1:2], cov = Sigma[1:2, 1:2])
    delta[1:2] <- multiBeta(Z[1:2])
    
    nuC ~ dinvgamma(11, 40)
    nuD ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    x0D ~ dunif(minD, maxD)
    
    # constrain deltas to sum to 1
    constrain_deltas ~ dconstraint(delta[1] + delta[2] <= 1)
    
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
        
        # compute alarms for each component
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, delta[1])
        alarmD[t] <- hillAlarm(smoothD[t], nuD, x0D, delta[2])
        
        # weighted sum of each component
        alarm[t] <- alarmC[t] + alarmD[t]
        
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
    probDetect ~ dbeta(detectA, detectB)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(lambdaShape, lambdaRate) # IH (mean 0.1)
    phi ~ dgamma(phiShape, phiRate)    # HD (mean 0.1)
    
    # alarm functions
    
    # delta prior using Gaussian copula
    Z[1:2] ~ dmnorm(zeros[1:2], cov = Sigma[1:2, 1:2])
    delta[1:2] <- multiBeta(Z[1:2])
    
    nuC ~ dinvgamma(11, 40)
    nuD ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    x0D ~ dunif(minD, maxD)
    
    # constrain deltas to sum to 1
    constrain_deltas ~ dconstraint(delta[1] + delta[2] <= 1)
    
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
        
        # compute alarms for each component
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, delta[1])
        alarmD[t] <- hillAlarm(smoothD[t], nuD, x0D, delta[2])
        
        # weighted sum of each component
        alarm[t] <- alarmC[t] + alarmD[t]
        
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
    
    ### compute alarms for summaries
    for (i in 1:n) {
        
        # compute alarms for each component
        yAlarmC[i] <- hillAlarm(xC[i], nuC, x0C, delta[1])
        yAlarmD[i] <- hillAlarm(xD[i], nuD, x0D, delta[2])
        
    }
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(lambdaShape, lambdaRate) # IH (mean 0.05)
    phi ~ dgamma(phiShape, phiRate)    # HD (mean 0.1)
    
    # alarm functions
    
    # delta prior using Gaussian copula
    Z[1:2] ~ dmnorm(zeros[1:2], cov = Sigma[1:2, 1:2])
    delta[1:2] <- multiBeta(Z[1:2])
    
    nuC ~ dinvgamma(11, 40)
    nuD ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    x0D ~ dunif(minD, maxD)
    
    # constrain deltas to sum to 1
    constrain_deltas ~ dconstraint(delta[1] + delta[2] <= 1)
    
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
        
        # compute alarms for each component
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, delta[1])
        alarmD[t] <- hillAlarm(smoothD[t], nuD, x0D, delta[2])
        
        # weighted sum of each component
        alarm[t] <- alarmC[t] + alarmD[t]
        
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
    lambda ~ dgamma(lambdaShape, lambdaRate) # IH (mean 0.1)
    phi ~ dgamma(phiShape, phiRate)    # HD (mean 0.1)
    
    # alarm functions
    
    # delta prior using Gaussian copula
    Z[1:2] ~ dmnorm(zeros[1:2], cov = Sigma[1:2, 1:2])
    delta[1:2] <- multiBeta(Z[1:2])
    
    nuC ~ dinvgamma(11, 40)
    nuD ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    x0D ~ dunif(minD, maxD)
    
    # constrain deltas to sum to 1
    constrain_deltas ~ dconstraint(delta[1] + delta[2] <= 1)
    
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
    probDetect ~ dbeta(detectA, detectB)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # transitions
    gamma1 ~ dgamma(20, 100) # IR (mean 0.2)
    gamma2 ~ dgamma(20, 100) # HR (mean 0.2)
    lambda ~ dgamma(lambdaShape, lambdaRate) # IH (mean 0.05)
    phi ~ dgamma(phiShape, phiRate)    # HD (mean 0.1)
    
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
    lambda ~ dgamma(lambdaShape, lambdaRate) # IH (mean 0.05)
    phi ~ dgamma(phiShape, phiRate)    # HD (mean 0.1)
    
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
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarms for each component
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, delta[1])
        alarmD[t] <- hillAlarm(smoothD[t], nuD, x0D, delta[2])
        
        # weighted sum of each component
        alarm[t] <- alarmC[t] + alarmD[t]
        
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
    
    
    ### compute alarms for summaries
    for (i in 1:n) {
        
        # compute alarms for each component
        yAlarmC[i] <- hillAlarm(xC[i], nuC, x0C, delta[1])
        yAlarmD[i] <- hillAlarm(xD[i], nuD, x0D, delta[2])
        
    }
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(detectA, detectB)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    # delta prior using Gaussian copula
    Z[1:2] ~ dmnorm(zeros[1:2], cov = Sigma[1:2, 1:2])
    delta[1:2] <- multiBeta(Z[1:2])
    
    nuC ~ dinvgamma(11, 40)
    nuD ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    x0D ~ dunif(minD, maxD)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.5)
    k ~ dgamma(100, 100)
    
    # constrain deltas to sum to 1
    constrain_deltas ~ dconstraint(delta[1] + delta[2] <= 1)
    
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
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarms for each component
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, delta[1])
        alarmD[t] <- hillAlarm(smoothD[t], nuD, x0D, delta[2])
        
        # weighted sum of each component
        alarm[t] <- alarmC[t] + alarmD[t]
        
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
    
    
    ### compute alarms for summaries
    for (i in 1:n) {
        
        # compute alarms for each component
        yAlarmC[i] <- hillAlarm(xC[i], nuC, x0C, delta[1])
        yAlarmD[i] <- hillAlarm(xD[i], nuD, x0D, delta[2])
        
    }
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    # delta prior using Gaussian copula
    Z[1:2] ~ dmnorm(zeros[1:2], cov = Sigma[1:2, 1:2])
    delta[1:2] <- multiBeta(Z[1:2])
    
    nuC ~ dinvgamma(11, 40)
    nuD ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    x0D ~ dunif(minD, maxD)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.5)
    k ~ dgamma(100, 100)
    
    # constrain deltas to sum to 1
    constrain_deltas ~ dconstraint(delta[1] + delta[2] <= 1)
    
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
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        probSI[t] <- 1 - exp(- beta * sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
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
    probDetect ~ dbeta(detectA, detectB)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.5)
    k ~ dgamma(100, 100)
    
    
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
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
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
    w0 ~ dnorm(5, sd = 0.5)
    k ~ dgamma(100, 100)
    
    
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
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # alarm is function of incidence only
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, deltaC)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarmC[t]) * 
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
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarmC[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    
    ### compute alarms for summaries
    for (i in 1:n) {
        
        # compute alarms for each component
        yAlarmC[i] <- hillAlarm(xC[i], nuC, x0C, deltaC)
        
    }
    
    ### Priors
    
    # detection probability (1/4 reported)
    probDetect ~ dbeta(detectA, detectB)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    deltaC ~ dunif(0, 1)
    nuC ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.5)
    k ~ dgamma(100, 100)
    
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
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    # time 1
    alarmC[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarmC[1]) * 
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
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, deltaC)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarmC[t]) * 
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
    probDetect ~ dbeta(detectA, detectB)
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    deltaC ~ dunif(0, 1)
    nuC ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.5)
    k ~ dgamma(100, 100)
    
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
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # alarm is function of incidence only
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, deltaC)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarmC[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        # SIR model
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarmC[1:tau]), 
                                 N = N, S = S[1:tau], maxInf = maxInf,
                                 iddCurve = idd_curve[1:maxInf])
    
    
    ### compute alarms for summaries
    for (i in 1:n) {
        
        # compute alarms for each component
        yAlarmC[i] <- hillAlarm(xC[i], nuC, x0C, deltaC)
        
    }
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    
    # alarm functions
    deltaC ~ dunif(0, 1)
    nuC ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.5)
    k ~ dgamma(100, 100)
    
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
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    # time 1
    alarmC[1] <- 0
    
    probSI[1] <- 1 - exp(- beta * (1 - alarmC[1]) * 
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
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, deltaC)
        
        probSI[t] <- 1 - exp(- beta * (1 - alarmC[t]) * 
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
    deltaC ~ dunif(0, 1)
    nuC ~ dinvgamma(11, 40)
    x0C ~ dunif(minC, maxC)
    
    # IDD Curve
    w0 ~ dnorm(5, sd = 0.5)
    k ~ dgamma(100, 100)
    
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

