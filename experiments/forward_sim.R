


dat <- read.csv('./Data/nycClean.csv')

library(nimble)

# Hill-Langmuir alarm function
hillAlarm <- nimbleFunction(     
    run = function(x = double(0), nu = double(0), x0 = double(0), delta = double(0)) {
        returnType(double(0))
        
        result <- delta / (1 + (x0 / x) ^ nu)
        
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


SIR_hill_multi <-  nimbleCode({
    
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
        
        smoothC[t] <- movingAverage(Istar[1:(t-1)], bw)[t-1]
        smoothH[t] <- movingAverage(Hstar[1:(t-1)], bw)[t-1]
        smoothD[t] <- movingAverage(Dstar[1:(t-1)], bw)[t-1]
        
        # compute alarms for each component
        alarmC[t] <- hillAlarm(smoothC[t], nuC, x0C, deltaC)
        alarmH[t] <- hillAlarm(smoothH[t], nuH, x0H, deltaH)
        alarmD[t] <- hillAlarm(smoothD[t], nuD, x0D, deltaD)
        
        # weighted sum of each component
        alarm[t] <- alpha[1] * alarmC[t] + 
            alpha[2] * alarmH[t] + 
            alpha[3] * alarmD[t]
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * I[t] / N)
        
        # SIHRD model
        # S -> I
        Istar[t] ~ dbin(probSI[t], S[t])
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
    
    ### Priors
    
    # transmission
    beta ~ dgamma(0.1, 0.1)
    alpha[1:3] ~ ddirch(priorWeights[1:3])
    
    # transitions
    gamma1 ~ dgamma(0.1, 0.1)
    gamma2 ~ dgamma(0.1, 0.1)
    lambda ~ dgamma(0.1, 0.1)
    phi ~ dgamma(0.1, 0.1)
    
    # alarm function
    deltaC ~ dbeta(1, 1)
    nuC ~ dunif(0, 50)
    x0C ~ dunif(0, 20000)
    deltaH ~ dbeta(1, 1)
    nuH ~ dunif(0, 50)
    x0H ~ dunif(0, 20000)
    deltaD ~ dbeta(1, 1)
    nuD ~ dunif(0, 50)
    x0D ~ dunif(0, 20000)
    
})



# simulate from this model


N <- 1e6
I0 <- 5
S0 <- N - I0
SIR_init <- c(S0, I0, 0)
tau <- 100



################################################################################

constantsList <- list(N = N, 
                      tau = tau,
                      S0 = S0,
                      I0 = I0,
                      bw = 60)


### 14-day Incidence
SIR_sim_model <- nimbleModel(code = SIR_hill_multi,
                            constants = constantsList)

# cModel <- compileNimble(SIR_sim_model)

dataNodes <- c('Istar', 'fromI', 'fromH')
dataNodes <- SIR_sim_model$expandNodeNames(dataNodes)
sim_R <- simulator(SIR_sim_model, dataNodes)

# sim_C <- compileNimble(sim_R, project = SIR_sim_model)

# https://onlinelibrary.wiley.com/doi/full/10.1002/rnc.5728
gamma1 <- 0.388  # I to Recovered
gamma2 <- 0.356  # H to Recovered
lambda <- 0.15 # I to H 
phi <- 0.15     # H to D 

# R0 = beta / (gamma1 + lambda) * S(t) / N
0.85 / (gamma1 + lambda)

trueVals <- c(beta = 0.85, 
              gamma1 = gamma1,
              gamma2 = gamma2,
              lambda = lambda,
              phi = phi,
              deltaC = 1,
              nuC = 6,
              x0C = 900,
              deltaH = 1,
              nuH = 6, 
              x0H = 200,
              deltaD = 1,
              nuD = 6,
              x0D = 200,
              alpha = c(0.2, 0.5, 0.3))


tmp <- sim_R$run(trueVals, 1)
names(tmp) <- SIR_sim_model$expandNodeNames(dataNodes,
                                            returnScalarComponents = TRUE)

par(mfrow = c(2,3))
plot(tmp[grep('Istar', names(tmp))], main = 'incidence', type = 'l', ylim = c(0, 8000))
lines(dat$dailyCases, col = 'red')
plot(tmp[grep('fromI.*1\\]', names(tmp))], main = 'hospitalizations', type = 'l', ylim = c(0, 2000))
lines(dat$dailyHosp, col = 'red')
plot(tmp[grep('fromH.*2\\]', names(tmp))], main = 'deaths', type = 'l', ylim = c(0, 800))
lines(dat$dailyDeaths, col = 'red')
plot(hillAlarm(1:2500, trueVals['nuC'], trueVals['x0C'], trueVals['deltaC']), 
     type = 'l', ylab = 'alarm', ylim = c(0,1))
plot(hillAlarm(1:600, trueVals['nuH'], trueVals['x0H'], trueVals['deltaH']), 
     type = 'l', ylab = 'alarm', ylim = c(0,1))
plot(hillAlarm(1:250, trueVals['nuD'], trueVals['x0D'], trueVals['deltaD']), 
     type = 'l', ylab = 'alarm', ylim = c(0,1))


# what does alarm function look like over time?

par(mfrow = c(1,3))
plot(SIR_sim_model$alarmC, type = 'l')
plot(SIR_sim_model$alarmH, type = 'l')
plot(SIR_sim_model$alarmD, type = 'l')

par(mfrow = c(1,1))
plot(SIR_sim_model$alarm)



