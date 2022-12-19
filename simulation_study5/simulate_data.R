################################################################################
# Simulate data under various scenarios:
# 
### deltaC = 0.6, deltaH = 0.1, deltaD = 0.1
### deltaC = 0.1, deltaH = 0.6, deltaD = 0.1
### deltaC = 0.1, deltaH = 0.1, deltaD = 0.6
################################################################################

### load libraries
library(nimble)

### source model code
source('./scripts/model_code.R')

### Number of simulated epidemics
nSim <- 50


### Constants for all models
N <- 1e6
I0 <- 5
S0 <- N - I0
tau <- 50
bw <- 30

### parameter values for all models 
# (from https://onlinelibrary.wiley.com/doi/full/10.1002/rnc.5728)
gamma1 <- 0.488  # I to Recovered
gamma2 <- 0.356  # H to Recovered
lambda <- 0.0279 # I to H 
phi <- 0.122     # H to D 

gamma1 <- 0.2  # I to Recovered
gamma2 <- 0.2  # H to Recovered
lambda <- 0.1 # I to H 
phi <- 0.1     # H to D 

################################################################################
### Set up simulator

constantsList <- list(N = N, 
                      tau = tau,
                      S0 = S0,
                      I0 = I0,
                      bw = bw)

SIR_sim_model <- nimbleModel(code = SIHRD_sim,
                             constants = constantsList)
SIR_sim_model_C <- compileNimble(SIR_sim_model)

dataNodes <- c('Istar', 'fromI', 'fromH')
dataNodes <- SIR_sim_model$expandNodeNames(dataNodes)
sim_R <- simulator(SIR_sim_model, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_sim_model)

################################################################################
### High incidence alarm

# delta values must sum to less than 1
deltaC <- 0.6
deltaH <- 0.1
deltaD <- 0.1

trueVals <- c(beta = 0.85, 
              gamma1 = gamma1,
              gamma2 = gamma2,
              lambda = lambda,
              phi = phi,
              deltaC = deltaC,
              nuC = 3,
              x0C = 100,
              deltaH = deltaH,
              nuH = 3,
              x0H = 25,
              deltaD = deltaD,
              nuD = 3,
              x0D = 10)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- SIR_sim_model$expandNodeNames(dataNodes,
                                                   returnScalarComponents = TRUE)

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

plot(toSave[1,grep('Istar', colnames(toSave))], type = 'l', col = 'grey', 
     ylim = c(0, 1000))
for (i in 2:nrow(toSave)) {
    lines(toSave[i,grep('Istar', colnames(toSave))], col = 'grey')
    lines(toSave[i,grep('fromI.*1\\]', colnames(toSave))], col = 'tomato')
    lines(toSave[i,grep('fromH.*2\\]', colnames(toSave))], col = 'lightblue')
}

peaks <- matrix(NA, ncol = 3, nrow = nSim)
for (i in 1:nrow(toSave)) {
    peaks[i,] <- c(which.max(toSave[i,grep('Istar', colnames(toSave))]),
                   which.max(toSave[i,grep('fromI.*1\\]', colnames(toSave))]),
                   which.max(toSave[i,grep('fromH.*2\\]', colnames(toSave))]))
}
colMeans(peaks)

saveRDS(toSave, './data/sim_inc.rds')

################################################################################
### High hospitalization alarm

# delta values must sum to less than 1
deltaC <- 0.1
deltaH <- 0.6
deltaD <- 0.1

trueVals <- c(beta = 0.85, 
              gamma1 = gamma1,
              gamma2 = gamma2,
              lambda = lambda,
              phi = phi,
              deltaC = deltaC,
              nuC = 3,
              x0C = 100,
              deltaH = deltaH,
              nuH = 3,
              x0H = 25,
              deltaD = deltaD,
              nuD = 3,
              x0D = 10)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- SIR_sim_model$expandNodeNames(dataNodes,                                                    returnScalarComponents = TRUE)

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]


plot(toSave[1,grep('Istar', colnames(toSave))], type = 'l', col = 'grey', ylim = c(0, 1000))
for (i in 2:nrow(toSave)) {
    lines(toSave[i,grep('Istar', colnames(toSave))], col = 'grey')
    lines(toSave[i,grep('fromI.*1\\]', colnames(toSave))], col = 'tomato')
    lines(toSave[i,grep('fromH.*2\\]', colnames(toSave))], col = 'lightblue')
}

peaks <- matrix(NA, ncol = 3, nrow = nSim)
for (i in 1:nrow(toSave)) {
    peaks[i,] <- c(which.max(toSave[i,grep('Istar', colnames(toSave))]),
                   which.max(toSave[i,grep('fromI.*1\\]', colnames(toSave))]),
                   which.max(toSave[i,grep('fromH.*2\\]', colnames(toSave))]))
}
colMeans(peaks)

saveRDS(toSave, './data/sim_hosp.rds')


################################################################################
### High deaths alarm

# delta values must sum to less than 1
deltaC <- 0.1
deltaH <- 0.1
deltaD <- 0.6

trueVals <- c(beta = 0.85, 
              gamma1 = gamma1,
              gamma2 = gamma2,
              lambda = lambda,
              phi = phi,
              deltaC = deltaC,
              nuC = 3,
              x0C = 100,
              deltaH = deltaH,
              nuH = 3,
              x0H = 25,
              deltaD = deltaD,
              nuD = 3,
              x0D = 10)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- SIR_sim_model$expandNodeNames(dataNodes,                                                    returnScalarComponents = TRUE)

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

plot(toSave[1,grep('Istar', colnames(toSave))], type = 'l', col = 'grey', ylim = c(0, 1000))
for (i in 2:nrow(toSave)) {
    lines(toSave[i,grep('Istar', colnames(toSave))], col = 'grey')
    lines(toSave[i,grep('fromI.*1\\]', colnames(toSave))], col = 'tomato')
    lines(toSave[i,grep('fromH.*2\\]', colnames(toSave))], col = 'lightblue')
}

peaks <- matrix(NA, ncol = 3, nrow = nSim)
for (i in 1:nrow(toSave)) {
    peaks[i,] <- c(which.max(toSave[i,grep('Istar', colnames(toSave))]),
                   which.max(toSave[i,grep('fromI.*1\\]', colnames(toSave))]),
                   which.max(toSave[i,grep('fromH.*2\\]', colnames(toSave))]))
}
colMeans(peaks)

saveRDS(toSave, './data/sim_death.rds')

################################################################################
### equal importance alarm

# delta values must sum to less than 1
deltaC <- 0.25
deltaH <- 0.25
deltaD <- 0.25

trueVals <- c(beta = 0.85, 
              gamma1 = gamma1,
              gamma2 = gamma2,
              lambda = lambda,
              phi = phi,
              deltaC = deltaC,
              nuC = 3,
              x0C = 100,
              deltaH = deltaH,
              nuH = 3,
              x0H = 25,
              deltaD = deltaD,
              nuD = 3,
              x0D = 10)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- SIR_sim_model$expandNodeNames(dataNodes,                                                    returnScalarComponents = TRUE)

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

plot(toSave[1,grep('Istar', colnames(toSave))], type = 'l', col = 'grey', ylim = c(0, 4000))
for (i in 2:nrow(toSave)) {
    lines(toSave[i,grep('Istar', colnames(toSave))], col = 'grey')
    lines(toSave[i,grep('fromI.*1\\]', colnames(toSave))], col = 'tomato')
    lines(toSave[i,grep('fromH.*2\\]', colnames(toSave))], col = 'lightblue')
}

peaks <- matrix(NA, ncol = 3, nrow = nSim)
for (i in 1:nrow(toSave)) {
    peaks[i,] <- c(which.max(toSave[i,grep('Istar', colnames(toSave))]),
                   which.max(toSave[i,grep('fromI.*1\\]', colnames(toSave))]),
                   which.max(toSave[i,grep('fromH.*2\\]', colnames(toSave))]))
}
colMeans(peaks)

saveRDS(toSave, './data/sim_equal.rds')



