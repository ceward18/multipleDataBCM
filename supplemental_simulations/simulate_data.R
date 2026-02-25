################################################################################
# Simulate data under various scenarios:
# 
################################################################################

### load libraries
library(nimble)

### source model code
source('./scripts/model_code.R')

### Number of simulated epidemics
nSim <- 50

### Constants for all models
N <- 1e6
C0 <- 2
I0 <- 8
S0 <- N - I0
tau <- 40
bw <- 30
maxInf <- 15

### parameter values for all models 

gamma1 <- 0.2  # I to Recovered
gamma2 <- 0.2  # H to Recovered
lambda <- 0.03 # I to H 
phi <- 0.1     # H to D 

# 1 in 4 cases are detected
probDetect <- 0.25

# overall alarm parameter
k <- 0.00003

# beta
beta <- 1.5

################################################################################
### Set up simulator

constantsList <- list(N = N, 
                      tau = tau,
                      S0 = S0,
                      I0 = I0,
                      C0 = C0,
                      bw = bw,
                      maxInf = maxInf)

SIR_sim_model <- nimbleModel(code = SIHRD_sim,
                             constants = constantsList)
SIR_sim_model_C <- compileNimble(SIR_sim_model)

dataNodes <- c('Istar', 'detectIstar', 'C', 'I', 'Hstar', 'RstarI', 'fromH', 'R0', 'alarm')
dataNodes <- SIR_sim_model$expandNodeNames(dataNodes)
sim_R <- simulator(SIR_sim_model, dataNodes)
sim_C <- compileNimble(sim_R, project = SIR_sim_model)

################################################################################
### High incidence alarm

trueVals <- c(probDetect = probDetect,
              beta = beta, 
              gamma1 = gamma1,
              gamma2 = gamma2,
              lambda = lambda,
              phi = phi,
              k = k,
              alpha = 0.85)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- SIR_sim_model$expandNodeNames(dataNodes,
                                                   returnScalarComponents = TRUE)

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('^Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

plot(toSave[1,grep('^Istar', colnames(toSave))], type = 'l', col = 'black', 
     ylim = c(0, 100), main = 'inc')
for (i in 2:nrow(toSave)) {
    lines(toSave[i,grep('^Istar', colnames(toSave))], col = 'black')
    lines(toSave[i,grep('detectIstar', colnames(toSave))], col = 'grey')
    lines(toSave[i,grep('Hstar', colnames(toSave))], col = 'gold')
    lines(toSave[i,grep('fromH.*2\\]', colnames(toSave))], col = 'lightblue')
}

saveRDS(toSave, './data/sim_inc.rds')

################################################################################
### High deaths alarm

trueVals <- c(probDetect = probDetect,
              beta = beta, 
              gamma1 = gamma1,
              gamma2 = gamma2,
              lambda = lambda,
              phi = phi,
              k = k,
              alpha = 0.15)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- SIR_sim_model$expandNodeNames(dataNodes,        
                                                   returnScalarComponents = TRUE)


# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('^Istar', colnames(epiSims))])
toSave <- epiSims[which(nInf > 100),][1:nSim,]

plot(toSave[1,grep('^Istar', colnames(toSave))], type = 'l', col = 'black', 
     ylim = c(0, 750), main = 'deaths')
for (i in 2:nrow(toSave)) {
    lines(toSave[i,grep('^Istar', colnames(toSave))], col = 'black')
    lines(toSave[i,grep('detectIstar', colnames(toSave))], col = 'grey')
    lines(toSave[i,grep('Hstar', colnames(toSave))], col = 'gold')
    lines(toSave[i,grep('fromH.*2\\]', colnames(toSave))], col = 'lightblue')
}

saveRDS(toSave, './data/sim_death.rds')

################################################################################
### equal importance alarm

trueVals <- c(probDetect = probDetect,
              beta = beta, 
              gamma1 = gamma1,
              gamma2 = gamma2,
              lambda = lambda,
              phi = phi,
              k = k,
              alpha = 0.5)

set.seed(1)
system.time(epiSims <- sim_C$run(trueVals, 1000))
colnames(epiSims) <- SIR_sim_model$expandNodeNames(dataNodes,                                                    returnScalarComponents = TRUE)

# save first 50 simulations with at least 100 infectious
nInf <- rowSums(epiSims[,grep('^Istar', colnames(epiSims))])
toSave <- epiSims[nInf > 100,][1:nSim,]

plot(toSave[1,grep('^Istar', colnames(toSave))], type = 'l', col = 'black', 
     ylim = c(0, 2000), main = 'equal')
for (i in 2:nrow(toSave)) {
    lines(toSave[i,grep('^Istar', colnames(toSave))], col = 'black')
    lines(toSave[i,grep('detectIstar', colnames(toSave))], col = 'grey')
    lines(toSave[i,grep('Hstar', colnames(toSave))], col = 'gold')
    lines(toSave[i,grep('fromH.*2\\]', colnames(toSave))], col = 'lightblue')
}

saveRDS(toSave, './data/sim_equal.rds')



