################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, alarmBase, N, I0, R0) {
    
    # constants that are the same for all models
    S0 <- N - I0 - R0
    tau <- length(incData)
    maxInf <- 10
    
    ### initial conditions probability
    initProb <- c(S0, I0, N - S0 - I0)/N
    SIR_init <- rmulti(1, N, initProb)
    
    ### constants
    n <- 10
    xAlarm <- seq(0, ceiling(max(alarmBase)), length.out = n)
    
    dists <- as.matrix(dist(matrix(xAlarm)))
    
    midDist <- getl(max(dists[lower.tri(dists)]))
    vals <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                         mid = midDist)$par, 2)
   
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          dists = dists,
                          mu0 = 1,
                          ones = logit(seq(0.0001, 0.9999, length.out = n)),
                          n = n,
                          xAlarm = xAlarm,
                          lShape = vals[1],
                          lScale = vals[2],
                          maxInf = maxInf)
    
    ### data
    dataList <- list(Istar = incData,
                     alarmBasis = alarmBase)
    
    ### inits 
    initsList <- list(beta = runif(1, 1/7, 1),
                      l = rinvgamma(1, vals[1], vals[2]),
                      sigma = rgamma(1, 100, 50),
                      w0 = rnorm(1, 2, 0.1),
                      k = rgamma(1, 100, 100),
                      SIR_init = SIR_init)
    
    
    ### MCMC specifications
    niter <- 1e6
    nburn <- 5e5
    nthin <- 10
  
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xAlarm = xAlarm)
    
    
}