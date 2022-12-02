################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, smoothC, smoothH, smoothD, N, I0, R0) {
    
    # constants that are the same for all models
    S0 <- N - I0 - R0
    tau <- length(incData)
    maxInf <- 10
    
    ### initial conditions probability
    initProb <- c(S0, I0, N - S0 - I0)/N
    SIR_init <- rmulti(1, N, initProb)
    
    ### constants
    n <- 20
    
    minC <- floor(min(smoothC))
    minH <- floor(min(smoothH))
    minD <- floor(min(smoothD))
    maxC <- ceiling(max(smoothC))
    maxH <- ceiling(max(smoothH))
    maxD <- ceiling(max(smoothD))
    
    xC <- seq(0, maxC, length.out = n)
    xH <- seq(0, maxH, length.out = n)
    xD <- seq(0, maxD, length.out = n)
    
    grid <- expand.grid(xC = xC, xH = xH, xD = xD)
    nn <- nrow(grid)
    
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          nn = nn,
                          xC_grid = grid[,1],
                          xH_grid = grid[,2],
                          xD_grid = grid[,3],
                          minC = minC,
                          minH = minH,
                          minD = minD,
                          maxC = maxC,
                          maxH = maxH,
                          maxD = maxD,
                          maxInf = 10)
    
    ### data
    dataList <- list(Istar = incData,
                     smoothC= smoothC,
                     smoothH= smoothH,
                     smoothD= smoothD,
                     constrain_deltas = 1)
    
    ### inits 
    initsList <- list(beta = runif(1, 1/7, 1),
                      delta1 = runif(1, 0, 0.3),
                      delta2 = runif(1, 0, 0.3),
                      delta3 = runif(1, 0, 0.3),
                      H1 = runif(1, 0, maxC/N/3),
                      H2 = runif(1, 0, maxH/N/3),
                      H3 = runif(1, 0, maxD/N/3),
                      w0 = rnorm(1, 2, 0.1),
                      k = rgamma(1, 100, 100),
                      SIR_init = SIR_init)
    
    
    
    ### MCMC specifications
    niter <- 5e5
    nburn <- 2e5
    nthin <- 10
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         grid = grid)
    
    
}