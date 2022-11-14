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
    n <- 10
    
    maxC <- ceiling(max(smoothC))
    maxH <- ceiling(max(smoothH))
    maxD <- ceiling(max(smoothD))
    
    xC <- seq(0, maxC, length.out = n)
    xH <- seq(0, maxH, length.out = n)
    xD <- seq(0, maxD, length.out = n)
    
    distC <- as.matrix(dist(matrix(xC)))
    distH <- as.matrix(dist(matrix(xH)))
    distD <- as.matrix(dist(matrix(xD)))

    midDist <- getl(max(distC[lower.tri(distC)]))
    valsC <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                         mid = midDist)$par, 2)

    midDist <- getl(max(distH[lower.tri(distH)]))
    valsH <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                         mid = midDist)$par, 2)

    midDist <- getl(max(distD[lower.tri(distD)]))
    valsD <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                         mid = midDist)$par, 2)
    
    ### constants
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          n = n,
                          distC = distC,
                          distH = distH,
                          distD = distD,
                          mu0 = 1,
                          ones = logit(seq(0.0001, 0.9999, length.out = n)),
                          xC = xC,
                          xH = xH,
                          xD = xD,
                          maxInf = 10,
                          priorWeights = rep(1, 3),
                          lCShape = valsC[1],
                          lCScale = valsC[2],
                          lHShape = valsH[1],
                          lHScale = valsH[2],
                          lDShape = valsD[1],
                          lDScale = valsD[2])
    
    ### data
    dataList <- list(Istar = incData,
                     smoothC= smoothC,
                     smoothH= smoothH,
                     smoothD= smoothD)
    
    ### inits 
    initsList <- list(beta = runif(1, 1/7, 1),
                      alpha = rdirch(1, c(1, 1, 1)),
                      sigma = rgamma(1, 100, 50),
                      lC = rinvgamma(1, valsC[1], valsC[2]),
                      lH = rinvgamma(1, valsH[1], valsH[2]),
                      lD = rinvgamma(1, valsD[1], valsD[2]),
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
         xC = xC,
         xH = xH,
         xD = xD)
    
    
}