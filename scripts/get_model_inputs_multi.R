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
    
    # distC <- as.matrix(dist(matrix(xC)))
    # distH <- as.matrix(dist(matrix(xH)))
    # distD <- as.matrix(dist(matrix(xD)))
    # 
    # midDist <- getl(max(distC[lower.tri(distC)]))
    # valsC <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
    #                      mid = midDist)$par, 2)
    # 
    # midDist <- getl(max(distH[lower.tri(distH)]))
    # valsH <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
    #                      mid = midDist)$par, 2)
    # 
    # midDist <- getl(max(distD[lower.tri(distD)]))
    # valsD <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
    #                      mid = midDist)$par, 2)
    # constantsList <- list(tau = tau,
    #                       N = N,
    #                       initProb = initProb,
    #                       n = n,
    #                       distC = distC,
    #                       distH = distH,
    #                       distD = distD,
    #                       mu0 = 1,
    #                       ones = logit(seq(0.0001, 0.9999, length.out= n)),
    #                       xC = xC,
    #                       xH = xH,
    #                       xD = xD,
    #                       maxInf = 10,
    #                       priorWeights = rep(1, 3),
    #                       cc = valsC[1],
    #                       dc = valsC[2],
    #                       ch = valsH[1],
    #                       dh = valsH[2],
    #                       cd = valsD[1],
    #                       dd = valsD[2],
    #                       zProb = c(1, 1, 1))
    
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          n = n,
                          nn = nn,
                          xC = xC,
                          xH = xH,
                          xD = xD,
                          xC_grid = grid[,1],
                          xH_grid = grid[,2],
                          xD_grid = grid[,3],
                          minC = minC,
                          minH = minH,
                          minD = minD,
                          maxC = maxC,
                          maxH = maxH,
                          maxD = maxD,
                          maxInf = 10,
                          grid = as.matrix(grid),
                          nb = 9)
    
    ### data
    dataList <- list(Istar = incData,
                     smoothC= smoothC,
                     smoothH= smoothH,
                     smoothD= smoothD,
                     constrain_min = 1,
                     constrain_max = 1)
    
    ### inits 
    
    ### inits (must satisfy constraint)
    repeat {
        initsList <- list(beta = runif(1, 1/7, 1),
                          b = runif(constantsList$nb, 0, 0.05),
                          knotsC = as.vector(quantile(xC, probs = sort(runif(2, 0.2, 0.8)))),
                          knotsH = as.vector(quantile(xH, probs = sort(runif(2, 0.2, 0.8)))),
                          knotsD = as.vector(quantile(xD, probs = sort(runif(2, 0.2, 0.8)))),
                          w0 = rnorm(1, 2, 0.1),
                          k = rgamma(1, 100, 100),
                          SIR_init = SIR_init)
        
        alarmInit <- splineAlarm(grid[,1], grid[,2], grid[,3],
                                 initsList$b, 
                                 initsList$knotsC, initsList$knotsH, initsList$knotsD)
        
        cond <- all(alarmInit >= 0) & all(alarmInit <= 1) 
        if (cond) break
    }
    
    
    
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
         grid = grid)
    
    
}