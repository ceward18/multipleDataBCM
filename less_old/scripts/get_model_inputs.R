################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, alarmFit, alarmBase, 
                          smoothC, smoothH, smoothD, N, I0, R0) {
    
    # constants that are the same for all models
    S0 <- N - I0 - R0
    tau <- length(incData)
    maxInf <- 10
    
    ### initial conditions probability
    initProb <- c(S0, I0, N - S0 - I0)/N
    SIR_init <- rmulti(1, N, initProb)
    
    if (modelType == 'uni') {
        
        xC <- xH <- xD <- grid <- NULL
        
        alarmBasis <- switch(alarmBase,
                             'inc' = smoothC,
                             'hosp' = smoothH,
                             'death' = smoothD)
        
        ### constants
        n <- 10
        xAlarm <- seq(0, ceiling(max(alarmBasis)), length.out = n)
        
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
                         alarmBasis = alarmBasis)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          l = rinvgamma(1, vals[1], vals[2]),
                          sigma = rgamma(1, 100, 50),
                          w0 = rnorm(1, 2, 0.1),
                          k = rgamma(1, 100, 100),
                          SIR_init = SIR_init)
        
        
        # end modelType == 'uni'
        
    } else if (modelType == 'multi') {
        
        xAlarm <- NULL
        
        if (alarmFit == 'hill') {
            
            ### constants
            n <- 10
            
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
                              nu1 = runif(1, 1, 5),
                              nu2 = runif(1, 1, 5),
                              nu3 = runif(1, 1, 5),
                              gamma1 = runif(1, 0, maxC/3),
                              gamma2 = runif(1, 0, maxH/3),
                              gamma3 = runif(1, 0, maxD/3),
                              delta1 = runif(1, 0, 0.3),
                              delta2 = runif(1, 0, 0.3),
                              delta3 = runif(1, 0, 0.3),
                              w0 = rnorm(1, 2, 0.1),
                              k = rgamma(1, 100, 100),
                              SIR_init = SIR_init)
            
            # end alarmFit == 'hill'
            
        } else if (alarmFit == 'spline') {
            
            ### constants
            n <- 10
            
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
            
            # end alarmFit == 'spline'
        }
        
    } # end modeltype == 'multi'
    
    
    ### MCMC specifications
    niter <- 6e5
    nburn <- 3e5
    nthin <- 10
    
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xAlarm = xAlarm,
         grid = grid,
         xC = xC,
         xH = xH,
         xD = xD)

}