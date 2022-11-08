################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(alarmFit, incData, smoothC, smoothH, smoothD,
                          infPeriod, N, I0, R0, Rstar0, lengthI) {
    
    # constants that are the same for all models
    S0 <- N - I0 - R0
    tau <- length(incData)
    
    # for exponential infectious period
    # puts 95% probability of mean infectious period between 6 and 8 days
    bb <- 1000
    aa <- 1/5*bb
    
   if (alarmFit == 'spline') {
        
        ### constants
        n <- 50
        maxC <- ceiling(max(smoothC))
        maxH <- ceiling(max(smoothH))
        maxD <- ceiling(max(smoothD))
     
        
        xC <- seq(0, maxC, length.out = n)
        xH <- seq(0, maxH, length.out = n)
        xD <- seq(0, maxD, length.out = n)
       
        nb <- 3
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              Rstar0 = Rstar0,
                              xC = xC,
                              xH = xH,
                              xD = xD,
                              n = n,
                              maxC = maxC,
                              maxH = maxH,
                              maxD = maxD,
                              lengthI = lengthI,
                              nb = nb,
                              priorWeights = rep(1, 3))
        
        ### data
        dataList <- list(Istar = incData,
                         smoothC= smoothC,
                         smoothH= smoothH,
                         smoothD= smoothD,
                         constrain_knotsC = 1,
                         constrain_knotsH = 1,
                         constrain_knotsD = 1)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          bC = rnorm(nb, 0, 4),
                          bH = rnorm(nb, 0, 4),
                          bD = rnorm(nb, 0, 4),
                          knotsC = as.vector(quantile(xC, probs = sort(runif(nb - 1, 0.1, 0.8)))),
                          knotsH = as.vector(quantile(xH, probs = sort(runif(nb - 1, 0.1, 0.8)))),
                          knotsD = as.vector(quantile(xD, probs = sort(runif(nb - 1, 0.1, 0.8)))),
                          alpha = rep(1/3, 3))
        
        ### MCMC specifications
        niter <- 1.2e6
        nburn <- 1e6
        nthin <- 10
        
    } else if (alarmFit == 'gp') {
        
        ### constants
        n <- 15
        
        xC <- c(0, quantile(smoothC, probs = seq(0, 1, length.out = n))[-1])
        
        xC <- seq(0, ceiling(max(smoothC)), length.out = n)
        xH <- seq(0, ceiling(max(smoothH)), length.out = n)
        xD <- seq(0, ceiling(max(smoothD)), length.out = n)
        
        distC <- as.matrix(dist(matrix(xC)))
        distH <- as.matrix(dist(matrix(xH)))
        distD <- as.matrix(dist(matrix(xD)))
        
        uniqueDists <- distC[lower.tri(distC)]
        midDist <- getl(max(uniqueDists))
        valsC <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                            mid = midDist)$par, 2)
        
        uniqueDists <- distH[lower.tri(distH)]
        midDist <- getl(max(uniqueDists))
        valsH <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                             mid = midDist)$par, 2)
        
        uniqueDists <- distD[lower.tri(distD)]
        midDist <- getl(max(uniqueDists))
        valsD <- round(optim(c(3, 2), myF, lower = c(2.001, 1.001), method = 'L-BFGS-B',
                             mid = midDist)$par, 2)
        
        constantsList <- list(tau = tau,
                              N = N,
                              S0 = S0,
                              I0 = I0,
                              Rstar0 = Rstar0,
                              lengthI = lengthI,
                              distC = distC,
                              distH = distH,
                              distD = distD,
                              mu0 = 1,
                              ones = logit(seq(0.0001, 0.9999, length.out= n)),
                              n = n,
                              xC = xC,
                              xH = xH,
                              xD = xD,
                              cc = valsC[1],
                              dc = valsC[2],
                              ch = valsH[1],
                              dh = valsH[2],
                              cd = valsD[1],
                              dd = valsD[2],
                              priorWeights = rep(1, 3))
        
        ### data
        dataList <- list(Istar = incData,
                         smoothC= smoothC,
                         smoothH= smoothH,
                         smoothD= smoothD)
        
        ### inits 
        initsList <- list(beta = runif(1, 1/7, 1),
                          lC = rinvgamma(1, valsC[1], valsC[2]),
                          lH = rinvgamma(1, valsH[1], valsH[2]),
                          lD = rinvgamma(1, valsD[1], valsD[2]),
                          sigmaC = rgamma(1, 100, 50),
                          sigmaH = rgamma(1, 100, 50),
                          sigmaD = rgamma(1, 100, 50),
                          alpha = rep(1/3, 3))
        
        
        ### MCMC specifications
        niter <- 1.2e6
        nburn <- 1e6
        nthin <- 10
        
    }
    
    # adjust specs if model is exponential infectious period
    if (infPeriod == 'exp') {
        
        # adjust constants
        constantsList$lengthI <- NULL
        constantsList$Rstar0 <- NULL
        constantsList$aa <- aa
        constantsList$bb <- bb
        
        # add initial value for rateI
        initsList$rateI <- rgamma(1, aa, bb)
        
        # add initial value for Rstar (everyone removed lengthI days later)
        initsList$Rstar = c(Rstar0, 
                            dataList$Istar[1:(tau-lengthI)])
        
    }
    
    
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