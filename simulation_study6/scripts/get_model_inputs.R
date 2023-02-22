################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, smoothC, smoothD,
                          hospData, deathData) {
    
    # constants that are the same for all models
    N <- 1e6
    I0 <- 5
    R0 <- 0
    S0 <- N - I0 - R0
    tau <- length(incData)
    maxInf <- 10
    
    ### initial conditions probability
    if (modelType == 'full') {
        
        # SIHRD model - start of epidemic so no H, R, D
        initProb <- c(S0, I0, 0, 0, 0)/N
        comp_init <- rmulti(1, N, initProb)
        
    } else if (modelType %in% c('simple', 'inc')) {
        
        # SIR models - start of epidemic so R[0] = 0
        initProb <- c(S0, rep(I0/maxInf, maxInf), R0)/N
        comp_init <- rmulti(1, N, initProb)
        
    }
    
    ### constants
    n <- 50
    
    minC <- floor(min(smoothC))
    minD <- floor(min(smoothD))
    maxC <- ceiling(max(smoothC))
    maxD <- ceiling(max(smoothD))
    
    xC <- seq(0, maxC, length.out = n)
    xD <- seq(0, maxD, length.out = n)
    
    rho <- -0.6
    corMat <- matrix(c(1, rho, 
                       rho, 1), ncol = 2, byrow = T)
    sds <- rep(0.7, 2)
    Sigma <- diag(sds) %*% corMat %*% diag(sds)
    
    constantsList <- list(tau = tau,
                          N = N,
                          initProb = initProb,
                          minC = minC,
                          minD = minD,
                          maxC = maxC,
                          maxD = maxD,
                          n = n,
                          xC = xC,
                          xD = xD,
                          maxInf = maxInf,
                          zeros = rep(0, 2),
                          Sigma = Sigma)
    
    ### data
    dataList <- list(detectIstar = incData,
                     smoothC = smoothC,
                     smoothD = smoothD,
                     constrain_deltas = 1)
    
    if (modelType == 'simple') {
        
        repeat {
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1, 25, 75),
                              beta = runif(1, 1/7, 1),
                              nuC = rinvgamma(1, 7, 26),
                              nuD = rinvgamma(1, 7, 26),
                              x0C = runif(1, minC + 1, maxC/3),
                              x0D = runif(1, minD + 1, maxD/3),
                              Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                              Istar = round(dataList$detectIstar * 5) + dataList$detectIstar,
                              w0 = rnorm(1, 5, 0.5),
                              k = rgamma(1, 100, 100))
            
            
            delta <- multiBeta(initsList$Z)
            
            if (sum(delta) <= 1) break
        }
        
        # end modeltype == 'simple'
        
    } else if (modelType == 'full') {
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1, 25, 75),
                              beta = runif(1, 1/7, 1),
                              gamma1 = rgamma(1, 2, 10), # IR
                              gamma2 = rgamma(1, 2, 10), # HR
                              lambda = rgamma(1, 1, 10), # IH
                              phi = rgamma(1, 1, 10) ,   # HD
                              nuC = rinvgamma(1, 7, 26),
                              nuD = rinvgamma(1, 7, 26),
                              x0C = runif(1, minC + 1, maxC/3),
                              x0D = runif(1, minD + 1, maxD/3),
                              Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                              Istar = round(dataList$detectIstar * 5) + dataList$detectIstar,
                              RstarI = round(0.3 * c(rep(0, 3), I0, dataList$detectIstar[1:(tau-4)])),
                              RstarH = round(0.3 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            
            delta <- multiBeta(initsList$Z)
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1) & sum(delta) <= 1) break
            
        } 
        
        names(initsList$RstarI) <- paste0('RstarI[', 1:tau, ']')
        names(initsList$RstarH) <- paste0('RstarH[', 1:tau, ']')
        
        # end modeltype == 'full'
        
    } else if (modelType == 'inc') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              initProb = initProb,
                              minC = minC,
                              maxC = maxC,
                              n = n,
                              xC = xC,
                              maxInf = maxInf,
                              zeros = rep(0, 2),
                              Sigma = Sigma)
        
        ### data
        dataList <- list(detectIstar = incData,
                         smoothC = smoothC)
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          probDetect = rbeta(1, 25, 75),
                          beta = runif(1, 1/7, 1),
                          nuC = runif(1, 1, 10),
                          x0C = runif(1, minC + 1, maxC/3),
                          deltaC = runif(1, 0, 1),
                          Istar = round(dataList$detectIstar * 5) + dataList$detectIstar,
                          w0 = rnorm(1, 3, 0.5),
                          k = rgamma(1, 100, 100))
        
        
        
        # end modeltype == 'inc'
        
    } 
    
    ### MCMC specifications
    niter <- 8e5
    nburn <- 4e5
    nthin <- 20
    
    
    list(constantsList = constantsList,
         dataList = dataList,
         initsList = initsList,
         niter = niter,
         nburn = nburn,
         nthin = nthin,
         xC = xC,
         xD = xD)
    
}
