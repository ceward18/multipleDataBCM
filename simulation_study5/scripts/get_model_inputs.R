################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, smoothC, smoothH, smoothD,
                          hospData, deathData) {
    
    # constants that are the same for all models
    N <- 1e6
    I0 <- 5
    R0 <- 0
    S0 <- N - I0 - R0
    tau <- length(incData)
    
    ### constants
    n <- 50
    
    minC <- floor(min(smoothC))
    minH <- floor(min(smoothH))
    minD <- floor(min(smoothD))
    maxC <- ceiling(max(smoothC))
    maxH <- ceiling(max(smoothH))
    maxD <- ceiling(max(smoothD))
    
    xC <- seq(0, maxC, length.out = n)
    xH <- seq(0, maxH, length.out = n)
    xD <- seq(0, maxD, length.out = n)
    
    # C, H, D
    corMat <- matrix(c(1, -0.5, -0.1,
                       -0.5, 1, -0.3,
                       -0.1, -0.3, 1), ncol = 3, byrow = T)
    sds <- c(0.5, 0.5, 0.5)
    Sigma <- diag(sds) %*% corMat %*% diag(sds)
    
    constantsList <- list(tau = tau,
                          N = N,
                          S0 = S0,
                          I0 = I0,
                          minC = minC,
                          minH = minH,
                          minD = minD,
                          maxC = maxC,
                          maxH = maxH,
                          maxD = maxD,
                          n = n,
                          xC = xC,
                          xH = xH,
                          xD = xD,
                          maxInf = 10,
                          zeros = rep(0, 3),
                          Sigma = Sigma)
    
    ### data
    dataList <- list(Istar = incData,
                     smoothC = smoothC,
                     smoothH = smoothH,
                     smoothD = smoothD,
                     constrain_deltas = 1)
    
    if (modelType == 'simple') {
        
        repeat {
            ### inits 
            initsList <- list(beta = runif(1, 1/7, 1),
                              nuC = runif(1, 1, 10),
                              nuH = runif(1, 1, 10),
                              nuD = runif(1, 1, 10),
                              x0C = runif(1, 0, maxC/3),
                              x0H = runif(1, 0, maxH/3),
                              x0D = runif(1, 0, maxD/3),
                              Z = rmnorm_chol(1, rep(0, 3),  chol(Sigma), prec_param = FALSE),
                              w0 = rnorm(1, 3, 0.5),
                              k = rgamma(1, 100, 100))
            
            delta <- multiBeta(initsList$Z, constantsList$Sigma)
            
            if (sum(delta) <= 1) break
        }
        
        
        
        # end modeltype == 'simple'
        
    } else if (modelType == 'full') {
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(beta = runif(1, 1/7, 1),
                              gamma1 = runif(1), # IR
                              gamma2 = runif(1), # HR
                              lambda = runif(1), # IH
                              phi = runif(1) ,   # HD
                              nuC = runif(1, 1, 5),
                              nuH = runif(1, 1, 5),
                              nuD = runif(1, 1, 5),
                              x0C = runif(1, 0, maxC/3),
                              x0H = runif(1, 0, maxH/3),
                              x0D = runif(1, 0, maxD/3),
                              Z = rmnorm_chol(1, rep(0, 3),  chol(Sigma), prec_param = FALSE),
                              RstarI = round(0.3 * c(rep(0, 3), I0, dataList$Istar[1:(tau-4)])),
                              RstarH = round(0.3 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            delta <- multiBeta(initsList$Z, constantsList$Sigma)
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1) & sum(delta) <= 1) break
            
        }
       
        
        names(initsList$RstarI) <- paste0('RstarI[', 1:tau, ']')
        names(initsList$RstarH) <- paste0('RstarH[', 1:tau, ']')
    }
    
    
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
         xC = xC,
         xH = xH,
         xD = xD)
    
}
