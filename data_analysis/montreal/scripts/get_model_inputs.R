################################################################################
# Get model inputs for nimble depending on which model should be fit
# used to run original model
# used for posterior prediction
# used for WAIC
################################################################################


getModelInput <- function(incData, modelType, assumeType, peak,
                          smoothC, smoothD,
                          hospData, deathData,
                          N, S0, I0, H0, D0, R0) {
    
    maxInf <- 10
    tau <- length(incData)
    
    #### initial conditions probability
    if (modelType %in% c('SIHRD_full', 'SIHRD_noAlarm')) {
        
        # SIHRD model
        initProb <- c(S0, I0, H0, R0, D0)/N
        comp_init <- rmulti(1, N, initProb)
        
        #  make sure I[1] has enough people to cover Hstar[1] and RstarI[1]
        # for undetected model only (otherwise I is known)
        if (assumeType == 'undetected') {
            
            if (I0 < H0) {
                initProb <- c(S0, I0 + H0, H0, R0, D0)/N
                comp_init <- rmulti(1, N, initProb)
            }
            
            
        }
        
    } else if (modelType %in% c('SIR_full', 'SIR_inc', 'SIR_noAlarm')) {
        
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
    
    
    # prior for probDetect depends on peak 
    #   (https://www.healthdata.org/sites/default/files/covid_briefs/101_briefing_Canada.pdf)
    # wave 1: Feb 25 - 11 July 2020          25% detection
    # wave 2: Aug 23, 2020 - March 20, 2021  40% detection
    # wave 3: March 21 - July 17, 2021       
    # wave 4: July 18 - Dec 4, 2021          20% detection
    if (peak == 1) {
        detectA <- 250
        detectB <- 750
    } else if (peak == 2) {
        detectA <- 400
        detectB <- 600
    } else if (peak == 4) {
        detectA <- 200
        detectB <- 800
    } 
    
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
                          Sigma = Sigma,
                          detectA = detectA,
                          detectB = detectB)
    
    ### data
    dataList <- list(detectIstar = incData,
                     smoothC = smoothC,
                     smoothD = smoothD,
                     constrain_deltas = 1)
    
    if (modelType == 'SIR_full') {
        
        repeat {
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1, constantsList$detectA, constantsList$detectB),
                              beta = runif(1, 1/7, 1),
                              nuC = rinvgamma(1, 7, 26),
                              nuD = rinvgamma(1, 7, 26),
                              x0C = runif(1, minC, maxC),
                              x0D = runif(1, minD, maxD),
                              Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                              w0 = rnorm(1, 5, 0.5),
                              k = rgamma(1, 100, 100))
            
            initsList$Istar <- round(incData / initsList$probDetect)
            
            delta <- multiBeta(initsList$Z)
            
            if (sum(delta) <= 1) break
        }
        
        # end modeltype == 'SIR_full'
        
    } else if (modelType == 'SIHRD_full') {
        
        dataList$Hstar <- hospData
        dataList$Dstar <- deathData
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1, constantsList$detectA, constantsList$detectB),
                              beta = runif(1, 1/7, 1),
                              gamma1 = rgamma(1, 2, 10), # IR
                              gamma2 = rgamma(1, 2, 10), # HR
                              lambda = rgamma(1, 1, 10), # IH
                              phi = rgamma(1, 1, 10) ,   # HD
                              nuC = rinvgamma(1, 7, 26),
                              nuD = rinvgamma(1, 7, 26),
                              x0C = runif(1, minC, maxC),
                              x0D = runif(1, minD, maxD),
                              Z = rmnorm_chol(1, rep(0, 2), chol(Sigma), prec_param = FALSE),
                              RstarI = round(0.3 * c(rep(0, 3), I0, dataList$detectIstar[1:(tau-4)])),
                              RstarH = round(0.01 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
            
            
            initsList$Istar <- round(incData / initsList$probDetect)
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            
            delta <- multiBeta(initsList$Z)
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1) & sum(delta) <= 1) break
            
        } 
        
        names(initsList$RstarI) <- paste0('RstarI[', 1:tau, ']')
        names(initsList$RstarH) <- paste0('RstarH[', 1:tau, ']')
        
        # end modeltype == 'SIHRD_full'
        
    } else if (modelType == 'SIR_inc') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              initProb = initProb,
                              minC = minC,
                              maxC = maxC,
                              n = n,
                              xC = xC,
                              maxInf = maxInf,
                              detectA = detectA,
                              detectB = detectB)
        
        ### data
        dataList <- list(detectIstar = incData,
                         smoothC = smoothC)
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          probDetect = rbeta(1, constantsList$detectA, constantsList$detectB),
                          beta = runif(1, 1/7, 1),
                          nuC = runif(1, 1, 10),
                          x0C = runif(1, minC, maxC),
                          deltaC = runif(1, 0, 1),
                          w0 = rnorm(1, 3, 0.5),
                          k = rgamma(1, 100, 100))
        
        initsList$Istar <- round(incData / initsList$probDetect)
        
        
        # end modeltype == 'SIR_inc'
        
    } else if (modelType == 'SIHRD_noAlarm') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              initProb = initProb,
                              maxInf = maxInf,
                              detectA = detectA,
                              detectB = detectB)
        
        ### data
        dataList <- list(detectIstar = incData,
                         Hstar = hospData,
                         Dstar = deathData)
        
        repeat {
            
            ### inits 
            initsList <- list(comp_init = comp_init,
                              probDetect = rbeta(1, constantsList$detectA, constantsList$detectB),
                              beta = runif(1, 1/7, 1),
                              gamma1 = rgamma(1, 2, 10), # IR
                              gamma2 = rgamma(1, 2, 10), # HR
                              lambda = rgamma(1, 1, 10), # IH
                              phi = rgamma(1, 1, 10) ,   # HD
                              RstarI = round(0.3 * c(rep(0, 3), I0, dataList$detectIstar[1:(tau-4)])),
                              RstarH = round(0.01 * c(rep(0, 4), dataList$Hstar[1:(tau-4)])))
            
            initsList$Istar <- round(incData / initsList$probDetect)
            
            probIH <- 1 - exp(-initsList$lambda)
            probIR <- 1 - exp(-initsList$gamma1)
            
            probHR <- 1 - exp(-initsList$gamma2)
            probHD <- 1 - exp(-initsList$phi)
            
            if ((probIH + probIR < 1) & (probHR + probHD < 1)) break
        }
        
        # end modeltype == 'SIHRD_noAlarm'
        
    } else if (modelType == 'SIR_noAlarm') {
        
        ### constants
        constantsList <- list(tau = tau,
                              N = N,
                              initProb = initProb,
                              maxInf = maxInf,
                              detectA = detectA,
                              detectB = detectB)
        
        ### data
        dataList <- list(detectIstar = incData)
        
        ### inits 
        initsList <- list(comp_init = comp_init,
                          probDetect = rbeta(1, constantsList$detectA, constantsList$detectB),
                          beta = runif(1, 1/7, 1),
                          w0 = rnorm(1, 5, 0.5),
                          k = rgamma(1, 100, 100))
        
        initsList$Istar <- round(incData / initsList$probDetect)
        
        # end modeltype == 'SIR_noAlarm'
        
    }
    
    
    if (assumeType == 'casesOnly') {
        constantsList$detectA <- NULL
        constantsList$detectB <- NULL
        initsList$probDetect <- NULL
        initsList$Istar <- NULL
        
        # change data
        dataList$detectIstar <- NULL
        dataList$Istar <- incData
        
    }
   
    
    ### MCMC specifications
    niter <- 800000
    nburn <- 400000
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
