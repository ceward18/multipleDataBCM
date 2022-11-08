# linear interpolation function to get alarm values for each observed incidence value
nim_approx <- nimbleFunction(     
    run = function(x = double(1), y = double(1), z = double(1),
                   c = double(1), 
                   xout = double(0), yout = double(0), zout = double(0)) {
        returnType(double(0))
        
        # if xout is > max(x), return the closest value
        if (xout >= max(x)) xout <- max(x) 
        if (yout >= max(y)) yout <- max(y)
        if (zout >= max(z)) zout <- max(z)
        
        # x values on either side of xout
        xPlacement <- 1 * (xout < x)
        yPlacement <- 1 * (yout < y)
        zPlacement <- 1 * (zout < z)
        
        x0 <- x[max(which(xPlacement == 0))]
        x1 <- x[min(which(xPlacement == 1))]
        
        y0 <- y[max(which(yPlacement == 0))]
        y1 <- y[min(which(yPlacement == 1))]
        
        z0 <- z[max(which(zPlacement == 0))]
        z1 <- z[min(which(zPlacement == 1))]
        
        xd <- (xout - x0) / (x1 - x0)
        yd <- (yout - y0) / (y1 - y0)
        zd <- (zout - z0) / (z1 - z0)
        
        # function values across cube
        c000 <- c[]
        
        # interpolate along x
        c00 <- c000 * (1 - xd) + c100 * xd
        c01 <- c001 * (1 - xd) + c101 * xd
        c10 <- c010 * (1 - xd) + c110 * xd
        c11 <- c011 * (1 - xd) + c111 * xd
        
        # interpolate along y
        c0 <- c00 * (1 - yd) + c10 * yd
        c1 <- c01 * (1 - yd) + c11 * yd
        
        # interpolate along z
        c <- c0 * (1 - zd) + c1 * zd
        
        return(out)
    })


x <- seq(1, 20, length.out = 10)
y <- seq(1, 50, length.out = 10)
z <- seq(1, 10, length.out = 10)

grid <- expand.grid(x, y, z)

xout <- 12

# x values on either side of xout
xPlacement <- 1 * (xout < x)

x0 <- x[max(which(xPlacement == 0))]
x1 <- x[min(which(xPlacement == 1))]




SIR_gp <-  nimbleCode({
    
    SIR_init[1:3] ~ dmulti(prob = initProb[1:3], size = N)
    S[1] <- SIR_init[1] - 1
    I[1, 1] <- SIR_init[2] + 1 # add 1 to ensure I0 > 0
    I[1, 2:maxInf] ~ dmulti(prob = Rstar0/length(Rstar0), size = sum(Rstar0))
    
    idd_curve[1:maxInf] <- logitDecay(1:maxInf, w0, k)
    
    ### rest of time points
    for(t in 1:tau) {
        
        # compute alarm
        alarm[t] <- nim_approx(xAlarm[1:n], yAlarm[1:n], alarmBasis[t])
        
        probSI[t] <- 1 - exp(- beta * (1 - alarm[t]) * 
                                 sum(idd_curve[1:maxInf] * I[t, 1:maxInf]) / N)
        
        Istar[t] ~ dbin(probSI[t], S[t])
        
        # update S and I
        S[t + 1] <- S[t] - Istar[t]
        I[t + 1, 2:maxInf] <- I[t, 1:(maxInf - 1)]  # shift current I by one day
        I[t + 1, 1] <- Istar[t]                     # add newly infectious
        
    }
    
    # estimated effective R0
    R0_update[1:(tau-maxInf)] <- get_R0(betat = beta * (1 - alarm[1:tau]), 
                                        N = N, S = S[1:tau], maxInf = maxInf,
                                        iddCurve = idd_curve[1:maxInf])
    
    yAlarm[1] <- 0
    mu[1:n] <- mu0 * ones[1:n]
    cov[1:n, 1:n] <- sqExpCov(dists[1:n, 1:n], sigma, l)
    logit(yAlarm[2:n]) ~ dmnorm(mu[2:n], cov = cov[2:n, 2:n])
    
    # priors
    beta ~ dgamma(0.1, 0.1)
    sigma ~ dgamma(150, 50)
    l ~ dinvgamma(c, d)
    w0 ~ dnorm(2, sd = 0.1)
    k ~ dgamma(100, 100)
    
})

ns(cbind(x, y, z), df = 5)


dists <- as.matrix(dist(grid))


