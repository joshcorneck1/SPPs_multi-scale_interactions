library(RandomFields)
library(spatstat)
library(glmnet)

## Load all necessary variables
# These are the variables from the the pre-ABC analysis
load("eps_simple.RData")
load("vars_simple.RData")
load("coeff_lm_simple.RData")
load("nonzero_simple.RData")
load("Tobs_simple.RData")

## Define the observation window
height.W <- 0.5
width.W <- 1
mod.W <- height.W * width.W

## Define a function for the summary statistics
Tx <- function(X){
  X <- ppp(x=as.numeric(X[,1]),y=as.numeric(X[,2]),marks = as.factor(X[,3]),window = owin(x=c(0,width.W),y=c(0,height.W)))
  #1 
  n <- X$n
  L <- Lcross(X,correction = "iso")
  #2
  Lmax <- max(L$iso - L$r)
  #3
  Lmin <- min(L$iso - L$r)
  #4 
  Largmin <- L$r[which((L$iso - L$r) == max(L$iso - L$r))]
  #5
  L2 <- Lcross(X,correction = "iso",r=seq(0,0.2*height.W,length.out = m))
  Lspace <- L2$iso - L2$r
  #6
  q <- 5
  count <- quadratcount(X,nx=5,ny=5)
  Cmax <- max(count)/n
  #7
  Cmin <- min(count)/n
  #8
  countlist <- as.vector(count)/n
  Clogvar <- log(var(countlist))
  return(c(n,Lmax,Lmin,Largmin,Lspace,Cmax,Cmin,Clogvar))
}

## x and y coordinate vectors (less points on y since half the length)
x <- seq(0,width.W,length.out = 2^8)
y <- seq(0,height.W,length.out = 2^7)

# Number of ABC runs
kABC <- 50

# Maximum number of iterations
maxit <- 20000

# Minimum number of observations
m <- 40

# Set the counts to 0
N1 <- 0
N2 <- 0

# The marks
marks <- c("A","B")

# Matrix for parameters
par_ABC <- matrix(nrow=kABC,ncol=5)

for (k in c(1:kABC)){
  # Acceptance
  ACCEPT <- 0
  count <- 0
  while(ACCEPT == 0){
    print(paste("count=",count))
    # Counts
    N1 <- 0
    N2 <- 0
    ## Sample parameters
    # The means
    ## Sample parameters
    # The means
    mu1 <- runif(1,5,7)
    mu2 <- mu1
    mu <- mu1
    # The scales for the covariance functions
    s <- runif(1,0.01,0.3)
    # The variance
    sig2 <- runif(1,0,4)
    # The interaction radius
    R <- runif(1,0,0.05)
    # The gamma values
    gam1 <- 1; gam2 <- 1; gam <- runif(1,0,1) ;gam12 <- c(gam1,gam2)
    # Store the sampled parameters
    pars <- c(mu,s,sig2,R,gam)
    count_sims <- 0
    ## Run the algorithm (note requirement of N>= m)
    while(N1 < m || N2 < m){
      ## If more than 10 runs have been done, resample
      if (count_sims > 10){
        ## Sample parameters
        # The means
        mu1 <- runif(1,5,7)
        mu2 <- mu1
        mu <- mu1
        # The scales for the covariance functions
        s <- runif(1,0.01,0.3)
        # The variance
        sig2 <- runif(1,0,4)
        # The interaction radius
        R <- runif(1,0,0.05)
        # The gamma values
        gam1 <- 1; gam2 <- 1; gam <- runif(1,0,1); gam12 <- c(gam1,gam2) #runif(1,0,1)
        # Store the sampled parameters
        pars <- c(mu,s,sig2,R,gam)
      }
      ## Sample the field
      # Use Ai to denote ith column of A
      A1 <- c(1,1)
      # If A has >1 column, add e.g. + RMmatrix(M=A2,RMexp(var=,scale=))
      model <- RMmatrix(M = A1, RMexp(var=sig2,scale=s))
      simu <- RFsimulate(model,x,y)
      # It has zero mean, so simply add vector of means to the values
      simu$variable1 <- simu$variable1 + mu1
      simu$variable2 <- simu$variable2 + mu2
      # Matrices for the field
      z1 <- exp(RFspDataFrame2conventional(simu)$data[,,1])
      z2 <- exp(RFspDataFrame2conventional(simu)$data[,,2])
      z <- list(z1,z2)
      
      ## Run the procedure
      # Initialise with the point that maxmimises the field
      mark <- sample(marks,1)
      rm <- ifelse(mark=="A",1,2)
      a <- which(z[[rm]] == max(z[[rm]]), arr.ind = TRUE)
      X <- matrix(c(x[a[1]],y[a[2]],mark),nrow=1,ncol=3)
      
      # Run the simulation
      for (i in 1:maxit){
        # Sample mark
        mark <- sample(marks,1)
        # Integer equivalent
        rm <- ifelse(mark=="A",1,2)
        p <- runif(1,0,1)
        if (p <= 0.5){
          rx <- sample(c(1:length(x)),size=1)
          ry <- sample(c(1:length(y)),size=1)
          # Sampled point (and mark)
          u <- c(x[rx],y[ry],mark)
          ## Points with:
          # Different mark
          diff_mark <- which(X[,3] != mark)
          # Same mark
          same_mark <- which(X[,3] == mark)
          ## Compute number of near-points of:
          # Different mark
          lR_diff <- sum(unlist(sapply(1:length(diff_mark), function(j) dist(rbind(X[diff_mark[j],1:2],u[1:2])) < R)))
          # Same mark
          lR_same <- sum(unlist(sapply(1:length(same_mark), function(j) dist(rbind(X[same_mark[j],1:2],u[1:2])) < R)))
          # Fix any NAs
          if (is.na(lR_diff)){
            lR_diff <- 0
          }
          if (is.na(lR_same)){
            lR_same <- 0
          }
          # Compute the acceptance ratio
          A <- z[[rm]][rx,ry]*(gam12[rm]^lR_same)*(gam^lR_diff)/(2*(nrow(X)+1))*mod.W
          if (runif(1,0,1) < A){
            X <- rbind(X,u[1:3])
          }
        }
        if (p > 0.5){
          r <- sample(1:nrow(X),1)
          u <- X[r,] 
          rx <- which(x==u[1])
          ry <- which(y==u[2])
          diff_mark <- which(X[,3] != u[3])
          same_mark <- which(X[,3] == u[3])
          lR_diff <- sum(unlist(sapply(1:length(diff_mark), function(j) dist(rbind(X[diff_mark[j],1:2],u[1:2])) < R)))
          lR_same <- sum(unlist(sapply(1:length(same_mark), function(j) dist(rbind(X[same_mark[j],1:2],u[1:2])) < R))) - 1
          if (is.na(lR_diff)){
            lR_diff <- 0
          }
          if (is.na(lR_same)){
            lR_same <- 0
          }
          A <- (2*nrow(X)/mod.W)*(z[[rm]][rx,ry]*(gam12[rm]^lR_same)*(gam^lR_diff))^(-1)
          if (runif(1,0,1) < A && nrow(X) > 1){
            X <- matrix(X[-r,],ncol=3)
          }
        }
      }
      N1 <- nrow(X[which(X[,3]=="A"),])
      N2 <- nrow(X[which(X[,3]=="B"),])
      count_sims <- count_sims + 1
    }
  
    # Compute the summary statistics and add an intercept
    TABC <- Tx(X) - Tobs
    
    # Compute the distance measure value
    distms <- sum(sapply(1:5, function(j) (t(as.matrix(TABC[nonzero[[j]]],nrow=1)) %*% as.matrix(coeff_lm[[j]][-1]))^2/vars[[j]]))

    # Check if we accept the sample
    if (distms < eps){
      ACCEPT <- 1
    }
    count <- count + 1
    print(distms)
    print(pars)
  }
  # Save the parameters
  par_ABC[k,] <- pars
  print(paste("k =",k))
}

# Save the parameters
save(par_ABC,file = paste("ABCparams.RData",sep = ""))

