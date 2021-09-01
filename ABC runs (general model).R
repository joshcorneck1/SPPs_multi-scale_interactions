library(RandomFields)
library(spatstat)
library(glmnet)

## Load all necessary variables
# These are the variables from the pre-ABC analysis
load("eps_gen.RData")
load("params_gen.RData")
load("vars_gen.RData")
load("Tobs_gen.RData")
load("coeff_lm_gen.RData")
load("nonzero_gen.RData")

# x and y coordinate vectors
x <- seq(0,1,length.out = 2^8)
y <- seq(0,1,length.out = 2^8)

# Function to compute summary statistics
Tx <- function(X){
  m <- 40
  if (is.ppp(X) == F){
    X <- ppp(x=as.numeric(X[,1]),y=as.numeric(X[,2]),marks = as.factor(X[,3]))
  }
  Xa <- split(X)$A
  Xb <- split(X)$B
  L <- Lcross(X,correction = "iso")
  La <- Lest(Xa,correction = "iso")
  Lb <- Lest(Xb,correction="iso")
  L2 <- Lcross(X,correction = "iso",r=seq(0,0.2,length.out = m))
  L2a <- Lest(Xa,correction = "iso",r=seq(0,0.2,length.out = m))
  L2b <- Lest(Xb,correction = "iso",r=seq(0,0.2,length.out = m))
  q <- 5
  count <- quadratcount(Xa,nx=5,ny=5)
  counta <- quadratcount(Xa,nx=5,ny=5)
  countb <- quadratcount(Xb,nx=5,ny=5)
  #1 
  n <- X$n
  na <- Xa$n
  nb <- Xb$n
  #2
  Lmax <- max(L$iso - L$r)
  Lmaxa <- max(La$iso - La$r)
  Lmaxb <- max(Lb$iso - Lb$r)
  #3
  Lmin <- min(L$iso - L$r)
  Lmina <- min(La$iso - La$r)
  Lminb <- min(Lb$iso - Lb$r)
  #4 
  Largmin <- L$r[which((L$iso - L$r) == max(L$iso - L$r))]
  Largmina <- La$r[which((La$iso - La$r) == max(La$iso - La$r))]
  Largminb <- Lb$r[which((Lb$iso - Lb$r) == max(Lb$iso - Lb$r))]
  #5
  Lspace <- L2$iso - L2$r
  Lspacea <- L2a$iso - L2a$r
  Lspaceb <- L2b$iso - L2b$r
  #6
  Cmax <- max(count)/n
  Cmaxa <- max(counta)/na
  Cmaxb <- max(countb)/nb
  #7
  Cmin <- min(count)/n
  Cmina <- min(counta)/na
  Cminb <- min(countb)/nb
  #8
  countlist <- as.vector(count)/n
  countlista <- as.vector(counta)/na
  countlistb <- as.vector(countb)/nb
  Clogvar <- log(var(countlist))
  Clogvara <- log(var(countlista))
  Clogvarb <- log(var(countlistb))
  #9
  nsign <- sum(diff(sign(L$iso - L$r)) != 0) - 1
  return(c(n,na,nb,Lmax,Lmaxa,Lmaxb,Lmin,Lmina,Lminb,Largmin,
           Largmina,Largminb,Lspace,Lspacea,Lspaceb,Cmax,Cmaxa,
           Cmaxb,Cmin,Cmina,Cminb,Clogvar,Clogvara,Clogvarb,nsign))
}

# Number of ABC runs
kABC <- 40

# Maximum number of iterations
maxit <- 20000

# Minimum number of observations
m <- 40

# Set the counts of each subprocess to 0
N1 <- 0
N2 <- 0

# The marks
marks <- c("A","B")

# Matrix for parameters
par_ABC_gen <- matrix(nrow=kABC,ncol=ncol(params_gen))

for (k in c(1:kABC)){
  # Acceptance
  ACCEPT <- 0
  count <- 0
  while(ACCEPT == 0){
    print(paste("count=",count))
    # Counts
    N1 <- 0
    N2 <- 0
    ## Run the algorithm (note requirement of N>= m)
    while(N1 < m || N2 < m){
      ## Sample parameters
      # The means
      mu1 <- runif(1,4,6)
      mu2 <- runif(1,4,6)
      
      # The scales for the covariance functions
      s1 <- runif(1,0.01,0.5)
      s2 <- runif(1,0.01,0.5)
      # The variance
      sig21 <- 1
      sig22 <- 1
      # The interaction radius
      R1 <- runif(1,0,0.05)
      R2 <- runif(1,0,0.05)
      R <- runif(1,0,0.05)
      # The gamma values
      gam1 <- runif(1,0.5,1); gam2 <- runif(1,0.5,1); gam <- runif(1,0,1); gam12 <- c(gam1,gam2)
      A1 <- runif(2,-2,2)
      A2 <- runif(2,-2,2)
      # Store the sampled parameters
      pars <- c(mu1,mu2,s1,s2,R1,R2,R,gam1,gam2,gam,A1,A2)
      ## Sample the field
      model <- RMmatrix(M = A1, RMexp(var=sig21,scale=s1)) + RMmatrix(M = A2, RMexp(var=sig22,scale=s2))
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
      X <- matrix(c(x[a[2]],y[a[1]],mark),nrow=1,ncol=3)
      
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
          if (rm==1){
            lR_same <- sum(unlist(sapply(1:length(same_mark), function(j) dist(rbind(X[same_mark[j],1:2],u[1:2])) < R1)))
          }else{
            lR_same <- sum(unlist(sapply(1:length(same_mark), function(j) dist(rbind(X[same_mark[j],1:2],u[1:2])) < R2)))
          }
          # Fix any NAs
          if (is.na(lR_diff)){
            lR_diff <- 0
          }
          if (is.na(lR_same)){
            lR_same <- 0
          }
          # Compute the acceptance ratio
          A <- z[[rm]][rx,ry]*(gam12[rm]^lR_same)*(gam^lR_diff)/(2*(nrow(X)+1))
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
          if (rm==1){
            lR_same <- sum(unlist(sapply(1:length(same_mark), function(j) dist(rbind(X[same_mark[j],1:2],u[1:2])) < R1))) - 1
          }else{
            lR_same <- sum(unlist(sapply(1:length(same_mark), function(j) dist(rbind(X[same_mark[j],1:2],u[1:2])) < R2))) - 1
          }
          if (is.na(lR_diff)){
            lR_diff <- 0
          }
          if (is.na(lR_same)){
            lR_same <- 0
          }
          A <- 2*nrow(X)*(z[[rm]][rx,ry]*(gam12[rm]^lR_same)*(gam^lR_diff))^(-1)
          if (runif(1,0,1) < A && nrow(X) > 1){
            X <- matrix(X[-r,],ncol=3)
          }
        }
      }
      N1 <- nrow(X[which(X[,3]=="A"),])
      N2 <- nrow(X[which(X[,3]=="B"),])
    }
  
    # Compute the summary statistics and add an intercept
    TABC <- Tx(X) - Tobs_gen
    
    # Compute the distance measure value
    distms <- sum(sapply(pos, function(j) (t(as.matrix(TABC[nonzero_gen[[j]]],nrow=1)) %*% as.matrix(coeff_lm_gen[[j]][-1]))^2/vars_gen[[j]]))
    
    # Check if we accept the sample
    if (distms < eps_gen){
      ACCEPT <- 1
    }
    count <- count + 1
  }
  # Save the parameters
  par_ABC_gen[k,] <- pars
  print(paste("k =",k))
}

# Save the parameters
save(par_ABC_gen,file = paste("ABCparams_gen",paste(run),".RData",sep = ""))

