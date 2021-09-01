library(RandomFields)
library(spatstat)
library(glmnet)

## Load all necessary variables
# The variables from the pre-ABC analysis
load("eps_strauss.RData")
load("vars_strauss.RData")
load("Tobs_strauss.RData")
load("coeff_lm_strauss.RData")
load("nonzero_strauss.RData")

# x and y coordinate vectors
x <- seq(0,1,length.out = 2^8)
y <- seq(0,1,length.out = 2^8)

# Summary statistics function
Tx <- function(X){
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
  L2 <- Lcross(X,correction = "iso",r=seq(0,0.2,length.out = m))
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

# Number of ABC runs
kABC <- 200

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
par_ABC <- matrix(nrow=kABC,ncol=3)

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
      # The mean
      mu <- runif(1,4,6)
      # The interaction radius
      R <- runif(1,0,0.05)
      # The gamma values
      gam <- runif(1,0,1)
      # Store the sampled parameters
      pars <- c(mu,R,gam)
      # Simulation matrices
      R_mat <- matrix(c(0,R,R,0),2,2)
      G_mat <- matrix(gam,2,2)
      B_mat <- c(exp(mu),exp(mu))
      M <- rmhmodel(cif="straussm",par=list(beta=B_mat,gamma=G_mat,radii=R_mat),
                    types = c("A","B"),w=square(1))
      X <- rmh(M)
      N1 <- sum(X$marks=="A")
      N2 <- sum(X$marks=="B")
    }
  
    # Compute the summary statistics and add an intercept
    TABC <- Tx(X) - Tobs
    
    # Compute the distance measure value
    distms <- sum(sapply(1:3, function(j) (t(as.matrix(TABC[nonzero[[j]]],nrow=1)) %*% as.matrix(coeff_lm[[j]][-1]))^2/vars[[j]]))
    
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
save(par_ABC,file = paste("ABCparams_strauss.RData",sep = ""))

