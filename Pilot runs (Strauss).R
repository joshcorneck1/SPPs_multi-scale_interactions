library(RandomFields)
library(spatstat)

## Define the observation window
height.W <- 1
width.W <- 1
mod.W <- height.W * width.W

## Define a function for the summary statistics
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
## x and y coordinate vectors
x <- seq(0,width.W,length.out = 2^8)
y <- seq(0,height.W,length.out = 2^8)

# Maxmium number of iterations
maxit <- 20000

# Minimum number of observations
m <- 40

# Number of pilot runs
kpilot <- 2000

# List of the simulated processes
Xlist <- list()
XlistA <- list()
XlistB <- list()

# Matrix for parameters
par_samp <- matrix(nrow=kpilot,ncol=3)

# Matrix for statistics
Tsamp_plus <- matrix(nrow=kpilot,ncol=47)

for (k in c(1:kpilot)){
  # The marks
  marks <- c("A","B")
  # Counts
  N1 <- 0
  N2 <- 0
  
  ## Sample parameters
  # The means
  mu <- runif(1,4,6)
  # The interaction radius
  R <- runif(1,0,0.05)
  # The gamma values
  gam <- runif(1,0,1)
  # Store the sampled parameters
  par_samp[k,] <- c(mu,R,gam)
  ## Run the algorithm (note requirement of N>= m)
  count <- 0
  while(N1 < m || N2 < m){
    ## If more than 10 runs have been done, resample 
    if (count > 10){
      ## Sample parameters
      # The mean
      mu <- runif(1,4,6)
      # The interaction radius
      R <- runif(1,0,0.05)
      # The gamma values
      gam <- runif(1,0,1)
      # Store the sampled parameters
      par_samp[k,] <- c(mu,R,gam)
    }
    R_mat <- matrix(c(0,R,R,0),2,2)
    G_mat <- matrix(gam,2,2)
    B_mat <- c(exp(mu),exp(mu))
    M <- rmhmodel(cif="straussm",par=list(beta=B_mat,gamma=G_mat,radii=R_mat),
                  types = c("A","B"),w=square(1))
    X <- rmh(M)
    N1 <- sum(X$marks=="A")
    N2 <- sum(X$marks=="B")
  }
  # Compute and save summary statistics
  Tsamp_plus[k,] <- Tx(X) 
  print(paste("k =",k))
  Xlist[[k]] <- X
}

## Save the data
save(Tsamp_plus,file = paste("summary_stats_strauss.RData",sep = ""))
save(par_samp,file = paste("params_strauss.RData",sep = ""))
save(Xlist,file=paste("Xlist_strauss.RData",sep = ""))
