library(RandomFields)
library(spatstat)

# x and y coordinate vectors
x <- seq(0,1,length.out = 2^8)
y <- seq(0,1,length.out = 2^8)

# Number of pilot runs
kpilot <- 500

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
par_pilot <- matrix(nrow=kpilot,ncol=14)

# Matrix for simulations
Xpilot <- list()

for (k in c(1:kpilot)){
  count <- 0
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
  # Save the parameters
  par_pilot[k,] <- pars
  Xpilot[[k]] <- ppp(x=as.numeric(X[,1]),y=as.numeric(X[,2]),marks=as.factor(X[,3]))
  print(paste("k =",k))
}

# Save the parameters
save(par_pilot,file = paste("pilotparams.RData",sep = ""))
save(Xpilot,file = paste("pilotsims.RData",sep = ""))
