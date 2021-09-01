setwd("C:/Users/joshu/OneDrive/Documents/Academics/Postgraduate/Imperial/Project/MSc Project Code/LGCP-Strauss ABC/Pilot runs for simple model")

## Function for summary statistics
Tx <- function(X){
  height.W <- 1
  width.W <- 1
  mod.W <- height.W * width.W
  m <- 40
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

## Combine parameters from pilot runs (last used on the Strauss pilots)
num <- c(2:5)
load("params_strauss1.RData")
params <- par_samp
for (j in num){
  load(paste("params_strauss",j,".RData",sep=""))
  params <- rbind(params,par_samp)
}
save(params, file = "params_comb_strauss.RData")

## Combine summary statistics from pilot runs
load("summary_stats_strauss1.RData")
statistics <- Tsamp_plus
for (j in num){
  load(paste("summary_stats_strauss",j,".RData",sep=""))
  statistics <- rbind(statistics,Tsamp_plus)
}
save(statistics, file = "stats_comb_strauss.RData")

## Compute the summary stats on the observed pattern and
## subtract from the pilot statistics
load("jap_pines.RData") # (last used for Japanese pines data)
Xobs <- X
Tobs <- Tx(Xobs)
Tsamp_sub <- statistics - t(replicate(10000,Tobs))

save(Tobs,file = "Tobs_strauss.RData")
save(Tsamp_sub,file="Tsamp_sub_strauss.RData")

## Fit the linear model
library(glmnet)
coeff_lm <- list()
nonzero <- list()
vars <- list()
for (i in 1:ncol(params)){
  cv.gnet <- cv.glmnet(x=Tsamp_sub,y=as.matrix(params[,i]),family="gaussian",alpha=1)
  Coeff <- coef(cv.gnet,s="lambda.min")
  # Extract lasso coefficients 
  coeff <- as.matrix(Coeff)
  # Which coefficients not equal to zero
  non0 <- which(coeff != 0)
  non0 <- non0[-1]
  if (length(non0)!=0){
    non0 <- non0 - 1
  }
  ## Fit the normal linear model
  # If the length is 0, only the intercept is significant, and so
  # we fit an intercept-only model.
  if (length(non0)==0){
    lm <- lm(as.matrix(params[,i]) ~ 1)
  }else{
    lm <- lm(as.matrix(params[,i]) ~ Tsamp_sub[,non0])
  }
  na <- which(is.na(lm$coefficients),arr.ind = T)
  if (length(na) != 0){
    lm$coefficients <- lm$coefficients[-na]
    non0 <- non0[-(na-1)]
  }
  coeff_lm[[i]] <- lm$coefficients
  nonzero[[i]] <- non0
  T_temp <- cbind(1,Tsamp_sub[,nonzero[[i]]])
  # Compute variance
  vars[[i]] <- var(T_temp%*%coeff_lm[[i]])
  print(i)
}

pos <- which(vars != 0)

# Save the relevant variables
save(vars,file="vars_strauss.RData")
save(coeff_lm,file="coeff_lm_strauss.RData")
save(nonzero,file="nonzero_strauss.RData")

# Vector of intercepts 
ints <- sapply(1:ncol(params), function(i) coeff_lm[[i]][1])

## Compute distance measure
kpilot <- 10000
distms <- seq(0, 1, length=kpilot)
for(i in 1:kpilot){
  distms[i] = sum(sapply(pos, function(j) (t(as.matrix(Tsamp_sub[i,nonzero[[j]]],nrow=1)) %*% as.matrix(coeff_lm[[j]][-1]))^2/vars[[j]]))
}

# The quantile
eps <- quantile(distms,probs=0.1)
save(eps, file="eps_strauss.RData")
