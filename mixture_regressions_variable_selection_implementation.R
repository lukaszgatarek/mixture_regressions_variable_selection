library("pscl")
library("MCMCpack")
library(mvtnorm)


change_value = function(xx, new_values){
  
  # function necessary for convertion of labels,
  # accordingly to the applied reordering principle

  zz <- xx
  for (ii in 1:max(xx)){
    zz[xx == ii] <- new_values[ii]}
  return(zz)
}


## generate the observations accordingly to the data generating process implied by the model 
# number of componets in the mixture
G <- 5

# length of time series is dependent on the lengths of individual components
# vector n records the lengths of the individual series; we draw them randomly
n <- round(runif(G) * 1000, 0)
n <- rep(300, G)
# the total number of observations is implied by the length of the individual components
N <- sum(n)

## generate explanatory variables Xs for each of the G components in the mixture of regression models 

# Parameter telling how many explanatory variable there are in each equation;
# that is the maximal possible number of parameters, in practice the algorithm would be pointing at how many there should be
numXs <- 1

dimG <- rep(numXs, G)

# this artificial variable is necessary to know which parameters to choose from the vector of betas
cumDimG <- cumsum(dimG)

# We simulate N observation of the explanatory variables 
simData <- list()
for (g in 1:G){
  simData[[g]] <- matrix(rnorm(n[g] * dimG[g], 0, 1), n[g], dimG[g])
}

## define data generating process

# vector of beta parameters
betas <- (1:sum(dimG) * 2 ) + rnorm(sum(dimG),0,1)
#betas <-   rnorm(sum(dimG),0,2)

# lets keep them ordered properly
reorderWay <-  order(betas)
betas <- betas[reorderWay]
#betas <- abs(rnorm(sum(dimG), 0, 5))

# vector of sigma_g group specific parameters
sigmas <- rep(1, G) 

## simulate from the data generating process
# for every component of the regression mixture we need to select the variables which enter into the data generating process;
# define a list which will keep the indexes of those variables among 1:numXs
variable_selector = {}

# randomly select the variables that would be entering each regression component;
# the necessary indexes are sampled from multinomial distribution with equal probability (for each external variable)
for (g in 1:G){
  variable_selector[[g]] <- which(rmultinom(1, numXs, rep(1/numXs, numXs) ) > 0) }

# those indexes are applied in regression to appropriately select the variables
for (g in 1:G){
  y <- NULL
 
  print(betas[ (max(cumDimG[g-1],0) + 1) : cumDimG[g] ])
  
  for (i in 1:n[g]){
    # take appropriate parameters from the beta vector, accordingly to the values in variable selector
    y[i] <- rnorm(1, simData[[g]][i, variable_selector[[g]]] %*%
                    betas[ (max(cumDimG[g-1],0) + 1) : cumDimG[g] ][variable_selector[[g]]], sqrt(sigmas[g]) )    
  }
  simData[[g]] <- cbind(y, simData[[g]])
} 

# stack all data together and introduce a label for the true observations
z_true <- rep(NA, N, 1)
data <- NULL
for (g in 1:G){
  # stack the data component by component
  data <- rbind(data, simData[[g]])
  
  # attach the true labels to the respective observations
  z_true[((g-1) * n[g] + 1) : (g * n[g])] <- g
}

# plot join scatterplot
for (g in G:1){
    if (g==G){
      plot(data[ (((g-1) * n[g]) + 1) : (g*n[g]), 2], data[(((g-1) * n[g]) + 1) : (g*n[g]), 1], col = g) 
    } else {
      points(data[ (((g-1) * n[g]) + 1) : (g*n[g]), 2], data[(((g-1) * n[g]) + 1) : (g*n[g]), 1], col = g)   
    }
}

# plot the final dependent variable data only
for (g in G:1){
  if (g==G){
    plot( data[(((g-1) * n[g]) + 1) : (g*n[g]), 1], col = g) 
  } else {
    points( data[(((g-1) * n[g]) + 1) : (g*n[g]), 1], col = g)   
  }
}


## Priors
# priors for beta: we assume G components, and number of coefficinets accordingly to the number of parameters in the regression
muBeta <- matrix(0, G, numXs)

# priors for sigmaG
ag <- 3
bg <- 1/2
shSigma <- ag # shape parameter for inverted gamma
raSigma <- bg^(-1)

## initialization of the mcmc
# number of iterations in the MCMC chain
nSim <- 3000

# container for sampled beta parameters
mBeta <- array(NA, dim = c(nSim, numXs, G))

# draw the first beta for every group
for (g in 1:G){
  #mBeta[1,g] <-  rnorm(1, muBeta[g, 1], VBeta[g,g])
  mBeta[1, ,g] <-  rmvnorm(1, (muBeta[g, ]), VBeta)
}

## probability of inclusion for each variable
probs <- array(NA, dim = c(nSim, numXs, G))
# at the beginning let us assume that all the variables are necessary as regressors
probs[1, , ] <- 1

probPro <- matrix(NA, numXs, 1)

# for sigma
mSigma2 <- matrix(NA, nSim, G)
# draw sigma, it is group specific
mSigma2[1,] <- rigamma(1, shSigma, raSigma)

# for pi
mPi <- matrix(NA, nSim, G)
alphaPrior <- rep(N/G, G)
# draw first value of p for every observation in the sample
mPi[1,] <- rdirichlet(1, alphaPrior)

# for zi
z <- matrix(NA, nSim, N*G)

# the group 
for (t in 1:(N*G)){
  z[1, t] <- which(rmultinom(1, 1, mPi[1,]) == 1)
}

# read in the data in the appropriate variables
X <- as.matrix(data[,2:(numXs+1)])
y <- as.matrix(data[,1])

# record the amounts of observations in each component
ng <- matrix(NA, nSim,G)
# for the first one we just assume they are equal
ng[1,] <- N * G / G

#hyperparameter for variable selection
tau <- sqrt(.0000001)
c <- sqrt(9/tau^2)
# prior probability for the value to be included in the regression
p <- 1/2

for (i in 2:nSim){
  if (i %% 100 == 0){print(i)}
  #print(i)
  
  for (t in 1:(N*G)){
    fig <- NULL
    for (g in 1:G){
      fig[g] <- dnorm(y[t,1], X[t,] %*% mBeta[i-1, , g], sqrt(mSigma2[i-1,g]) ) * mPi[i-1,g] 
    }
    if (all(fig) == 0){fig <- fig + 1/G}
    
    z[i, t] <- which( rmultinom(1, 1, fig/sum(fig) ) == 1)
  }
  

  
  for (g in 1:G){
    
    for (j in 1 : numXs){
      VBeta[j, j] <- probs[i - 1, j, g] * c^2 * tau^2 + (1 - probs[i - 1, j, g]) * tau^2
    } 
    
    DBeta <- solve( t( X[ z[i,] == g, ] ) %*% X[ z[i,] == g, ] / mSigma2[i-1,g] + solve(VBeta) )
    dBeta <- t(X[ z[i,] == g, ]) %*% y[z[i,] == g, 1] / mSigma2[i-1,g] + solve(VBeta) %*% muBeta[g,]
  
    mBeta[i,,g] <- rmvnorm(1, DBeta %*% dBeta, DBeta)
    ng[i,g] <- sum(z[i,] == g)
    
    mSigma2[i,g] <- rinvgamma(1, ng[i, g]/2 + shSigma, raSigma + 1/2 * t( y[z[i,] == g, 1] - (X[ z[i,] == g, ] %*% mBeta[i,,g]) ) %*% ( y[z[i,] == g, 1] - (X[ z[i,] == g, ] %*% mBeta[i,,g]) )  )
   # mSigma2[i,g] <- rinvgamma(1, ng[i, g]/2 + shSigma, 1/b + .5 * t(y-X %*% betas[t,]) %*% (y-X %*% betas[t,]) )
  }  
  
  reorderWay <-  order(mBeta[i,1,]) # ordering by one input only is wrong as it might be actually excluded from the explanatory variables, 
  # and then the parameter would be zero..., where in fact the parameter of the real betas was different than zero
  reorderWay <-  order(apply(mBeta[i,,],2,sum) )
  mBeta[i, ,] <- mBeta[i, ,reorderWay ]
  #probs[i, ,] <- probs[i, ,reorderWay ]
  ng[i,] <- ng[i, reorderWay]
  mSigma2[i,] <- mSigma2[i, reorderWay]
  
  mPi[i,] <- rdirichlet(1, alphaPrior + ng[i,])  

  # change numbers in z accordingly to the proper labeling
  z[i, ] <- change_value(z[i, ], reorderWay)
  
  
  
  # check the inclusion probability for each variable in every component of the mixture
  for (g in 1:G){  
    for (j in 1:numXs){
      
      if (p * dnorm(mBeta[i, j, g], 0, sqrt(c^2 * tau^2) ) == 0){
        probPro[j] <- 0  
      } else {
        probPro[j] <- ( p * dnorm(mBeta[i, j, g], 0, sqrt(c^2 * tau^2) ) 
                        / (p * dnorm(mBeta[i, j, g], 0, sqrt(c^2 * tau^2) ) + (1 - p) * dnorm(mBeta[i, j, g], 0, sqrt(tau^2) ) ) )
      }
      
      probs[i, j, g] <- rbinom(1, 1, probPro[j])
    }
  }
  
    
}

# print the variable inclusion probability
for (g in 1:G){ 
  print(  apply(probs[,,g], 2, mean) ) }

# print the variables' indexes which has been used in the regressions
print(variable_selector)

# print the average sampled betas
for (g in 1:G){ 
  print(  apply(mBeta[,,g], 2, mean) ) }
s
betasMat <- t(matrix(betas, numXs, G))
for (g in 1:G){
  # mask the initially simulated values with zero
  
  betasMat[g, setdiff((1: numXs), variable_selector[[g]]) ] <- 0
}

reorderWayBetas <-  order(apply(betasMat,1,sum) )

betasTrueReordered <- vec(t(betasMat[ reorderWayBetas, ]))
print(betasTrueReordered)

# apart from that the true labels need to be reordered
z_true_reordered <- change_value(z_true, reorderWayBetas)

# check how good it fits 
par(mfrow=c(5,2))
betas_counter = 1
for (g in 1:G){
  for (v in 1:numXs){
    plot(mBeta[, v, g])
    abline( h = betasTrueReordered[betas_counter], col = 'red') 
    betas_counter <- betas_counter + 1
  }
}


par(mfrow=c(5,2))
betas_counter = 1
for (g in 1:G){
  for (v in 1:numXs){
    plot(probs[, v, g])
    betas_counter <- betas_counter + 1
  }
}


## plot the true labeling vs the expected one
# based on the last draw from the sampler
plot(z[nSim, ], pch=40, main = 'The true vs estimated labeling', ylab = 'labels')
points(z_true_reordered, col = 'red')

# based on the mode of the last 500 samples from Gibbs
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

estimated_z <-apply(z[(nSim-500+1):nSim, ], 2, getmode)

plot(estimated_z, pch=40, main = 'The true vs estimated labeling', ylab = 'labels')
points(z_true_reordered, col = 'red')

sum(z_true_reordered != estimated_z) / length(estimated_z)
