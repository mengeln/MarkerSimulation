
# Model 1: Cells distributed by gaussian distribution centered
# a bin, with no movement between cells
pointinput <- function(sourcepoint = 50, dist=0.5,
                       quantity=10^9, bins = 100, sd=1){
  quantiles <- lapply(sourcepoint, function(s){
    quantileVec <- dist * (1:bins - s)
    distribution <- dnorm(quantileVec, sd=sd)
    rmultinom(1, quantity, distribution)[, 1]
  })  
  Reduce('+', quantiles)
}

# Model 2: Cells put into a specific bin, and then spread
# between bins via a diffusion matrix transformation
pointinputExact <- function(sourcepoint = 50, 
                            quantity = 10^6,
                            bins = 100){
  inputs <- rep(0, bins)
  inputs[sourcepoint] <- quantity
  inputs
}

# Distribute cells among all bins with a uniform distribution
nonpointinput <- function(quantity=10^3, bins = 100){
  freq <- table(sample.int(bins, quantity, replace=TRUE))
  blank <- rep(0, bins)
  names(blank) <- as.character(1:bins)
  blank[names(freq)] <- freq
  as.vector(blank)
}

# Stocastically remove cells from all bins
dilute <- function(proportion=0.95, sd=0.1, bins=100){
  dilution <- replicate(5, rnorm(bins, proportion, sd))
  dilution[dilution > 1] <- 1
  dilution
}

# Proportion of markers in each host type
markerMatrix <- matrix(c(1, 0, 0, 0, 0,
                         0, 1, 0, 0, 0,
                         0, 0, 1, 0, 0,
                         0, 0, 0, 1, 0,
                         0, 0, 0, 0, 1),
                       ncol=5)
colnames(markerMatrix) <- c("cow", "human", "bird", "dog", "deer")
rownames(markerMatrix) <- c("cow", "human", "bird", "dog", "deer")

# Proportion of cells transported to adjacent bin each timepoint
# Movement is fixed in one direction (ascending bin order)
movementMatrix <- function(bins, prcntMov, sd){
  if(is.null(prcntMov))return(diag(bins))
  movement <- rnorm(bins, prcntMov, sd)
  movement[movement > 1] <- 1
  movement[movement < 0] <- 0
  m <- diag(bins)
  m[head(which(m == 1) + 1, -1)] <- head(movement, -1)
  m[m==1] <- 1 - movement
  m
}

# random diffusion of cells across bins used in model 2
diffusionMatrix <- function(bins, diffusion, sd){
  if(is.null(diffusion))return(diag(bins))
  remain <- rnorm(bins, diffusion, sd)
  remain[remain > 1] <- 1
  remain[remain < 0] <- 0
  movement <- (1 - remain)/2
  m <- diag(bins)
  m[head(which(m == 1) + 1, -1)] <- head(movement, -1)
  m[(which(m == 1) - 1)[-1]] <- tail(movement, -1)
  m[m==1] <- remain
  m
}

# Creates clean environment
init <- function(bins=100){
  m <- matrix(rep(0, 5 * bins), ncol=5)
  colnames(m) <- c("cow", "human", "bird", "dog", "deer")
  m
}

# Take results from simulation, divide environment into 
# n number of segments, and randomly pick one bin from
# each segment from which to sample
sampleEnv <- function(env, divisions = 3){

  segments <- split(as.data.frame(env),
                    cut(1:nrow(env), 3,
                        labels=FALSE))
  
  samples <- sapply(segments, function(x){
    x[sample(nrow(x), 1), ]
  })
  samples
}

# Summary of an environment
summaryEnv <- function(env){
  summ <- data.frame(mean = apply(env, 2, mean),
             median = apply(env, 2, median),
             max = apply(env, 2, max),
             sd = apply(env, 2, sd))
  summ$cov <- (summ$sd^2)/summ$mean
  summ
}

