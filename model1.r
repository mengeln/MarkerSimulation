source("simulation.r")

# Model 1: Human, cow, and dog are point sources; bird and
# deer are nonpoint sources. Point sources are distributed to bin
# via normal distribution. No diffusion. Movement matrix (i.e., current) is
# allowed but defaults to the identity matrix (no movement)
simulate1 <- function(nbins=500, 
                      time=100, # number of timepoints
                      saveAll = FALSE, # retain the environment at every timepoint
                      
                      humanPoints = c(50, 400), # which bins are the sources
                      humanQuant = 10^9, # how many cells per time point per source
                      cowPoints = c(400),
                      cowQuant = 10^9,
                      dogPoints = c(250, 300),
                      dogQuant = 10^9,
                      birdQuant = 10^3, # bird and deer are non-point sources
                      deerQuant = 10^3,
                      
                      markerPrev = markerMatrix,
                      movement = NULL, # proportion of cells to move to adjacent bin each
                      # time point. If NULL, then no movement
                      movementSD = 0.01, # SD of movement
                      
                      binDist = 0.2, # distance between bins on a normal distribution with
                      # mean of 0 and sd of 1. For example, if binDist = 0.2,
                      # then 10 bins away is approximately 2 SD (1.95/0.2).
                      
                      dilutionPrnt = 0.9, # mean percent remaining after each time point
                      dilutionSD = 0.1 # SD of the normal distribution from which the dilution
                      # value is drawn
){
  environment <- init(nbins)
  envRecord <- list(environment)
  
  for(i in 1:time){
    inputMatrix <- cbind("cow" = pointinput(cowPoints, binDist, quantity=cowQuant, nbins),
                         "human" = pointinput(humanPoints, binDist, quantity=humanQuant, nbins),
                         "bird" = nonpointinput(birdQuant, nbins),
                         "dog" = pointinput(dogPoints, binDist, quantity=dogQuant, nbins),
                         "deer"= nonpointinput(deerQuant, nbins))
    
    dilutionMatrix <- dilute(dilutionPrnt, dilutionSD, nbins)
    
    movementM <- movementMatrix(nbins, movement, movementSD)
    
    environment <- (movementM %*% environment + inputMatrix) * dilutionMatrix
    
    if(saveAll)
      envRecord <-  c(list(environment), envRecord)
  }
  markerCellMatrix <- environment %*% markerPrev
  
  if(saveAll){
    envRecord <- lapply(envRecord, '%*%', markerPrev)
    list(markerCellMatrix, rev(envRecord))
  }
  else
    markerCellMatrix
}