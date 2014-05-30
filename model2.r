source("simulation.r")

# Model 2: Human, cow, and dog are point sources; bird and
# deer are nonpoint sources. Cells from point sources are placed into
# a single bin, which are then allowed to mix with adjacent cells. Movement matrix 
#(i.e., current) is allowed but defaults to the identity matrix (no movement)
simulate2 <- function(nbins=500, 
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
                      
                      diffusion = 0.8, # mean proportion of cells to remain in bin after diffusion
                      diffusionSD = 0.05, # standard deviation of normal distribution from which
                                          # diffusion is drawn
                      
                      dilutionPrnt = 0.9, # mean percent remaining after each time point
                      dilutionSD = 0.1 # SD of the normal distribution from which the dilution
                      # value is drawn
){
  environment <- init(nbins)
  envRecord <- list(environment)
  
  for(i in 1:time){
    inputMatrix <- cbind("cow" = pointinputExact(cowPoints, quantity=cowQuant, nbins),
                         "human" = pointinputExact(humanPoints, quantity=humanQuant, nbins),
                         "bird" = nonpointinput(birdQuant, nbins),
                         "dog" = pointinputExact(dogPoints, quantity=dogQuant, nbins),
                         "deer"= nonpointinput(deerQuant, nbins))
    
    dilutionMatrix <- dilute(dilutionPrnt, dilutionSD, nbins)
    
    movementM <- movementMatrix(nbins, movement, movementSD)
    diffusionM <- diffusionMatrix(nbins, diffusion, diffusionSD)
    
    environment <- (movementM %*% (diffusionM %*% (environment + inputMatrix))) * dilutionMatrix
    
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