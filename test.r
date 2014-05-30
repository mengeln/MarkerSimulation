library(reshape2)
library(ggplot2)

source("model1.r")
source("model2.r")

# test simulation
test <- simulate1(nbins = 200,
                  time = 200,
                  humanPoints = 10, # which bins are the sources
                  humanQuant = 10^3, # how many cells per time point per source
                  cowPoints = 1,
                  cowQuant = 0,
                  dogPoints = 1,
                  dogQuant = 0,
                  binDist = 0.2,
                  saveAll = TRUE)


# test <- simulate2(nbins = 200,
#                  time = 200,
#                  humanPoints = 10, # which bins are the sources
#                  humanQuant = 10^3, # how many cells per time point per source
#                  cowPoints = 1,
#                  cowQuant = 0,
#                  dogPoints = 1,
#                  dogQuant = 0,
#                  diffusion = 0.8,
#                  saveAll = TRUE)


# extract human samples from all runs
human <- sapply(test[[2]], function(x)x[, 2])
humanm <- melt(human)
colnames(humanm) <- c("bin", "timepoint", "cells")

# plot the bin at the source over time, as well as another bin 5 away
ggplot(humanm[humanm$bin %in% c(10, 15), ], aes(timepoint, cells, colour=as.factor(bin))) +
  geom_line()

# Simulate sampling from the end point of the simulation
sampleEnv(test[[1]])
summaryEnv(test[[1]])
