# load simulation output for median parameter values and parameter envelopes, and plot summary statistic envelopes together with observations

library(ggplot2)
library(gridExtra)

rm(list=ls())

setwd("put your working directory here")
load(file = 'miceTicksBbBab_ABCRej0005_Envelopes30_BI.RData')
xStarNew50 <- xStarNew  # output of trials with parameters values from the 30% credible intervals of the posterior distributions
load(file = 'miceTicksBbBab_ABCRej0005_Envelopes10_BI.RData')
xStarNew10 <- xStarNew  # output of trials with parameters values from the 10% credible intervals of the posterior distributions

# function to plot envelopes based on quantiles together with observations
plot_envelopes = function(obsDf, df, varName) {
  ggplot() + theme_grey(base_size = 22) +
    geom_ribbon(data = df, aes(x = day, ymin = xStarMin10, ymax = xStarMax10), fill="gray30", alpha=0.5) + 
    geom_ribbon(data = df, aes(x = day, ymin = xStarMin30, ymax = xStarMax30), fill="gray30", alpha=0.2) + 
    geom_ribbon(data = df, aes(x = day+365, ymin = xStarMin10, ymax = xStarMax10), fill="gray30", alpha=0.5) + 
    geom_ribbon(data = df, aes(x = day+365, ymin = xStarMin30, ymax = xStarMax30), fill="gray30", alpha=0.2) + 
    geom_ribbon(data = df, aes(x = day+730, ymin = xStarMin10, ymax = xStarMax10), fill="gray30", alpha=0.5) + 
    geom_ribbon(data = df, aes(x = day+730, ymin = xStarMin30, ymax = xStarMax30), fill="gray30", alpha=0.2) + 
    
    geom_line(data = df, aes(x = day, y = xStarMedian), color="skyblue2", size = 1, alpha=1) +
    geom_line(data = df, aes(x = day+365, y = xStarMedian), color="skyblue2", size = 1, alpha=1) +
    geom_line(data = df, aes(x = day+730, y = xStarMedian), color="skyblue2", size = 1, alpha=1) +
    
    geom_point(data = obsDf, aes(x = obsDays, y = obs), color="salmon", size = 2) +
    xlab('Day') + ylab(varName) + xlim(90,1000)
}

# -------------- main body ---------------------------

# construct aggregated variables to match data
xStarMedian; # simulation output with median parameter values
xStarMin10 = c(); xStarMax10 = c() # max and min values across all trials of each state variables at each time point 
xStarMin30 = c(); xStarMax30 = c()

threshViable <- 1e-2  # minimum number of mouse infections required for a simulation to be considered 'viable' i.e. R0 > 1 for both pathogens

# generate trials from parameters sampled from the 10% (cut 2) and 30% (cut 1) credible intervals of the posteriors
for (iCut in 1:2) {
  if (iCut == 1) {xStarNew <- xStarNew30}
  else {xStarNew <- xStarNew10}
  
  xStarSummary = matrix(NA, nrow = nrow(xStarNew), ncol = 217) # summary stats produced by combining some state variables
  xStarViable = matrix(FALSE, nrow = nrow(xStarNew), ncol = 1) # indicates whether a trial was viable, see above
  
  for (iTrial in 1:nrow(xStarNew)) {
      # for each trial, mouse density is the sum of elements 1-4
      denMice <- apply( matrix(xStarNew[iTrial,], 31, 15)[,1:4], 1, sum)
      # mouse Bb prevalence is sum of elements 2&4, over sum of elements 1-4
      prevBbMice <- apply(matrix(xStarNew[iTrial,], 31, 15)[, c(2,4)], 1, sum)/denMice
      # mouse Bb prevalence is sum of elements 3&4, over sum of elements 1-4
      prevBabMice <- apply(matrix(xStarNew[iTrial,], 31, 15)[, c(3,4)], 1, sum)/denMice
      # total number nymphs is sum of elements 10-13
      totalNymphs <- apply(matrix(xStarNew[iTrial,], 31, 15)[,10:13], 1, sum)
      # nymph Bb prevalence is sum of elements 11&13 over sum of elements 10-13
      prevBbNymphs <- apply(matrix(xStarNew[iTrial,], 31, 15)[,c(11,13)], 1, sum)/totalNymphs
      # nymph Bab prevalence is sum of elements 12&13 over sum of elements 10-13
      prevBabNymphs <- apply(matrix(xStarNew[iTrial,], 31, 15)[,c(12,13)], 1, sum)/totalNymphs
      # larval burden, nymph burden are as given
      burdenLarvae <- matrix(xStarNew[iTrial,], 31, 15)[,14]
      burdenNymphs <- matrix(xStarNew[iTrial,], 31, 15)[,15]
  
      # construct single vectors with the summary state variables
      xStarSummary[iTrial,] <- matrix(c(as.vector(denMice), as.vector(burdenLarvae), as.vector(burdenNymphs), as.vector(prevBbMice),  
                                     as.vector(prevBbNymphs), as.vector(prevBabMice), as.vector(prevBabNymphs)) )
      # check if the trial was viable
      if (max(prevBbMice) > threshViable & max(prevBabMice) > threshViable) {
        xStarViable[iTrial] <- TRUE
      }
  }
  # find the max and min values across all trials for each state variable at each time point
  if (iCut == 1) {
    # for parameters from the 30% CI
    for (i in 1:217) {
      xStarMin30[i] <- ifelse(any(is.nan(xStarSummary[, i])), NaN, min(xStarSummary[xStarViable==TRUE, i], na.rm=TRUE))
      xStarMax30[i] <- ifelse(any(is.nan(xStarSummary[, i])), NaN, max(xStarSummary[xStarViable==TRUE, i], na.rm=TRUE))
    }
  }
  else {
    # for parameters from the 10% CI
    for (i in 1:217) {
      xStarMin10[i] <- ifelse(any(is.nan(xStarSummary[, i])), NaN, min(xStarSummary[xStarViable==TRUE, i], na.rm=TRUE))
      xStarMax10[i] <- ifelse(any(is.nan(xStarSummary[, i])), NaN, max(xStarSummary[xStarViable==TRUE, i], na.rm=TRUE))
    }
  }
}
      
#  trial <- matrix(iTrial, 31, 1)
day <- seq(90, 240, 5)
 
# prepare the summary state variable at each time point across all trials for the posterior median parameter values
denMiceMedian <- apply( matrix(xStarMedian, 31, 15)[,1:4], 1, sum)
prevBbMiceMedian <- apply(matrix(xStarMedian, 31, 15)[, c(2,4)], 1, sum)/denMiceMedian
prevBabMiceMedian <- apply(matrix(xStarMedian, 31, 15)[, c(3,4)], 1, sum)/denMiceMedian
totalNymphsMedian <- apply(matrix(xStarMedian, 31, 15)[,10:13], 1, sum)
prevBbNymphsMedian <- apply(matrix(xStarMedian, 31, 15)[,c(11,13)], 1, sum)/totalNymphsMedian
prevBabNymphsMedian <- apply(matrix(xStarMedian, 31, 15)[,c(12,13)], 1, sum)/totalNymphsMedian
burdenLarvaeMedian <- matrix(xStarMedian, 31, 15)[,14]
burdenNymphsMedian <- matrix(xStarMedian, 31, 15)[,15]
xStarMedian.df <- data.frame(denMiceMedian, burdenLarvaeMedian, burdenNymphsMedian, prevBbMiceMedian, prevBbNymphsMedian,  prevBabMiceMedian, prevBabNymphsMedian)

# read in observations and process
obsBI <- read.csv(file = "allDataBI.csv", header =TRUE, sep = ",")
obsBI2014 <- subset(obsBI, year == 2014); obsBI2015 <- subset(obsBI, year == 2015); obsBI2016 <- subset(obsBI, year == 2016)
# construct a single vector of observation days
obsDays <- c(obsBI2014$day, obsBI2015$day+365, obsBI2016$day+730)

#obsCT <- read.csv(file = "allDataCT.csv", header =TRUE, sep = ",")
#obsCT2014 <- subset(obsCT, year == 2014); obsCT2015 <- subset(obsCT, year == 2015); obsCT2016 <- subset(obsCT, year == 2016)
# construct a single vector of observation days
#obsDays <- c(obsCT2014$day, obsCT2015$day+365, obsCT2016$day+730)

# structure observations into single array for comparison with model in ABC algorithm, in following order:
# mouse density, larval burden, nymph burden, prevalence Bb in mice, prevalence Bb in nymphs, prevalence Bab in mice, prevalence Bab in nymphs
# note - infection prevalence in nymphs only recorded in some sessions, remaining sessions should be omitted, not set to 0, to avoid biasing outcome
obs <- c(obsBI2014$denMice, obsBI2015$denMice, obsBI2016$denMice, obsBI2014$meanBurdenLarvae, obsBI2015$meanBurdenLarvae, obsBI2016$meanBurdenLarvae, 
                  obsBI2014$meanBurdenNymphs, obsBI2015$meanBurdenNymphs, obsBI2016$meanBurdenNymphs, obsBI2014$prevBbMice, obsBI2015$prevBbMice, obsBI2016$prevBbMice,
                  obsBI2014$prevBbNymphs, obsBI2015$prevBbNymphs, obsBI2016$prevBbNymphs, obsBI2014$prevBabMice, obsBI2015$prevBabMice, obsBI2016$prevBabMice, 
                  obsBI2014$prevBabNymphs, obsBI2015$prevBabNymphs, obsBI2016$prevBabNymphs)

#obs <- c(obsCT2014$denMice, obsCT2015$denMice, obsCT2016$denMice, obsCT2014$meanBurdenLarvae, obsCT2015$meanBurdenLarvae, obsCT2016$meanBurdenLarvae, 
#         obsCT2014$meanBurdenNymphs, obsCT2015$meanBurdenNymphs, obsCT2016$meanBurdenNymphs, obsCT2014$prevBbMice, obsCT2015$prevBbMice, obsCT2016$prevBbMice,
#         obsCT2014$prevBbNymphs, obsCT2015$prevBbNymphs, obsCT2016$prevBbNymphs, obsCT2014$prevBabMice, obsCT2015$prevBabMice, obsCT2016$prevBabMice, 
#         obsCT2014$prevBabNymphs, obsCT2015$prevBabNymphs, obsCT2016$prevBabNymphs)

# for each variable, plot the envelope and data
myplots = c()

# list of variable and initial index of that variable in the simulation data and observation data
varNames <- c("mouse density", "larva burden", "nymph burden", "mouse Bb prev", "nymph Bb prev", "mouse Bab prev",  "nymph Bab prev")
dataIndex <- c(1, 31+1, 2*31+1, 3*31+1, 4*31+1, 5*31+1, 6*31+1, 7*31+1)
obsIndex <- c(1, 21+1, 2*21+1, 3*21+1, 4*21+1, 5*21+1, 6*21+1, 7*21+1)

for (i in c(1,2,3,4,6,5,7)) {
  # construct a dataframe with time series of min,max for variable i, the non-sequeitial ordering of i places the nymph prevalences at the end
  df <- data.frame(day, "xStarMin10" = xStarMin10[dataIndex[i]:(dataIndex[i+1]-1)], "xStarMax10" = xStarMax10[dataIndex[i]:(dataIndex[i+1]-1)], 
                        "xStarMin30" = xStarMin30[dataIndex[i]:(dataIndex[i+1]-1)], "xStarMax30" = xStarMax30[dataIndex[i]:(dataIndex[i+1]-1)], 
                        "xStarMedian" = unlist(xStarMedian.df[i])) 
  # construct dataframe with time series of observations for variable i
  obsDf <- data.frame(obsDays, "obs" = obs[obsIndex[i]:(obsIndex[i+1]-1)] )
  # call plotting function to construct plot
  myplots[[varNames[i]]] <- plot_envelopes(obsDf, df, varNames[i])  
}

# construct a grid of plots
x11()
#png(file="mouseTickBbBab_ABCRejEnvelope_BI.png", width=1560, height=960)
grid.arrange(grobs=myplots, nrows=2, ncols=3)
#dev.off()