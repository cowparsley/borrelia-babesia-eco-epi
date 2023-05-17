# read output from simulations produced using mouseTickBbBabABC and use easyabc to perform ABC rejection with respect to observational data

library(abc)

# clear the workspace
rm(list=ls())

# load in the data file
setwd("put your working directory here")
# ubuntu partition

# if first use - simulation output is initially in 100 files - load, combine, save in single file
# iPart <- 1
# load(file = paste('miceTicksBbBab_ABC', iPart, '_BI.RData', sep=''))
# load(file = paste('miceTicksBbBab_ABC', iPart, '_CT.RData', sep=''))

# d<- miceTicksBbBab_ABC
# for (iPart in c(2:100)) {
#   load(file = paste('miceTicksBbBab_ABC', iPart, '_BI.RData', sep=''))
#    load(file = paste('miceTicksBbBab_ABC', iPart, '_CT.RData', sep=''))
#    dPart <- miceTicksBbBab_ABC
#    d$param <- rbind(d$param, dPart$param)
#    d$stats <- rbind(d$stats, dPart$stats)
#    d$weights <- rbind(d$weights, dPart$weights)
#    d$stats_normalization <- rbind(d$stats_normalization, dPart$stats_normalization)
#    d$nsim <- d$nsim + dPart$nsim
#    d$computime <- d$computime + dPart$computime
#}

# miceTicksBbBab_ABC <- d
# save(miceTicksBbBab_ABC, file = 'miceTicksBbBab_ABC500k_BI.RData')
# save(miceTicksBbBab_ABC, file = 'miceTicksBbBab_ABC500k_CT.RData')

# if simulation output files already combined, just load the master files
load(file = 'miceTicksBbBab_ABC500k_BI.RData')
miceTicksBbBab_ABC_BI <- miceTicksBbBab_ABC

load(file = 'miceTicksBbBab_ABC500k_CT.RData')
miceTicksBbBab_ABC_CT <- miceTicksBbBab_ABC

# --- process the observation data --------------------
# read in the empirical data, extract observations, separate into years
obsBI <- read.csv(file = "allDataBI.csv", header =TRUE, sep = ",")
obsBI2014 <- subset(obsBI, year == 2014); obsBI2015 <- subset(obsBI, year == 2015); obsBI2016 <- subset(obsBI, year == 2016)

obsCT <- read.csv(file = "allDataCT.csv", header =TRUE, sep = ",")
obsCT2014 <- subset(obsCT, year == 2014); obsCT2015 <- subset(obsCT, year == 2015); obsCT2016 <- subset(obsCT, year == 2016)

# construct a single vector of observation days
daysBI <- c(obsBI2014$day, obsBI2015$day, obsBI2016$day)
daysCT <- c(obsCT2014$day, obsCT2015$day, obsCT2016$day)

# structure observations into single array for comparison with model in ABC algorithm, in following order:
# mouse density, larval burden, nymph burden, prevalence Bb in mice, prevalence Bb in nymphs, prevalence Bab in mice, prevalence Bab in nymphs
# note - infection prevalence in nymphs only recorded in some sessions, remaining sessions should be omitted, not set to 0, to avoid biasing outcome
observationsBI <- c(obsBI2014$denMice, obsBI2015$denMice, obsBI2016$denMice, obsBI2014$meanBurdenLarvae, obsBI2015$meanBurdenLarvae, obsBI2016$meanBurdenLarvae, 
                  obsBI2014$meanBurdenNymphs, obsBI2015$meanBurdenNymphs, obsBI2016$meanBurdenNymphs, obsBI2014$prevBbMice, obsBI2015$prevBbMice, obsBI2016$prevBbMice,
                  obsBI2014$prevBbNymphs, obsBI2015$prevBbNymphs, obsBI2016$prevBbNymphs, obsBI2014$prevBabMice, obsBI2015$prevBabMice, obsBI2016$prevBabMice, 
                  obsBI2014$prevBabNymphs, obsBI2015$prevBabNymphs, obsBI2016$prevBabNymphs)

observationsCT <- c(obsCT2014$denMice, obsCT2015$denMice, obsCT2016$denMice, obsCT2014$meanBurdenLarvae, obsCT2015$meanBurdenLarvae, obsCT2016$meanBurdenLarvae, 
                  obsCT2014$meanBurdenNymphs, obsCT2015$meanBurdenNymphs, obsCT2016$meanBurdenNymphs, obsCT2014$prevBbMice, obsCT2015$prevBbMice, obsCT2016$prevBbMice,
                  obsCT2014$prevBbNymphs, obsCT2015$prevBbNymphs, obsCT2016$prevBbNymphs, obsCT2014$prevBabMice, obsCT2015$prevBabMice, obsCT2016$prevBabMice, 
                  obsCT2014$prevBabNymphs, obsCT2015$prevBabNymphs, obsCT2016$prevBabNymphs)

# remove na's
obsBI <- observationsBI[!is.na(observationsBI)]
obsCT <- observationsCT[!is.na(observationsCT)]
  
# --- process the simulation data -----------------------------
# set the number of trials we wish to extract from the data file
nTrials <- 500000
  
# extract the state variables for each trial, each variable at 3 x 7 = 21 time points except nymphs where fewer time points
# are used in some years - see mouseTickBbBabABC.R for details

# 1:21 - denMice, 22:42 - larval burden, 43:63 - nymph burden, 64:84 - prevBbMice, 85:100 - prevBbNymphs, 
# 101:121 - prevBabMice, 122:147 - prevBabNymphs

# record the parameter names
paramNames <- c("r", "mu", "K", "omegaM", "tauE", "tauL", "tauN", "etaE", "etaL", "etaN","lambda", "Omega", "omegaL", "D", "delta", "nu", "betaML1", "betaNM1",
                "gamma1", "betaML2", "betaNM2", "sigma", "xi", "alpha")

for (iSite in 1:2) {
   # select the site and use observations and abc output for that site
   if (iSite == 1) { miceTicksBbBab_ABC <- miceTicksBbBab_ABC_BI; obs <- obsBI }
   else { miceTicksBbBab_ABC <- miceTicksBbBab_ABC_CT; obs <- obsCT }

   # remove all trials that did not lead to positive Bb and Bab prevalence in mice
   maxPrevBbMice <-  apply(miceTicksBbBab_ABC$stats[, 64:84], 1, max)       # get the max Bb prevalence in each trial
   maxPrevBabMice <- apply(miceTicksBbBab_ABC$stats[, 101:121], 1, max)     # same for Bab
   threshViable <- 1e-2                                                     # minimum max prevalence for which trial considered viable i.e. R0 > 1
   iViable <- maxPrevBbMice > threshViable &  maxPrevBabMice > threshViable # indices of viable trials - where Bb and Bab max prevalence higher than threshold

   # extract the parameters of the viable trials only
   parameters.df <- data.frame(miceTicksBbBab_ABC$param[iViable,1:24])     
   colnames(parameters.df) <- paramNames;
  
   # extract the summary statistics (ie state variables) for the viable trials only
   sumstats.df <- data.frame(miceTicksBbBab_ABC$stats[iViable,])
  
   # a random sample of parameters from all the trials will generate the prior parameter distributions
   mysample <- sample(nTrials, 100000)    # indices  
   parameterSample.df <- data.frame(miceTicksBbBab_ABC$param[mysample,1:24])   # parameter sets
   colnames(parameterSample.df) <- paramNames;

   # run abc rejection sampler to select trials closest to observed data
   # target is the observed data - records of populations sizes and prevalences at 7 time points each year
   # sumstats.df are the corresponding model state variables generated by simulation with the parameter values listed in parameters.df
   tol = 0.005  # proportion of trials to keep
   abc.var <-  abc(target=obs, param=parameters.df, sumstat=sumstats.df, tol=tol, method ="rejection")

   # save the result in an informatively named variable
   if (iSite == 1) { 
      assign(paste('rej', gsub("\\.", "", tol), '_BI', sep=''), abc.var) # output of abc rejection algorithm
      assign('parameterSample.df_BI', parameterSample.df )   # prior parameters
      assign('parameters.df_BI', parameters.df )             # parameters of viable tgrials only
      
      }
   else { 
      assign(paste('rej', gsub("\\.", "", tol), '_CT', sep=''), abc.var) }
      assign('parameterSample.df_CT', parameterSample.df )
      assign('parameters.df_CT', parameters.df )
   }

# save all the output for plotting etc elsewhere
save(parameters.df_BI, parameters.df_CT, parameterSample.df_BI, parameterSample.df_CT, obsBI, obsCT, rej0005_BI, rej0005_CT, file="miceTicksBbBab_ABCRej0005_BICT.RData")

