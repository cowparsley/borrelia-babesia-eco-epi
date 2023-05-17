# run parallel simulations of mouse-tick-borrelia-babesia model with samples from best estimate parameter distributions
# to get detailed time series for constructing envelope

library(foreach)
library(doParallel)
library(deSolve)                       # ode solvers
library(edfun)

# clear the memory
rm(list=ls())

# set up the parralel cores
registerDoParallel(cores=12)

# set working directory
setwd("put your working directory here")

# specify function that simulates model and returns a single vector composed of each variable at each time point
# will be run in parallel
simulateData = function(theta) {
  
  # function that will solve ode for given parameter set theta, to be called from foreach loop
  library(deSolve)                       # ode solvers
  source("./mouseTickBbBabODE.R")        # specifies the model for the ode solver 
  
  tauM <- 90                                       # start time, day of the year
  tauF <- 240                                      # end time, day of the year
  tauW <- 365 - tauF + tauM                        # duration of winter period, days
  day <- seq(from = tauM+5, to = tauF-5, by = 5)   # output daily
  tSpan <- c(tauM, day, tauF)  
  
  # initial conditions, based on values from asymptotic solution with parameter values that look about right for BI2016
  # initM0 = 15.48; initM1 = 3.87e-5; initM2 = 5.57; initM12 = 1.65e-5; initE = 8993.4; 
  # initN0 = 16002.0; initN1 = 550.6; initN2 = 3197.5; initN12 = 643.9; 
  
  # initial conditions for CT
  initM0 = 9.372; initM1 = 8.269e-4; initM2 = 1.864e-4; initM12 = 1.648e-8; initE = 2.652e4; 
  initN0 = 1.651e4; initN1 = 9.222e2; initN2 = 1.8525e3; initN12 = 3.681e2; 
  
  # 50 years spin-up, output last year
  for (year in c(1:50)) {
    # input the initial conditions    
    initConditions <- c(M0=initM0, M1=initM1, M2=initM2, M12=initM12, L=0, W0=0, W1=0, W2=0, W12=0, N0=0, N1=0, N2=0, N12=0, LB=0, NB=0)
    
    # prepare parameters to be passed to ode solver
    params <- c(  r = theta[2], mu = theta[3], K = theta[4], tauE = theta[6], tauL = theta[7], tauN = theta[8], 
                  etaE = theta[9], etaL = theta[10], etaN = theta[11], lambda = theta[12], Omega = theta[13], 
                  D = theta[15], delta = theta[16], nu = theta[17], betaML1 = theta[18], betaNM1 = theta[19],  gamma1 = theta[20],
                  betaML2 = theta[21], betaNM2 = theta[22], sigma = theta[23], xi = theta[24], 
                  alpha = theta[25], initE=initE, initN0=initN0, initN1=initN1, initN2=initN2, initN12=initN12 ) 
    
    # call the ode solver
    xStar <- ode(initConditions, tSpan, mouseTickBbBabODE, params, method = "lsode")
    
    # ready for next year, construct inital conditions and initial values for implicit variables to account for overwinter mortality 
    # and, in mice, recovery
    initM0 =  theta[5]*( xStar[[dim(xStar)[1],"M0"]] + xStar[[dim(xStar)[1],"M1"]]*(1-exp(-theta[20]*tauW)) ) 
    initM1 =  theta[5]*( xStar[[dim(xStar)[1],"M1"]]*exp(-theta[20]*tauW) )
    initM2 =  theta[5]*( xStar[[dim(xStar)[1],"M2"]] +  xStar[[dim(xStar)[1],"M12"]]*(1-exp(-theta[20]*tauW)) )
    initM12 = theta[5]*( xStar[[dim(xStar)[1],"M12"]]*exp(-theta[20]*tauW) )
    initE =   theta[14]*xStar[[dim(xStar)[1],"L"]]    
    initN0 =  theta[14]*xStar[[dim(xStar)[1],"W0"]]   
    initN1 =  theta[14]*xStar[[dim(xStar)[1],"W1"]]   
    initN2 =  theta[14]*xStar[[dim(xStar)[1],"W2"]]   
    initN12 = theta[14]*xStar[[dim(xStar)[1],"W12"]]   
  }
  
  # return single vector composed of each variable at each time points
  return( t(as.vector(xStar[,2:dim(xStar)[2]])) )
}

#-------------- main body -----------------------

# load the best fit parameter values from previous analysis with tolerance 0.005
load(file = "miceTicksBbBab_ABCRej0005_BICT.RData")

# extract the parameters of the best fits into a dataframe
parameters.df <- data.frame(rej0005_BI$unadj.values)
# parameters.df <- data.frame(rej0005_CT$unadj.values)

paramNames <- c("r", "mu", "K", "omegaM", "tauE", "tauL", "tauN", "etaE", "etaL", "etaN","lambda", "Omega", "omegaL", "D", "delta", "nu", "betaML1", "betaNM1",
                "gamma1", "betaML2", "betaNM2", "sigma", "xi", "alpha")
colnames(parameters.df) <- paramNames

# for each parameter find the median value of the posterior
medianParams <- apply(rej0005_BI$unadj.values, 2, median, na.rm = TRUE)
# medianParams <- apply(rej0005_CT$unadj.values, 2, median, na.rm = TRUE)

# for each parameter, separate out values that fall between the 25% and 75% quantiles then generate 1000 entirely new parameter 
# sets with all values drawn from this range of the posterior distributions

# set up a dataframe for the new parameter sets, same column names as dataframe for original parameters
new_params.df <- parameters.df[nrow(parameters.df)+1:1000,]

# for each parameter, get the points in the posterior that fall between 35% and 65% quantiles using cut,
# generate a distribution from these points using edfun, generate 1000 random samples from this distribution using edfun,
# store the result in the dataframe

nCl <- ncol(parameters.df)                 # number of parameters to be processed
for (iCol in 1:nCl) {
  params <- parameters.df[, iCol]          # extract all values of the given parameter
  params.qcuts <- cut(params, breaks = quantile(params, probs = c(0, 0.25, 0.75, 1)), include.lowest = TRUE, labels = c("q35l", "q30", "q65u")) # find the 25 and 75% breakpoints and label values according to which interval they fall in
  params.cut <- params[params.qcuts=="q30"]  # get all the parameters in the middle interval, i.e. exclusing the top and bottom 25%
  params_dist.cut <- edfun(params.cut)       # estimate a distribution for these parameters
  new_params <- params_dist.cut$rfun(1000)   # generate a 1000 samples from this distribution, and use these for simulations
  new_params.df[,iCol] <- new_params
}

# run simulations in parallel for each parameter set generated from the central 50% posterior
out <- foreach (iParam = 1:nrow(new_params.df)) %dopar% {
     # construct a parameter vector theta that mimics one used by easyabc, first element is NaN - not used
     simulateData( c(NaN, as.numeric(new_params.df[iParam,]) ) )
}
# get results into a matrix  
xStarNew <-  matrix(unlist(out), ncol=465, byrow=TRUE)

# simulate using median values of parameter distributions, and best fit parameter values
xStarMedian <- simulateData( as.numeric( c(NaN, medianParams ) ) )

save(xStarNew, xStarMedian, file= 'miceTicksBbBab_ABCRej0005_Envelopes30_BI.RData')
#save(xStarNew, xStarMedian, file= 'miceTicksBbBab_ABCRej0005_Envelopes30_CT.RData')


