# use easyabc package to runs simulations of mouse-tick-borrelia-babesia model with parameters drawn from prior distributions and
# store results for ABC rejection post-processing. No rejection sampling is done here. 

library(EasyABC)

# clear the memory
rm(list=ls())

# change the working directory
setwd("put your working directory here")

# load the simulation function that will be called by the ABC function
source('./mouseTickBbBabSimulation.R')

# set up uniform priors for each parameter - BI and CT
priorParamsUnif.df <- data.frame("r" = c(0.07, 0.2), "mu" = c(0,0.05), "K" = c(10, 70), "omegaM" = c(0,1), "tauE" = c(90, 140), 
                                 "tauL" = c(170, 220), "tauN" = c(90, 135), "etaE" = c(0.05, 0.2), "etaL" = c(0.05, 0.2), "etaN" = c(0.05, 0.2),
                                 "lambda" = c(1e-4, 1e-3), "Omega" = c(1e4, 1e5), "omegaL"= c(0.2, 0.8), "D" = c(0, 100), "delta" = c(0.2, 0.4),
                                  "nu" = c(0, 1), "betaML1" = c(0, 1), "betaNM1" = c(0,1), "gamma1" = c(0, 0.05), "betaML2" = c(0,1), "betaNM2"= c(0,1),
                                  "sigma" = c(0,3), "xi" = c(0, 1), "alpha" = c(1,2), "dummy0" = c(0, 0.0001), "dummy1" = c(0.9999, 1) )
 
# set up lognormal priors for each parameter based on mouse-tick only parameter estimation - BI
priorParamsLNorm.df <- data.frame("r" = c(-2.18, 0.473), "mu" = c(-4.044,1.014), "K" = c(3.790, 0.216), "omegaM" = c(-0.903, 0.870), "tauE" = c(4.751, 0.124), 
                                  "tauL" = c(5.298, 0.049), "tauN" = c(4.693, 0.112), "etaE" = c(-2.160, 0.385), "etaL" = c(-2.201, 0.393), "etaN" = c(-2.144, 0.379),
                                  "lambda" = c(-7.860, 0.563), "Omega" = c(10.428, 0.528), "omegaL"= c(-0.904, 0.383), "D" = c(3.760, 0.833), "delta" = c(-1.188, 0.190),
                                  "nu" = c(log(0.74),0.2), "betaML1" = c(log(0.83),0.2), "betaNM1" = c(log(0.83),0.2), "betaML2" = c(log(0.37),0.2), 
                                  "betaNM2"= c(log(0.83),0.2), "gamma1" = c(log(0.05),0.2), "sigma" = c(log(1.54),0.2), 
                                  "xi" = c(log(0.87),0.2), "alpha" = c(log(1),0.2) )

# set up lognormal priors for each parameter based on mouse-tick only parameter estimation - CT
#priorParamsLNorm.df <- data.frame("r" = c(-2.30, 0.538), "mu" = c(-3.973,0.974), "K" = c(3.401, 0.376), "omegaM" = c(-0.941, 0.909), "tauE" = c(4.733, 0.124), 
#                                  "tauL" = c(5.311, 0.059), "tauN" = c(4.686, 0.114), "etaE" = c(-2.151, 0.383), "etaL" = c(-2.219, 0.399), "etaN" = c(-2.125, 0.376),
#                                  "lambda" = c(-7.847, 0.638), "Omega" = c(9.970, 0.495), "omegaL"= c(-0.960, 0.382), "D" = c(4.000, 0.670), "delta" = c(-1.164, 0.184),
#                                  "nu" = c(log(0.74),0.2), "betaML1" = c(log(0.83),0.2), "betaNM1" = c(log(0.83),0.2), "betaML2" = c(log(0.37),0.2), 
#                                  "betaNM2"= c(log(0.83),0.2), "gamma1" = c(log(0.05),0.2), "sigma" = c(log(1.54),0.2), 
#                                  "xi" = c(log(0.87),0.2), "alpha" = c(log(1),0.2) )

# set up the actual priors to be used, a mix of uniform and lognormal depending on how good the lognormal fits were in the mouse-tick only estimation stage 
# the dummy0 and dummy1 priors are used in the prior comparison component of easyabc to ensure certain other parameters are in (0,1)
# BI
priorsMixed_miceTicksBbBab = list( c("unif",priorParamsUnif.df[,"r"]), c("unif",priorParamsUnif.df[,"mu"]), c("lognormal",priorParamsLNorm.df[,"K"]),
                                   c("unif",priorParamsUnif.df[,"omegaM"]), c("unif",priorParamsUnif.df[,"tauE"]), c("lognormal",priorParamsLNorm.df[,"tauL"]), 
                                   c("lognormal",priorParamsLNorm.df[,"tauN"]), c("unif",priorParamsUnif.df[,"etaE"]), c("unif",priorParamsUnif.df[,"etaL"]), 
                                   c("unif",priorParamsUnif.df[,"etaN"]), c("lognormal",priorParamsLNorm.df[,"lambda"]), c("lognormal",priorParamsLNorm.df[,"Omega"]), 
                                   c("lognormal",priorParamsLNorm.df[,"omegaL"]), c("unif",priorParamsUnif.df[,"D"]), c("unif",priorParamsUnif.df[,"delta"]),
                                   c("lognormal",priorParamsLNorm.df[,"nu"]), c("lognormal",priorParamsLNorm.df[,"betaML1"]), c("lognormal",priorParamsLNorm.df[,"betaNM1"]), 
                                   c("unif",priorParamsUnif.df[,"gamma1"]), c("lognormal",priorParamsLNorm.df[,"betaML2"]), c("lognormal",priorParamsLNorm.df[,"betaNM2"]), 
                                   c("lognormal",priorParamsLNorm.df[,"sigma"]), c("lognormal",priorParamsLNorm.df[,"xi"]), 
                                   c("unif",priorParamsUnif.df[,"alpha"]), c("unif",priorParamsUnif.df[,"dummy0"]), c("unif",priorParamsUnif.df[,"dummy1"]) )

# CT
# priorsMixed_miceTicksBbBab = list( c("unif",priorParamsUnif.df[,"r"]), c("unif",priorParamsUnif.df[,"mu"]), c("lognormal",priorParamsLNorm.df[,"K"]),
#                                   c("unif",priorParamsUnif.df[,"omegaM"]), c("unit",priorParamsUnif.df[,"tauE"]), c("lognormal",priorParamsLNorm.df[,"tauL"]), 
#                                   c("lognormal",priorParamsLNorm.df[,"tauN"]), c("unif",priorParamsUnif.df[,"etaE"]), c("unif",priorParamsUnif.df[,"etaL"]), 
#                                    c("unif",priorParamsUnif.df[,"etaN"]), c("lognormal",priorParamsLNorm.df[,"lambda"]), c("lognormal",priorParamsLNorm.df[,"Omega"]), 
#                                   c("lognormal",priorParamsLNorm.df[,"omegaL"]), c("unif",priorParamsUnif.df[,"D"]), c("unif",priorParamsUnif.df[,"delta"]),
#                                   c("lognormal",priorParamsLNorm.df[,"nu"]), c("lognormal",priorParamsLNorm.df[,"betaML1"]), c("lognormal",priorParamsLNorm.df[,"betaNM1"]), 
#                                   c("unif",priorParamsUnif.df[,"gamma1"]), c("lognormal",priorParamsLNorm.df[,"betaML2"]), c("lognormal",priorParamsLNorm.df[,"betaNM2"]), 
#                                   c("lognormal",priorParamsLNorm.df[,"sigma"]), c("lognormal",priorParamsLNorm.df[,"xi"]), 
#                                   c("unif",priorParamsUnif.df[,"alpha"]), c("unif",priorParamsUnif.df[,"dummy0"]), c("unif",priorParamsUnif.df[,"dummy1"]) )

# run rejection abc_algorithm without any observation data to generate simulation output only
# in case of crashes, run in 100 parts, each of 5000 trials, and save output after each part
for (iPart in c(1:100)) {
  miceTicksBbBab_ABC <-ABC_rejection(model=mouseTickBbBabSimulationCompact, prior = priorsMixed_miceTicksBbBab, 
                                   prior_test="x01>x02 & x07<x06 & x05<x06 & x16>x25 & x16<x26 & x17>x25 & x17<x26 & x18>x25 & x18<x26 & x20>x25 & x20<x26 & x21>x25 & x21<x26", 
                                   nb_simul = 5000, summary_stat_target=NULL, use_seed = TRUE, n_cluster = 8, verbose = FALSE)
  # construct the filename for output by assigning number iPart to root, and save
  assign(paste('miceTicksBbBab_ABC', iPart, sep=''), miceTicksBbBab_ABC); 
  save(miceTicksBbBab_ABC, file = paste('miceTicksBbBab_ABC', iPart, '_BI.RData', sep=''))
# save(miceTicksBbBab_ABC, file = paste('miceTicksBbBab_ABC', iPart, '_CT.RData', sep=''))
}


