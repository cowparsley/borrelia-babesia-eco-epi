# borrelia-babesia-eco-epi

Codes for the manuscript Ecological interactions driving population dynamics of two tick-borne pathogens, Borrelia burgdorferi and Babesia microti

miceRMark.R

This script processes the field data and then uses the Mark program, via the RMark interface, to estimate survival, observation and state transition probabilities for a multistate Markov model.

miceMSM.R

This script processes the field data and then uses the msm package, to estimate continuous-time state transition rates for a multistate Markov model. These are then used to calculate discrete-time state transition probabilities.

mouseTickBbBabSimluation.R

This function solves the semidiscrete model (specificed in mouseTickBbBabODE.R) for a given parameter set theta. The ode model specified in mouseTickBbBabODE is used in the active season, discrete time update is used for inactive season. The initial condition is specified in the function. Starting from that point, the system is solved for 50 years. At this point it is close to steady state. It is then solved for a further 3 years to generate 3 years of output for comparison with data.

mouseTickBbBabODE.R

This function specifies the ode system that describes the dynamics during the active season.

mouseTickBbBabABC.R

This code calls the easyabc package to run simulations of the mouse-tick-borrelia-babesia model with parameters drawn from specified prior distributions. Simulation results are stored for ABC rejection post-processing. No rejection sampling is done here.

mouseTickBbBabABCRej.R

This code reads output from simulations produced using mouseTickBbBabABC and uses the abc package to perform ABC rejection with respect to the observational data

mouseTickBbBabABCRej_Plot.R

This code reads in output of mouseTickBBBabABCRej and process it to plot prior and posterior parameter distributions

mouseTickBbBabABCRejEnvelopes.R

This code generates simulation output for plotting best estimate trajectories and envelopes. It generates parameter values from the centre of the posterior distributions and uses them to run parallel simulations of mouse-tick-borrelia-babesia model

mouseTickBbBabABCRejEnvelopesPlot.R

This code load simulation output for median parameter values and parameter envelopes, and plots summary statistic envelopes together with observations
