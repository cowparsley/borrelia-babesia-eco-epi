# code to read in output of mouseTickBBBabABCRej and process it to plot prior and posterior parameter distributions

library(abc)
library(gridExtra)
library(ggplot2)
library(knitr)

# clear the workspace
rm(list=ls())

# specify a function to plot densities of the prior and rejection posterior parameter distributions 
plot_densities = function (parameterSample.df.BI, parameterSample.df.CT, rej0005.BI, rej0005.CT, paramName) {
  # for parameter paramName, plot priors for BI and CT, and rejection posteriors for BI and CT
  ggplot(data = as.data.frame(rej0005.BI$unadj.values), aes_string(x=paramName)) + geom_density(color='skyblue2', size = 1) + 
    geom_density(data = as.data.frame(rej0005.CT$unadj.values), aes_string(x=paramName), color='salmon', size = 1) +
    # plot prior
    geom_density(data = parameterSample.df.BI, aes_string(x=paramName), color='skyblue2', linetype='dotted', size = 1) + 
    geom_density(data = parameterSample.df.CT, aes_string(x=paramName), color='salmon', linetype='dotted', size = 1) +
    if (paramName == 'lambda') {xlim(0, 0.0015)} else if (paramName == 'Omega') {xlim(0, 1e5)} else if (paramName == 'omegaL') {xlim(0, 1)}
    else if (paramName == 'D') {xlim(0, 80)} else if (paramName == 'K') {xlim(0, 100)} else if (paramName == 'betaML2') {xlim(0, 0.7)} 
    else if (paramName == 'sigma') {xlim(0, 3)} else if (paramName == 'xi') {xlim(0, 1.6)} 
}

# load in the data file
setwd("put your working directory here")

# load file with observations and trials selected by abc rejection, generated using mouseTickBbBabABCRej.R 
load(file="miceTicksBbBab_ABCRej0005_BICT.RData")

# for each parameter generate density plots, store in a list
myplots = c()
for (paramName in names(parameters.df_BI)) {
  myplots[[paramName]] <- plot_densities(parameterSample.df_BI, parameterSample.df_CT, rej0005_BI, rej0005_CT, paramName)  
}

# construct a grid of plots
# png(file="miceTicksBbBab_ABCRej0005_BI.png")  # use to plot direct to png instead of window
x11() # use to plot to new window
nGrid = ceiling(sqrt(length(parameters.df_BI)))
grid.arrange(grobs=myplots, nrows=nGrid, ncols=nGrid)
#dev.off()  # use to stop plotting to png

# calculate and print summary statistics of parameter distributions
quants.BI <- apply(rej0005_BI$unadj.values, 2, quantile, probs = c(0.1, 0.5, 0.9))
quants.CT <- apply(rej0005_CT$unadj.values, 2, quantile, probs = c(0.1, 0.5, 0.9))

print(kable(quants.BI[1,]))
print(kable(quants.CT[1,]))


