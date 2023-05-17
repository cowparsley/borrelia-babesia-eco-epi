# Function that computes solution of semi-discrete time model for mouse and tick population dynamics with
# borrelia and babesia epidemiology using parameters specified in theta; first element of theta is the 
# seed value for the pseudo random number generator - this is not used here); all parameters are the same each year

# author: Ben Adams, University of Bath
# updated: Nov 2022

mouseTickBbBabSimulationCompact <- function(theta) {
  # this function will be called on multiple cores, so need to load libraries and helper 
  # functions on each core
  library(deSolve)                       # ode solvers

  source("./mouseTickBbBabODE.R")        # specifies the model for the ode solver 
  
  seedVal <- theta[1]                    # not used here

  # write to file to monitor progress within multiple cores
  conn <- file( sprintf("./output_%d.txt" , Sys.getpid()) , open = "a" )
  write.table(t(theta), conn , append = TRUE , row.names = FALSE, col.names = FALSE )
  close( conn )
  
  # model is solved in continuous time between tauM and tauF each year, with discrete time update to 
  # state variables between tauF in year t and tauM in year t+1
  
  tauM <- 90                                       # start time, day of the year
  tauF <- 240                                      # end time, day of the year
  tauW <- tauM + (365 - tauF)                      # duration of winter period
  
  # output times in years 2014, 2015, 2016 and initial conditions appropriate for comparison with BI data
  # day2014 <- c(142, 155, 168, 182, 197, 211, 225)  
  # day2015 <- c(147, 160, 174, 189, 203, 217, 231) 
  # day2016 <- c(145, 159, 174, 188, 202, 216, 230) 
  # tSpan <- c(tauM, day2014, tauF)                   # timespan of 1 year for use during spin up, any year will do
  
  # initial conditions, based on values from asymptotic solution with parameter values that produce trajectories that look about right for BI2016
  # initM0 = 15.48; initM1 = 3.87e-5; initM2 = 5.57; initM12 = 1.65e-5; initE = 8993.4; 
  # initN0 = 16002.0; initN1 = 550.6; initN2 = 3197.5; initN12 = 643.9; 
  
  # output times in years 2014, 2015, 2016 and initial conditions appropriate for comparison with CT data
  day2014 <- c(133, 148, 161, 175, 189, 203, 217)   
  day2015 <- c(139, 153, 167, 181, 196, 209, 223) 
  day2016 <- c(145, 159, 173, 187, 201, 215, 229) 
  
  tSpan <- c(tauM, day2014, tauF)                   # timespan of 1 year for use during spin up, any year will do
  
  # initial conditions for CT
  initM0 = 9.372; initM1 = 8.269e-4; initM2 = 1.864e-4; initM12 = 1.648e-8; initE = 2.652e4; 
  initN0 = 1.651e4; initN1 = 9.222e2; initN2 = 1.8525e3; initN12 = 3.681e2; 
  
  #  day <- seq(from = tauM+1, to = tauF-1, by = 1)     #  use if daily output required
  
  # 50 years spin-up, hopefully to steady state
  for (year in c(1:50)) {
    # input initial conditions    
    initConditions <- c(M0=initM0, M1=initM1, M2=initM2, M12=initM12, L=0, W0=0, W1=0, W2=0, W12=0, N0=0, N1=0, N2=0, N12=0, LB=0, NB=0)
    
    # prepare parameters to be passed to ode solver i.e associate names to the parameters in the vector theta, and add some initial population sizes
    params <- c(  r = theta[2], mu = theta[3], K = theta[4], tauE = theta[6], tauL = theta[7], tauN = theta[8], 
                  etaE = theta[9], etaL = theta[10], etaN = theta[11], lambda = theta[12], Omega = theta[13], 
                  D = theta[15], delta = theta[16], nu = theta[17], betaML1 = theta[18], betaNM1 = theta[19], gamma1 = theta[20],
                  betaML2 = theta[21], betaNM2 = theta[22], sigma = theta[23], xi = theta[24], 
                  alpha = theta[25], initE=initE, initN0=initN0, initN1=initN1, initN2=initN2, initN12=initN12 ) 
    
    # call the ode solver
    xStar <- ode(initConditions, tSpan, mouseTickBbBabODE, params, method = "lsoda")
    
    # ready for the next year, construct initial conditions and initial values for implicit variables to account for overwinter mortality 
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
  
  # another 3 years, for output
  # input the initial conditions, for overwintered populations use final year of 50 year spin up above   
  for (year in c(2014:2016)) {
    initConditions <- c(M0=initM0, M1=initM1, M2=initM2, M12=initM12, L=0, W0=0, W1=0, W2=0, W12=0, N0=0, N1=0, N2=0, N12=0, LB=0, NB=0)
  
    params <- c(  r = theta[2], mu = theta[3], K = theta[4], tauE = theta[6], tauL = theta[7], tauN = theta[8], 
                  etaE = theta[9], etaL = theta[10], etaN = theta[11], lambda = theta[12], Omega = theta[13], 
                  D = theta[15], delta = theta[16], nu = theta[17], betaML1 = theta[18], betaNM1 = theta[19], gamma1 = theta[20],  
                  betaML2 = theta[21], betaNM2 = theta[22], sigma = theta[23], xi = theta[24], 
                  alpha = theta[25], initE=initE, initN0=initN0, initN1=initN1, initN2=initN2, initN12=initN12 ) 
  
    # set the output times for the solver so they match up with data
    if (year == 2014) {
      tSpan <- c(tauM, day2014, tauF);  
      xStar <- ode(initConditions, tSpan, mouseTickBbBabODE, params, method = "lsoda")
      xStar2014 <- xStar      
      } 
    else if (year == 2015) {
      tSpan <- c(tauM, day2015, tauF);  
      xStar <- ode(initConditions, tSpan, mouseTickBbBabODE, params, method = "lsoda")
      xStar2015 <- xStar      
      }
    else if (year == 2016) {
      tSpan <- c(tauM, day2016, tauF);  
      xStar <- ode(initConditions, tSpan, mouseTickBbBabODE, params, method = "lsoda")
      xStar2016 <- xStar
      }
    
    # ready for next year, construct initial conditions and initial values for implicit variables to account for overwinter mortality 
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

  # the model output supplied to the ABC algorithm needs to correspond with sample data; so it must consist only of population sizes at interior time points 
  # i.e. we need to remove the time array and the outputs at tauM and tauF of each year produced by the ode solver. The solver also resets the time to 0 at 
  # the beginning of each year, but the data records time as total days since Jan 1 2014, so need to correct that here.  
  
  # then ABC algorithm requires state variables as single array
  denMice2014 <- xStar2014[2:8,"M0"]+xStar2014[2:8,"M1"]+xStar2014[2:8,"M2"]+xStar2014[2:8,"M12"]  # total mouse population at times 1 - 7 in year 2014
  denMice2014[denMice2014 < 1e-10] = 0                                                             # set very small populations to 0 
  denMice2015 <- xStar2015[2:8,"M0"]+xStar2015[2:8,"M1"]+xStar2015[2:8,"M2"]+xStar2015[2:8,"M12"]  # same for 2015
  denMice2015[denMice2015 < 1e-10] = 0
  denMice2016 <- xStar2016[2:8,"M0"]+xStar2016[2:8,"M1"]+xStar2016[2:8,"M2"]+xStar2016[2:8,"M12"]  # same for 2016
  denMice2016[denMice2016 < 1e-10] = 0
  denMice <- c(denMice2014, denMice2015, denMice2016)                                              # combine all three years into single array

  prevBbMice2014 <- (xStar2014[2:8,"M1"]+xStar2014[2:8,"M12"])/denMice2014    # total prevalence of Bb and (next line) Bab in mice in year 2014      
  prevBabMice2014 <- (xStar2014[2:8,"M2"]+xStar2014[2:8,"M12"])/denMice2014  
  prevBbMice2014[denMice2014 == 0] = 0                                        # if mouse population size 0, set prevalence to 0
  prevBabMice2014[denMice2014 == 0] = 0
  
  prevBbMice2015 <- (xStar2015[2:8,"M1"]+xStar2015[2:8,"M12"])/denMice2015    # same for 2015
  prevBabMice2015 <- (xStar2015[2:8,"M2"]+xStar2015[2:8,"M12"])/denMice2015
  prevBbMice2015[denMice2015 == 0] = 0
  prevBabMice2015[denMice2015 == 0] = 0
  
  prevBbMice2016 <- (xStar2016[2:8,"M1"]+xStar2016[2:8,"M12"])/denMice2016    # same for 2016
  prevBabMice2016 <- (xStar2016[2:8,"M2"]+xStar2016[2:8,"M12"])/denMice2016
  prevBbMice2016[denMice2016 == 0] = 0
  prevBabMice2016[denMice2016 == 0] = 0
  
  prevBbMice <- pmax(c(prevBbMice2014, prevBbMice2015, prevBbMice2016),0)       # combine all three years into single array
  prevBabMice <- pmax(c(prevBabMice2014, prevBabMice2015, prevBabMice2016),0)

  # BI data only records nymph infection prevalence at first 5 time points in 2014 and first 4 in 2016 - need to account for this
  # CT data only records nymph infection prevalence at first 2 time points in 2014 and first 5 in 2016 - need to account for this
  
#  i2014 = 6 # last nymph index for 2014 (BI)
#  i2016 = 5 # last nymph index for 2016 (BI)
  i2014 = 3 # last nymph index for 2014 (CT)
  i2016 = 6 # last nymph index for 2016 (CT)
  
  
  totalNymphs2014 <- (xStar2014[2:i2014,"N0"]+xStar2014[2:i2014,"N1"]+xStar2014[2:i2014,"N2"]+xStar2014[2:i2014,"N12"])  # total nymph population 2014
  totalNymphs2014[totalNymphs2014 < 1e-10] = 0                                                                           # set anything very small to 0
  totalNymphs2015 <- (xStar2015[2:8,"N0"]+xStar2015[2:8,"N1"]+xStar2015[2:8,"N2"]+xStar2015[2:8,"N12"])                  # same for 2015
  totalNymphs2015[totalNymphs2015 < 1e-10] = 0
  totalNymphs2016 <- (xStar2016[2:5,"N0"]+xStar2016[2:i2016,"N1"]+xStar2016[2:i2016,"N2"]+xStar2016[2:i2016,"N12"])      # same for 2016
  totalNymphs2016[totalNymphs2016 < 1e-10] = 0
  totalNymphs <- c(totalNymphs2014, totalNymphs2015, totalNymphs2016)                                                    # all years combined into one array
  
  prevBbNymphs2014 <- (xStar2014[2:i2014,"N1"]+xStar2014[2:i2014,"N12"])/totalNymphs2014    # prevalence of Bb and (next line) Bab in nymphs, 2014
  prevBabNymphs2014 <- (xStar2014[2:i2014,"N2"]+xStar2014[2:i2014,"N12"])/totalNymphs2014
  prevBbNymphs2014[totalNymphs2014 == 0] = 0
  prevBabNymphs2014[totalNymphs2014 == 0] = 0
  
  prevBbNymphs2015 <- (xStar2015[2:8,"N1"]+xStar2015[2:8,"N12"])/totalNymphs2015            # same for 2015
  prevBabNymphs2015 <- (xStar2015[2:8,"N2"]+xStar2015[2:8,"N12"])/totalNymphs2015
  prevBbNymphs2015[totalNymphs2015 == 0] = 0
  prevBabNymphs2015[totalNymphs2015 == 0] = 0
  
  prevBbNymphs2016 <- (xStar2016[2:i2016,"N1"]+xStar2016[2:i2016,"N12"])/totalNymphs2016    # same for 2016
  prevBabNymphs2016 <- (xStar2016[2:i2016,"N2"]+xStar2016[2:i2016,"N12"])/totalNymphs2016
  prevBbNymphs2016[totalNymphs2016 == 0] = 0
  prevBabNymphs2016[totalNymphs2016 == 0] = 0
  
  prevBbNymphs <- pmax(0, c(prevBbNymphs2014, prevBbNymphs2015, prevBbNymphs2016))          # prevalences for all years combined in single arrays
  prevBabNymphs <- pmax(0, c(prevBabNymphs2014, prevBabNymphs2015, prevBabNymphs2016))             

  burdenLarvae2014 <- xStar2014[2:8,"LB"]              # burdens of nymphs and larvae on mice, 2014
  burdenNymphs2014 <- xStar2014[2:8,"NB"]    
  burdenLarvae2015 <- xStar2015[2:8,"LB"]              # 2015
  burdenNymphs2015 <- xStar2015[2:8,"NB"]
  burdenLarvae2016 <- xStar2016[2:8,"LB"]              # 2016
  burdenNymphs2016 <- xStar2016[2:8,"NB"]
  burdenLarvae <- pmax(0, c(burdenLarvae2014, burdenLarvae2015, burdenLarvae2016))  # combined into one array
  burdenNymphs <- pmax(0, c(burdenNymphs2014, burdenNymphs2015, burdenNymphs2016))
    
  # all output combined into array corresponding exactly to data, ready for direct comparison in abc algorithm
  out <- c(denMice, burdenLarvae, burdenNymphs, prevBbMice, prevBbNymphs, prevBabMice, prevBabNymphs) 
  
  # for general purposes output three years solution, all times and populations as matrix
  # year <- c( rep(2014, length(xStar[,"time"])), rep(2015, length(xStar[,"time"])), rep(2016, length(xStar[,"time"])))
  # day <- rep( xStar[,"time"], 3 )
  # denMice <- rep( xStar[,"M0"]+xStar[,"M1"]+xStar[,"M2"]+xStar[,"M12"], 3 )
  # prevBbMice <- rep( xStar[,"M1"]+xStar[,"M12"], 3 )/denMice
  # prevBabMice <- rep( xStar[,"M2"]+xStar[,"M12"], 3 )/denMice
  # prevBbNymphs <- rep( (xStar[,"N1"]+xStar[,"N12"])/(xStar[,"N0"]+xStar[,"N1"]+xStar[,"N2"]+xStar[,"N12"]), 3 )
  # prevBabNymphs <- rep( (xStar[,"N2"]+xStar[,"N12"])/(xStar[,"N0"]+xStar[,"N1"]+xStar[,"N2"]+xStar[,"N12"]), 3 )
  # burdenLarvae <- rep(xStar[,"LB"], 3)
  # burdenNymphs <- rep(xStar[,"NB"], 3)
  # out <- cbind(year, day, denMice, prevBbMice, prevBabMice, prevBbNymphs, prevBabNymphs, burdenLarvae, burdenNymphs)
  
  return(out)
}
