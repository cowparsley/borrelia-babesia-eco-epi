# function that computes values of derivatives at time t for ode model of mouse and tick population dynamics
# inputs are time, state variables and parameters

mouseTickBbBabODE <- function(t, stateVars, params) {
  # t is current time, stateVars is current value of state variables, params is list of parameters
  
  with(as.list(c(stateVars, params)), {   # use the state variables stateVars and parameters params

    params[params < 0] = 0
    
    # precompute some composite parameters
    xiBetaML1 = min(xi*betaML1, 1)
    sigBetaML2 = min(sigma*betaML2, 1)
    alphaBetaNM1 = min(alpha*betaNM1, 1) 
    
    # set switches indicating whether emergence is occurring
    if(t > tauE) { switchE = 1} else { switchE = 0 }                # switch on/off larval emergence from overwintered questing larvae
    if(t > tauL) { switchL = 1} else { switchL = 0 }                # switch on/off larval emergence from eggs
    if(t > tauN) { switchN = 1} else { switchN = 0 }                # switch on/off nymph emergence from overwintered engorged larvae
    
    # specify the differential equations for continuous time dyanmics during the active season
    # mice
    allMice = M0 + M1 + M2 + M12
    dM0 = r*(M0 + M1 + (1 - nu)*(M2 + M12))*(1 - allMice/K) - mu*M0 + gamma1*M1 - lambda*M0*(betaNM1*N1 + betaNM2*N2 + (betaNM1 + betaNM2 - betaNM1*betaNM2)*N12) # uninfected 
    dM1 = lambda*M0*(betaNM1*N1 + betaNM1*(1 - betaNM2)*N12) - lambda*M1*betaNM2*(N2 + N12) - (gamma1 + mu)*M1                                     # borrelia infected 
    dM2 = r*nu*(M2 + M12)*(1 - allMice/K) + lambda*M0*(betaNM2*N2 + (1-betaNM1)*betaNM2*N12) + gamma1*M12 - lambda*M2*alphaBetaNM1*(N1 + N12) - mu*M2                                          # babesia infected 
    dM12 = lambda*(M0*betaNM1*betaNM2*N12 + M1*betaNM2*(N2 + N12) + M2*alphaBetaNM1*(N1 + N12)) - (gamma1 + mu)*M12                                                                # coinfected     
    
    # questing larvae (all uninfected)
    dL = switchE*initE*etaE*exp(-etaE*(t - tauE)) + switchL*Omega*etaL*exp(-etaL*(t-tauL)) - lambda*L*(allMice + D)  
    
    # engorged overwintering larvae
    dW0 = lambda*L*(M0 + (1-betaML1)*M1 + (1-betaML2)*M2 + (1-xiBetaML1)*(1-sigBetaML2)*M12 + D ); # uninfected 
    dW1 = lambda*L*(betaML1*M1 + xiBetaML1*(1-sigBetaML2)*M12)                                     # borrelia infected
    dW2 = lambda*L*(betaML2*M2 + (1-xiBetaML1)*sigBetaML2*M12)                                     # babesia infected
    dW12 = lambda*L*xiBetaML1*sigBetaML2*M12                                                       # co-infected
    
    # questing nymphs
    dN0 = switchN*initN0*etaN*exp(-etaN*(t - tauN)) - lambda*N0*(allMice + D)                       # uninfected
    dN1 = switchN*initN1*etaN*exp(-etaN*(t - tauN)) - lambda*N1*(allMice + D)                       # borrelia infected
    dN2 = switchN*initN2*etaN*exp(-etaN*(t - tauN)) - lambda*N2*(allMice + D)                       # babesia infected
    dN12 = switchN*initN12*etaN*exp(-etaN*(t - tauN)) - lambda*N12*(allMice + D)                    # co-infected
        
    # burdens   
    dLB = -delta*LB + lambda*L                                              # larvae
    dNB = -delta*NB + lambda*(N0 + N1 + N2 + N12)                           # nymphs
    
    
    # return value is list of derivatives, in same order as state variables in y
    return(list(c(dM0, dM1, dM2, dM12, dL, dW0, dW1, dW2, dW12, dN0, dN1, dN2, dN12, dLB,dNB)))
  })
}