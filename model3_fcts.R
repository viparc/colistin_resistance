# Model 3 functions

eq_mod=function(t, state, param){
  
  beta = param[[1]]
  gamma = param[[2]]
  alpha = param[[3]]
  ab_expo_farm = param[[4]]
  
  S = state[1]
  R = state[2]
  
  dS = -(beta * (alpha^ab_expo_farm[t]) * S * R /(S+R)) + gamma * R
  dR = (beta * (alpha^ab_expo_farm[t]) * S * R /(S+R)) - gamma * R
  
  return(list(c(dS, dR)))
}

obs_process=function(state, param){
  n_chick_samp = param[[1]] # Number of individuals whose feces are collected
  pS = param[[2]] # Proportion of resistant bacteria in "S" chickens
  pR = param[[3]] # Proportion of resistant bacteria in "R" chickens
  
  S = state[1]
  R = state[2]
  
  n_samp_R = sum(sample(c(rep(1, R), rep(0, S)), n_chick_samp, replace=F)) # Number of "R" individuals whose feces are collected
  n_samp_S = n_chick_samp - n_samp_R # Number of "S" individuals whose feces are collected
  
  # Re: Proportion of resistant bacteria in the pooled sample collected
  Re = pS * n_samp_S/(n_samp_S + n_samp_R) + pR * n_samp_R/(n_samp_S + n_samp_R)
  
  return(Re)
}
