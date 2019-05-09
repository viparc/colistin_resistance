# Model 3 functions

eq_mod_3=function(t, state, param){
  
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

pred_mod_3 = function(farm, param){
  beta = param[[1]]
  gamma = param[[2]]
  alpha = param[[3]]
  ab_expo_farm = param[[4]]
  pS = param[[5]] # Proportion of resistant bacteria in "S" chickens
  pR = param[[6]] # Proportion of resistant bacteria in "R" chickens
  
  init_prev = round(obs[farm, 1]*200)
  pred_indiv = ode(c(200-init_prev, chick_prev), 1:n_weeks, eq_mod_3, list(beta, gamma, alpha, ab_expo[farm,]))
  
  # Observation process:
  pred = rep(NA, n_weeks)
  for(t in 1:n_weeks){
    pred[t] = obs_process(pred_indiv[t,2:3], list(50, 0, 1))
  }
  
  return(pred)
}
