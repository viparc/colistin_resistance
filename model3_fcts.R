# Model 3 functions

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


library(adaptivetau)

mod_transitions <- list(
  c(S = -1, R = 1), # infection by resistant bacteria
  c(R = -1, S = 1), # recovery from resistant bacteria
  c(S = -1, R = 0), # death of "S" individuals
  c(R = -1, S = 0) # death of "R" individuals
)

mod_rateFunc <- function(state, parameters, t) {
  
  beta <- parameters[[1]]
  gamma <- parameters[[2]]
  alpha <- parameters[[3]]
  nu <- parameters[[4]]
  ab_expo_farm <- parameters[[5]]
  
  S <- state["S"] 
  R <- state["R"] 
  N <- S + R
  
  t_day = max(ceiling(t), 1)
  
  return(c(
    (alpha ^ ab_expo_farm[t_day]) * beta * S * R / N,
    gamma * R,
    nu * S,
    nu * R
  ))
}


pred_mod_3 = function(this_farm_size, ab_expo_farm, param){
  beta = param[[1]]
  gamma = param[[2]]
  alpha = param[[3]]
  nu = param[[4]]
  eta = param[[5]]
  pS = param[[6]] # Proportion of resistant bacteria in "S" chickens
  pR = param[[7]] # Proportion of resistant bacteria in "R" chickens
  
  init_prev = round(eta*this_farm_size)
  
  pred_ev = data.frame(ssa.adaptivetau(init.values = c(S = this_farm_size-init_prev, R = init_prev),
                               transitions = mod_transitions,
                               rateFunc = mod_rateFunc,
                               params = list(beta, gamma, alpha, nu, ab_expo[farm,]),
                               tf = n_weeks))
  
  pred_S = approx(x = pred_ev$time, y = pred_ev$S, xout = 1:n_weeks, method = "constant")
  pred_R = approx(x = pred_ev$time, y = pred_ev$R, xout = 1:n_weeks, method = "constant")

  # Observation process:
  pred = rep(NA, n_weeks)
  for(t in 1:n_weeks){
    pred[t] = obs_process(c(pred_S$y[t], pred_R$y[t]), list(50, 0, 1))
  }
  
  return(pred)
}
