dat = obtab[,c("amu", "auc", "init", "res")]
dat = na.omit(dat)
# load("C:/Users/Jonathan/Desktop/Data EM.rdata")

# E-step:

e.step <- function(dat, params, type_mod, sigma = 0.1) {

  lin_intro = params[["mu"]] + params[["alpha"]] * dat$amu + params[["gamma"]] * dat$init + params[["eta"]]
  lin_nointro = params[["mu"]] + params[["alpha"]] * dat$amu + params[["gamma"]] * dat$init
  
  if(type_mod == "logistic"){
    
    lin_intro = pmin(200, lin_intro)
    lin_nointro = pmin(200, lin_nointro)
    
    p_res_if_nointro = 10^lin_nointro/(10^lin_nointro + 1)
    p_res_if_nointro = pmin(0.9999, p_res_if_nointro)
    p_res_if_nointro = pmax(0.0001, p_res_if_nointro)
    
    p_res_if_intro = 10^lin_intro/(10^lin_intro + 1)
    p_res_if_intro = pmin(0.9999, p_res_if_intro)
    p_res_if_intro = pmax(0.0001, p_res_if_intro)
    
    p_nores_if_nointro = 1-p_res_if_nointro
    p_nores_if_intro = 1-p_res_if_intro
    
    p_intro = params[["p_intro"]]
    p_intro_if_res = p_res_if_intro * p_intro / (p_res_if_intro * p_intro + p_res_if_nointro * (1-p_intro))
    p_intro_if_nores = p_nores_if_intro * p_intro / (p_nores_if_intro * p_intro + p_nores_if_nointro * (1-p_intro))
    
    prob_intro = dat$res * p_intro_if_res + (1-dat$res) * p_intro_if_nores
    introd = rbinom(n = length(prob_intro), size = 1, prob = prob_intro)
    introd = as.logical(introd)
    
  }else if(type_mod == "linear"){
    
    r_intro = abs(lin_intro - dat$auc) # residuals
    r_nointro = abs(lin_nointro - dat$auc)
    
    exp_intro = exp(- r_intro^2/sigma^2)
    exp_nointro = exp(- r_nointro^2/sigma^2)
    
    exp_intro = pmin(0.9999, exp_intro)
    exp_intro = pmax(0.0001, exp_intro)
    exp_nointro = pmin(0.9999, exp_nointro)
    exp_nointro = pmax(0.0001, exp_nointro)
    
    prob_intro <- exp_intro / (exp_intro + exp_nointro)

    introd = rbinom(n = length(prob_intro), size = 1, prob = prob_intro)
    introd = as.logical(introd)
  }

  return(introd)
}

# M-step:

m.step <- function(dat, introd, type_mod) {
  dat_intr = dat
  dat_intr$intro = introd
  
  prop_intro = sum(dat_intr$intro)/length(dat_intr$intro)
  
  if(type_mod == "logistic"){
    
    modstep = glm(data = dat_intr, res ~ amu + init + intro, family="binomial")
    llmod = as.numeric(logLik(modstep))
    
  }else if(type_mod == "linear"){
    modstep = lm(data = dat_intr, auc ~ amu + init + intro)
    llmod = as.numeric(logLik(modstep))
  }
  
  est_param = as.list(c(summary(modstep)$coefficients[,"Estimate"]))
  if(length(est_param) < 4){est_param[[4]] = 0}
  est_param = c(est_param, prop_intro)
  names(est_param) = c("mu", "alpha", "gamma", "eta", "p_intro")
  
  return(list(est_param, llmod, modstep))
}

# Algo:

em.2lines <- function(dat, tol=1e-6, max.step=1e3, type_mod, init_params = list(mu = runif(1), alpha = runif(1), gamma = runif(1), eta = runif(1), p_intro = runif(1))) {
  step = 0
  loglik = 10^6
  
  params = runif(n=5) # init_params
  names(params) = c("mu", "alpha", "gamma", "eta", "p_intro")
  
  repeat {
    introd = e.step(dat, params, type_mod=type_mod)
    fit_mod = m.step(dat, introd, type_mod=type_mod)
    params = fit_mod[[1]]
    old.loglik = loglik
    loglik = fit_mod[[2]]
    modstep = fit_mod[[3]]

    if (abs(loglik - old.loglik) < tol){
      counting_for_convergence = counting_for_convergence +1
    }else{
      counting_for_convergence = 0
    }
    
    if(counting_for_convergence == 5){break}
    
    step = step +1
    if (step > max.step){
      break
    }
  }
  print(paste0("Final step: ", step))
  
return(list(modstep, params, introd, loglik))
}

rand = round(100*runif(n=1))
set.seed(rand)
results = em.2lines(dat, type_mod = "linear") #, init_params = list(mu = 0.05, alpha = 0.0025, gamma = 0, eta = 1, p_intro = sum(dat$res[dat$amu == 0]) /length(dat$res[dat$amu == 0])))
summary(results[[1]])
# results[[2]]
print(paste0("Probability of introduction: ", results[[2]][["p_intro"]]))
print(paste0("Loglikelihood:", results[[4]]))

# Plot results:

dat_intr = dat
dat_intr$intro = results[[3]]

ggplot(dat_intr, aes(x=amu, y=auc, col=intro, shape=res)) + geom_point()
