# dat = obtab[,c("amu", "auc", "init", "res")]
# dat = na.omit(dat)
# load("C:/Users/Jonathan/Desktop/Data EM.rdata")

# E-step:

e.step <- function(dat=dat, params, type_mod, sigma = 0.15) {
  
  if(type_mod == "logistic"){
    
    lin_intro = params[["mu"]] + params[["alpha"]] * dat$amu + params[["gamma"]] * dat$init + params[["eta"]]
    lin_nointro = params[["mu"]] + params[["alpha"]] * dat$amu + params[["gamma"]] * dat$init
    
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
    
    lin_intro = params[["mu"]] + params[["alpha"]] * dat$amu + params[["gamma"]] * dat$init + params[["eta"]]
    lin_nointro = params[["mu"]] + params[["alpha"]] * dat$amu + params[["gamma"]] * dat$init
    
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
    
  }else if(type_mod == "linear2"){
    
    lin_intro = run2(aim="pred", intr=1, mu=params[["mu"]], gamma=params[["gamma"]], alpha=params[["alpha"]], teta=params[["teta"]], delta=params[["delta"]], beta=params[["beta"]], eta=params[["eta"]])[["pred"]]
    lin_nointro = run2(aim="pred", intr=0, mu=params[["mu"]], gamma=params[["gamma"]], alpha=params[["alpha"]], teta=params[["teta"]], delta=params[["delta"]], beta=params[["beta"]], eta=logit(0))[["pred"]]
    
    r_intro = abs(lin_intro - obtab$auc) # residuals
    r_nointro = abs(lin_nointro - obtab$auc)
    
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

m.step <- function(dat=dat, introd, type_mod) {
  
  prop_intro = sum(introd)/length(introd)

  if(type_mod == "logistic"){
    
    dat_intr = dat
    dat_intr$intro = introd
    
    modstep = glm(data = dat_intr, res ~ amu + init + intro)
    llmod = as.numeric(logLik(modstep))
    
    est_param = as.list(c(summary(modstep)$coefficients[,"Estimate"]))
    if(length(est_param) < 4){est_param[[4]] = 0}
    
    est_param = c(est_param, prop_intro)
    names(est_param) = c("mu", "alpha", "gamma", "eta", "p_intro")
    true_coef = est_param
    
  }else if(type_mod == "linear"){
    
    dat_intr = dat
    dat_intr$intro = introd
    
    modstep = lm(data = dat_intr, auc ~ amu + init + intro)
    llmod = as.numeric(logLik(modstep))
    
    est_param = as.list(c(summary(modstep)$coefficients[,"Estimate"]))
    if(length(est_param) < 4){est_param[[4]] = 0}
    
    est_param = c(est_param, prop_intro)
    names(est_param) = c("mu", "alpha", "gamma", "eta", "p_intro")
    true_coef = est_param
    
  }else if(type_mod == "linear2"){
    
    Sd = 0.15
    run2 = function(aim="est", intr=introd, mu, gamma, alpha, teta, delta, beta, eta){ #, lambda){
      
      mu = sigmoid(mu)
      gamma = gamma
      alpha = alpha
      teta = exp(teta) +1
      delta = 10 * sigmoid(delta) +1
      beta = 10 * sigmoid(beta) +1
      eta = sigmoid(eta)
      # lambda = sigmoid(lambda)
      
      obsdat = aucamu(nwkseff=beta, quantuse=T, tempo_amu="exp_decay", thresh=0.4, logtrans=F, ignore.nul=F, decay_init=0, teta=teta, delta=delta, obj="est")
      obsdat$introd = intr
      
      predmod = mu + alpha * obsdat$amu + eta * obsdat$introd + gamma * obsdat$init #+ lambda^(obsdat$age <= delta)
      
      if(aim == "est"){
        return(-sum(dnorm(x=obsdat$auc, mean=predmod, sd=Sd, log=T), na.rm=T))
        # return(-sum(dbeta(x=obsdat$auc, shape1=shape, shape2=shape*(1.00001-predmod)/(predmod+0.00001), log=T), na.rm=T))
      }else if (aim == "pred"){
        return(list(obs = obsdat$auc, pred = predmod, minll = -sum(dnorm(x=obsdat$auc, mean=predmod, sd=Sd, log=T), na.rm=T)))
      }
    }
    
    fit_mle2 = list()
    llmod = list()
    for(rep in 1:3){
      print(paste0("M-step repeat ", rep))
      fit_mle2[[rep]] = mle2(minuslogl = run2, method = "Nelder-Mead", start = list(mu = logit(runif(1)), gamma = runif(1), alpha = runif(1), teta = log(runif(1)), delta = logit(runif(1)), beta = logit(runif(1)), eta = logit(runif(1))))
      true_coef = coef(fit_mle2[[rep]])
      llmod[[rep]] = run2(aim="pred", mu=true_coef[1], gamma=true_coef[2], alpha=true_coef[3], teta=true_coef[4], delta=true_coef[5], beta=true_coef[6], eta=true_coef[7])[[3]]
    }
    fit_mle2 = fit_mle2[[which.min(llmod)]]
    llmod = llmod[[which.min(llmod)]]
    
    # fit_mle2 = mle2(minuslogl = run2, method = "Nelder-Mead", start = list(mu = logit(runif(1)), gamma = runif(1), alpha = runif(1), teta = log(runif(1)), delta = logit(runif(1)), beta = logit(runif(1)), eta = logit(runif(1))))
    # true_coef = coef(fit_mle2)
    # llmod = run2(aim="pred", mu=true_coef[1], gamma=true_coef[2], alpha=true_coef[3], teta=true_coef[4], delta=true_coef[5], beta=true_coef[6], eta=true_coef[7])[[3]]
    
    modstep = fit_mle2
    
    true_coef = coef(fit_mle2)
    true_coef = c(true_coef, p_intro = prop_intro)
    est_param = true_coef
    
    true_coef[1] = sigmoid(true_coef[1])
    true_coef[4] = exp(true_coef[4]) +1
    true_coef[5] = 10 * sigmoid(true_coef[5]) +1
    true_coef[6] = 10 * sigmoid(true_coef[6]) +1
    true_coef[7] = sigmoid(true_coef[7])
    
  }
  
  return(list(est_param, llmod, modstep, true_coef))
}

# Algo:

em.2lines <- function(dat=0, tol=1e-6, max.step=1e3, type_mod, init_params) {
  step = 0
  loglik = 10^6
  
  params = init_params

  repeat {
    print(paste0("Step: ", step))
    introd = e.step(dat, params, type_mod=type_mod)
    fit_mod = m.step(dat, introd, type_mod=type_mod)
    params = fit_mod[[1]]
    old.loglik = loglik
    loglik = fit_mod[[2]]
    print(paste0("MinusLL: ", loglik))
    modstep = fit_mod[[3]]
    param_true_format = fit_mod[[4]]

    if (abs(loglik - old.loglik) < tol){
      counting_for_convergence = counting_for_convergence +1
    }else{
      counting_for_convergence = 0
    }
    
    if((counting_for_convergence == 1) | (loglik < -30)){break}
    
    step = step +1
    if (step > max.step){
      break
    }
  }
  print(paste0("Final step: ", step))
  
return(list(modstep, params, introd, loglik, param_true_format))
}

rand = round(100*runif(n=1))
set.seed(rand)
results = em.2lines(type_mod = "linear2", init_params = list(mu = logit(runif(1)), gamma = runif(1), alpha = runif(1), teta = log(runif(1)), delta = logit(runif(1)), beta = logit(runif(1)), eta = logit(runif(1)), p_intro = runif(1)))
summary(results[[1]])
results[[5]]
print(paste0("Probability of introduction: ", results[[2]][["p_intro"]]))
print(paste0("Loglikelihood:", results[[4]]))

# Plot results:

obs_vs_pred = run2(aim="pred", intr=results[[3]], mu=results[[2]][1], gamma=results[[2]][2], alpha=results[[2]][3], teta=results[[2]][4], delta=results[[2]][5], beta=results[[2]][6], eta=results[[2]][7]) #, lambda=coef(fit_mle2)[8])
obs_vs_pred = data.frame(obs = obs_vs_pred[[1]], pred = obs_vs_pred[[2]])

p = ggplot(obs_vs_pred, aes(x = obs, y = pred))
p = p + xlab("Obs") + ylab("Pred")
p = p + ggtitle(paste0("Observed VS Predicted AUC (minus LogLikelihood =", round(results[[4]], 2), ")"))
p = p + geom_segment(aes(xend = obs, y = pmax(0, obs_vs_pred[[2]] - 1.96*Sd), yend = pmin(1, obs_vs_pred[[2]] + 1.96*Sd)), col = "grey")
p = p + geom_abline(intercept = 0, slope = 1, linetype = "dashed")
p = p + geom_point(aes(col = as.logical(c(rep(0, 4), 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, rep(0, 10), 1, 0, 0, 1, rep(0, 8), 1, rep(0,3))))) + labs(col = "Introduction")
p

dat_intr = aucamu(obj = "est", nwkseff=results[[5]][6], teta=results[[5]][4], delta=results[[5]][5])
dat_intr$intro = results[[3]]

ggplot(dat_intr, aes(x=amu, y=auc, col=intro, shape=res)) + geom_point()
