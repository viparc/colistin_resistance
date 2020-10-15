mod_fitted = resulpar[[best_model]][[1]]
valpar = mod_fitted@fullcoef
fixedvar = names(mod_fitted@fullcoef) [! names(mod_fitted@fullcoef) %in% names(mod_fitted@coef)]

conf_int = confint(mod_fitted, level = 0.95)
conf_int["beta",] = 5*sigmoid(conf_int["beta",])
print(conf_int)

prof = profile(mod_fitted)
plot(prof)
prof@profile

plot(type="l", prof@profile$beta$par.vals[,"beta"], abs(prof@profile$beta$z))


AA = prof


plot(type="l", 5*sigmoid(AA@profile$beta$par.vals[,"beta"]), abs(AA@profile$beta$z))


AA@profile$beta = AA@profile$beta[abs(AA@profile$beta$par.vals[,"beta"]) < 4.8,]
plot(AA)

AA@profile$beta$par.vals[,"beta"] = 5*sigmoid(AA@profile$beta$par.vals[,"beta"])
plot(AA)
