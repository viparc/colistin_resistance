rm(list=ls(all=TRUE))

library(deSolve)
library(bbmle)
library(ggplot2)
library(ggrepel)
library(grid)
library(reshape2)
library(readxl)
library(pracma)
library(agrmt)
library(shiny)
library(parallel)
library(plotly)
library(sp)
library(ape)
library(rgdal)
library(ggsn)
library(ggforce)
library(dendextend)
library(mixtools)
library(ggpubr)

resrap <- function(conc_used, almost_one){
  
  definition_amu = "allab"
  amu_model_is_quanti = F
  effect_recent_amu = "step"
  
  load(paste0("C:/Users/Jonathan/Desktop/Resultats simulations/colistin_resistance/", definition_amu, " - ", ifelse(amu_model_is_quanti, "Quantitative", "Qualitative"), " - ", effect_recent_amu, "/eta=", almost_one, "/c=", conc_used, "/resulpar_all_", definition_amu, "_isquanti", amu_model_is_quanti, "_", effect_recent_amu, "_c", conc_used, ".rdata"))
  
  l_Q_fct = l_H_fct = l_numb_par = c()
  
  for(m in 1:12){
    l_Q_fct[m] = -as.numeric(resulpar[[m]][[4]])
    
    l_numb_par[m] = length(resulpar[[m]][[1]]@coef)
    
    vec_p = resulpar[[m]][[3]]
    
    l_H_fct[m] = sum(vec_p * log(vec_p) + (1-vec_p) * log(1-vec_p))
  }
  
  l_IC_HQ = -2*l_Q_fct + 2*l_H_fct + 2*l_numb_par
  l_LL_obs = l_Q_fct - l_H_fct
  
  tab_comp_mod = data.frame(mod = as.factor(1:12), Q_fct = l_Q_fct, H_fct = l_H_fct, LL_obs = l_LL_obs, numb_par = l_numb_par, IC_HQ = l_IC_HQ)
  
  best_model = which.min(tab_comp_mod$IC_HQ)
  print(paste0("Best model: model ", best_model))
  tab_comp_mod$IC_HQ <- round(tab_comp_mod$IC_HQ, 2)
  print(tab_comp_mod[,c("mod", "IC_HQ")])
  
  print(c(conc_used, almost_one, round(resulpar[[best_model]][[5]][c("mu","alpha","beta")], 2)))
  
  # val_ICHQ = data.frame(mod = 1:12, IC_HQ = l_IC_HQ, rec_amu = NA, init_amu = NA, first_samp = NA, prev_samp = NA)
  # 
  # gap_minmax = max(l_IC_HQ) - min(l_IC_HQ)
  # val_ICHQ[7:12, "rec_amu"] = min(l_IC_HQ) - gap_minmax/20
  # val_ICHQ[c(4:6, 10:12), "init_amu"] = min(l_IC_HQ) - 2 * gap_minmax/20
  # val_ICHQ[seq(2,11,3), "first_samp"] = min(l_IC_HQ) - 3 * gap_minmax/20
  # val_ICHQ[seq(3,12,3), "prev_samp"] = min(l_IC_HQ) - 4 * gap_minmax/20
  # 
  # val_ICHQ = val_ICHQ[order(val_ICHQ$IC_HQ, decreasing = F),]
  # val_ICHQ$mod = factor(val_ICHQ$mod, levels = val_ICHQ$mod)
  # 
  # val_ICHQ = melt(val_ICHQ, id.vars=c("mod", "IC_HQ"), variable.name="var", value.name="position")
  # 
  # p = ggplot(data = val_ICHQ, aes(x = mod, y = IC_HQ, group = 1, na.rm = TRUE))
  # p = p + xlab("Models") + ylab(bquote(IC['H,Q']))
  # p = p + theme_bw()
  # p = p + geom_line() + geom_point(size = 2)
  # p = p + geom_segment(aes(x = as.numeric(mod)-0.5, xend = as.numeric(mod)+0.5, y = position, yend = position, col = var), size = 3)
  # p = p + scale_colour_brewer(type="qual", palette=1,
  #                             name = "Variables included:",
  #                             labels = c("Recent AMU", "Initial AMU", "First sample", "Previous sample"))
  # p
}

resrap(1, 0.999)
resrap(2, 0.999)
resrap(4, 0.999)
resrap(1, 0.99)
resrap(2, 0.99)
resrap(4, 0.99)

###################################################

# Plots

load("C:/Users/Jonathan/Desktop/Resultats simulations/colistin_resistance/maps_robustness/plot_eta=0.999_c=2.R")
p_0.999_2 = pB + theme(legend.position = "none")

load("C:/Users/Jonathan/Desktop/Resultats simulations/colistin_resistance/maps_robustness/plot_eta=0.999_c=4.R")
p_0.999_4 = pB + theme(legend.position = "none")

load("C:/Users/Jonathan/Desktop/Resultats simulations/colistin_resistance/maps_robustness/plot_eta=0.99_c=1.R")
com_leg = get_legend(pB)
p_0.99_1 = pB + theme(legend.position = "none")

load("C:/Users/Jonathan/Desktop/Resultats simulations/colistin_resistance/maps_robustness/plot_eta=0.99_c=2.R")
p_0.99_2 = pB + theme(legend.position = "none")

load("C:/Users/Jonathan/Desktop/Resultats simulations/colistin_resistance/maps_robustness/plot_eta=0.99_c=4.R")
p_0.99_4 = pB + theme(legend.position = "none")


plot.list <- lapply(list(p_0.999_2, p_0.999_4, p_0.99_1, p_0.99_2), 
                    function(p) p + theme(plot.background = element_rect(color = "black")))

all_plot = ggarrange(p_0.999_2, p_0.999_4, p_0.99_1, p_0.99_2, p_0.99_4, com_leg, nrow = 3, ncol = 2, font.label = list(size = 16), labels = c("(2; 0.999)", "(4; 0.999)", "(1; 0.99)", "(2; 0.99)", "(4; 0.99)")) +
  annotation_custom(grid.polygon(c(0, 0.5, 1, 0, 0.5, 1, 0.5, 0.5, 0.5),
                                   c(0.666, 0.666, 0.666, 0.333, 0.333, 0.333, 0, 0.5, 1), 
                                   id = c(1,1,1,2,2,2,3,3,3), 
                                   gp = gpar(lwd = 1.5)))


plot(all_plot)

ggsave("C:/Users/Jonathan/Desktop/maps_robustness.png", all_plot, width = 20, height = 32, units = "cm")
