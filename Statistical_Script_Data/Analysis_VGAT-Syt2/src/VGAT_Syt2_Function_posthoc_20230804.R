###for report----------------------------------------------------
library(dplyr)
library(emmeans)
library(parameters)
library(bayestestR)
library(ggplot2)
library(knitr)


##defince function
# extract lm bootstrapped coefficients ------------------------------------

lm_boot_extract = function(x) {
  
  # extract simpleboot object
  mod = lm_boot_list[[x]]
  
  # extract bootstrapped coefs
  boot_ls = lapply(1:length(mod$boot.list), 
                   function(x) mod$boot.list[[x]]$coef)
  boot_rep = as.data.frame(do.call(rbind, boot_ls))
  
  # keep complete cases
  boot_rep = boot_rep[complete.cases(boot_rep), ]
  
  # assign attributes and class to the replicates
  attr(boot_rep, "original_model") = mod$orig.lm
  class(boot_rep) = append(class(boot_rep), "bootstrap_model")
  class(boot_rep) = append(class(boot_rep), "see_bootstrap_model")
  
  # return replicates
  return(boot_rep)
}


# extract the replicates and assign attributes and class ------------------

lmer_boot_extract = function(x) {
  # extract lmeresamp objects
  mod = case_boot_list[[x]]
  
  # save replicates
  boot_rep = mod$replicates
  
  # keep complete cases
  boot_rep = boot_rep[complete.cases(boot_rep), ]
  
  # assign attributes and class to the replicates
  attr(boot_rep, "original_model") = mod$model
  class(boot_rep) = append(class(boot_rep), "bootstrap_model")
  class(boot_rep) = append(class(boot_rep), "see_bootstrap_model")
  
  # return replicates
  return(boot_rep)
}


boot_model_extract_raw = function(x) {
  # extract replicates
  if (x %in% names(case_boot_list)) {
    boot_rep = lmer_boot_extract(x)
  } else {
    boot_rep = lm_boot_extract(x)
  }
  
  return(boot_rep)
}



plot_raw_boot_facet_layer = function(x, var_name = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1),
                                     title = x,
                                     legend_position = "none") {
  # define label
  #layer.label = c("Layer 2/3", "Layer 4")
  layer.label = c("L2/3", "L4")
  names(layer.label) = c("L2/3", "L4")
  
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)
  
  # use original data set instead
  # https://tibble.tidyverse.org/reference/subsetting.html
  # the subsetting behavior is somewhat different in tibbles lol
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  
  # marginal effect point estimates and CI 
  emms = emmip(mod_boot, ~ sex | layer, CIs = TRUE, plotit = FALSE)
  emm_plot = ggplot(emms, aes(x = sex, y = yvar, color = sex)) +
    geom_pointrange(aes(ymin = LCL, ymax = UCL), size = 1,
                    position = position_dodge (width =0.5),
                    shape = 16, fatten = 1, alpha = 1, linewidth = 1) +
    geom_beeswarm(data = df_lite, dodge.width = 0.5, alpha = 0.4,
                  shape = 1, size = 1) +  
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ layer,labeller = labeller(layer =layer.label),
               strip.position = "bottom") +
    labs(x = "", y= "", color = "sex") +
    ggtitle(title)+
    theme(strip.text = element_text(size = 10), 
          strip.background = element_blank(),
          #strip.background = element_rect(color = "black", linewidth = 1), 
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position, legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(1, 0, 0, 0), 'cm'))->p
    #scale_color_nejm() 
  
  return(p)
  
}

plot_raw_boot_facet_sex = function(x, var_name = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1),
                                     title = x,
                                     legend_position = "none") {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("L2", "L4")
  
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)
  
  # use original data set instead
  # https://tibble.tidyverse.org/reference/subsetting.html
  # the subsetting behavior is somewhat different in tibbles lol
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  
  # marginal effect point estimates and CI 
  emms = emmip(mod_boot, ~ layer | sex, CIs = TRUE, plotit = FALSE)
  emm_plot = ggplot(emms, aes(x = layer, y = yvar, color = layer)) +
    geom_pointrange(aes(ymin = LCL, ymax = UCL), size = 1,
                    position = position_dodge (width =0.5),
                    shape = 16, fatten = 1, alpha = 1, linewidth = 1) +
    geom_beeswarm(data = df_lite, dodge.width = 0.5, alpha = 0.4,
                  shape = 1, size = 1) +  
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ sex,labeller = labeller(sex =sex.label),
               strip.position = "bottom") +
    labs(x = "", y= "", color = "layer") +
    ggtitle(title) +
    theme(strip.text = element_text(size = 10), 
          strip.background = element_blank(),
          #strip.background = element_rect(color = "black", linewidth = 1), 
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position, legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'cm')) ->p
    
  return(p)
  
}



# report main effect ------------------------------------------------------
emm_var_main = function(var) {
  emm = emmeans(mod_boot, var)
  emm_contr = contrast(emm, method = "consec")
  emm_param = as.data.frame(model_parameters(
    emm_contr, centrality = c("median", "mean"), test = "pd"))
  emm_param$pval = pd_to_p(emm_param$pd)
  return(emm_param)
}



## report main of sex, layer, age to compare ANOVA table before boots -------------
# x = "area"
report_main_effect = function(x, var_name = x) {
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  # anova of selected model-before bootstrapping
  anova = anova_list_select[[x]]
  
  # emmeans: sex
  emm_sex = emmeans(mod_boot, "sex")
  emm_contr_sex = contrast(emm_sex, method = "pairwise")
  emm_param_sex = as.data.frame(model_parameters(
    emm_contr_sex, centrality = c("median", "mean"), test = "pd"))
  emm_param_sex$pval = pd_to_p(emm_param_sex$pd)
  
  # emmeans: layer
  emm_layer = emmeans(mod_boot, "layer")
  emm_contr_layer = contrast(emm_layer, method = "pairwise")
  emm_param_layer = as.data.frame(model_parameters(
    emm_contr_layer, centrality = c("median", "mean"), test = "pd"))
  emm_param_layer$pval = pd_to_p(emm_param_layer$pd)
  
  table1 = kable(anova,caption = "ANOVA table before bootstrapping", digits = 3)
  print(table1)
  
  table2 = kable(list(emm_param_sex, emm_param_layer), digits = 3)
  print(table2)
}



# report interaction effect -----------------------------------------------
# x = "duration"
report_interaction = function(x, var_int1 = "layer", var_int2 = "sex") {
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  f = paste(var_int1, "*", var_int2)
  emm = emmeans(mod_boot, eval(bquote(~ .(f))))
  emm_contr = contrast(emm, "pairwise")
  emm_param = as.data.frame(model_parameters(
    emm_contr, centrality = c("median", "mean"), test = "pd"))
  emm_param$pval = pd_to_p(emm_param$pd)
  
  table = kable(emm_param, caption = "Complete comparison", digits = 3)
  print(table)
}

