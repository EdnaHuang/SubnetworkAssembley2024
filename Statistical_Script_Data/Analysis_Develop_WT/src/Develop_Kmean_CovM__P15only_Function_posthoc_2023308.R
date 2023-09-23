###for report----------------------------------------------------
library(dplyr)
library(emmeans)
library(parameters)
library(bayestestR)
library(ggplot2)
library(knitr)


# extract the replicates and assign attributes and class ------------------

# extract lmer bootstrapped coefficients ------------------------------------

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


# extract lm bootstrapped coefficients ------------------------------------

lm_boot_extract = function(x) {
  
  # extract simpleboot object
  mod = lm_boot_list[[x]]
  
  # extract bootstrapped coefs
  boot_ls = lapply(1:length(mod$boot.list), 
                   function(x) mod$boot.list[[x]]$coef)
  boot_rep = do.call(rbind, boot_ls)
  
  # keep complete cases
  boot_rep = boot_rep[complete.cases(boot_rep), ]
  
  # assign attributes and class to the replicates
  attr(boot_rep, "original_model") = mod$orig.lm
  class(boot_rep) = append(class(boot_rep), "bootstrap_model")
  class(boot_rep) = append(class(boot_rep), "see_bootstrap_model")
  
  # return replicates
  return(boot_rep)
}


##extract the replicates and assign attributes and class -------
# boot model extract raw function returns two data frames for plotting figure: 
# (1) boot sample and (2) original data

boot_model_extract_raw = function(x) {
  # extract replicates
  if (x %in% names(case_boot_list)) {
    mod = case_boot_list[[x]]
  } else {
    mod = lm_boot_list[[x]]
  }
  
  # save replicates
  mod_rep = mod$replicates
  
  # save original datas
  mod_data = mod$model@frame
  
  # keep complete cases
  mod_rep = mod_rep[complete.cases(mod_rep), ]
  
  # assign attributes and class to the replicates
  attr(mod_rep, "original_model") = mod$model
  class(mod_rep) = append(class(mod_rep), "bootstrap_model")
  class(mod_rep) = append(class(mod_rep), "see_bootstrap_model")
  
  # return replicates
  return(list(mod_rep, mod_data))
}



##plot raw data and bootstrapping result -Facet by layer ------------------

# x = "density_number_100_um_3_vgat_number"

plot_raw_boot_facet_layer = function(x, var_name = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1), 
                                     legend_position = c(0.4, 0.85)) {
  # define label
  layer.label <- c("Layer 2/3", "Layer 4")
  names(layer.label) <- c("2", "4")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)[[1]]
  
  # extract data
  mod_data = boot_model_extract_raw(x)[[2]]
  
  mod_data$yvar = mod_data[, x]
  mod_data$xvar = mod_data$sex
  mod_data$tvar = mod_data$f.layer
  
  emms = emmip(mod_boot, ~ sex | f.layer, CIs = TRUE, plotit = FALSE)
  
  emm_plot = ggplot(emms, aes(x = sex, y = yvar, color = sex)) +
    geom_pointrange(aes(ymin = LCL, ymax = UCL), linewidth = 1) +
    geom_beeswarm(data = mod_data, dodge.width = 0.5, alpha = 0.5, shape = 1) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ f.layer,labeller = labeller(f.layer =layer.label)) +
    labs(x = "", y= "", color = "Sex") +
    theme(strip.text = element_text(size = 10), 
          strip.background = element_rect(color = "black", linewidth = 1), 
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position, legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank())+
     scale_color_nejm()
  
  print(emm_plot)
}


##plot raw data and bootstrapping result -Facet by sex ------------------

plot_raw_boot_facet_sex = function(x, var_name = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1), 
                                     legend_position = c(0.4, 0.85)) {
  # define label
  sex.label <- c("Female", "Male")
  names(sex.label) <- c("F", "M")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)[[1]]
  
  # extract data
  mod_data = boot_model_extract_raw(x)[[2]]
  
  mod_data$yvar = mod_data[, x]
  mod_data$xvar = mod_data$f.layer
  mod_data$tvar = mod_data$ sex
  
  emms = emmip(mod_boot, ~ f.layer | sex, CIs = TRUE, plotit = FALSE)
  
  emm_plot = ggplot(emms, aes(x = f.layer, y = yvar, color = f.layer)) +
    geom_pointrange(aes(ymin = LCL, ymax = UCL), linewidth = 1) +
    geom_beeswarm(data = mod_data, dodge.width = 0.5, alpha = 0.5, shape = 1) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ sex,labeller = labeller(sex =sex.label)) +
    labs(x = "", y= "", color = "layer") +
    theme(strip.text = element_text(size = 10), 
          strip.background = element_rect(color = "black", linewidth = 1), 
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position, legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank())
  
  print(emm_plot)
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
  
  # extract replicates
  #mod_boot = boot_model_extract(x)
  
  # anova of selected model-before bootstrapping
  anova = anova_list_select[[x]]
  
  # emmeans: sex
  emm_sex = emmeans(mod_boot, "sex")
  emm_contr_sex = contrast(emm_sex, method = "pairwise")
  emm_param_sex = as.data.frame(model_parameters(
    emm_contr_sex, centrality = c("median", "mean"), test = "pd"))
  emm_param_sex$pval = pd_to_p(emm_param_sex$pd)
  
  # emmeans: layer
  emm_layer = emmeans(mod_boot, "f.layer")
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
report_interaction = function(x, var_int1 = "f.layer", var_int2 = "sex") {
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

