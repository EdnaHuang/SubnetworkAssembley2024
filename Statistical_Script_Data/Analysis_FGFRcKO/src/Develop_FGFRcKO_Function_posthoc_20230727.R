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

##plot raw data and bootstrapping result -Facet by layer ------------------

plot_raw_boot_facet_layer = function(x, var_name = x,title = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1), 
                                     legend_position = "none") {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
   # extract replicates
  mod_boot = boot_model_extract_raw(x)
  
  # use original data set instead
  # https://tibble.tidyverse.org/reference/subsetting.html
  # the subsetting behavior is somewhat different in tibbles lol
  df_lite$yvar = df_lite[, x, drop = TRUE]
 
  # marginal effect point estimates and CI 
  emms = emmip(mod_boot, group ~ f.age | f.layer, CIs = TRUE, plotit = FALSE)
  emm_plot = ggplot(emms, aes(x = f.age, y = yvar, color = group)) +
    geom_pointrange(aes(ymin = LCL, ymax = UCL), size = 1,
                    position = position_dodge (width =0.5),
                    shape = 16, fatten = 1) +
    geom_beeswarm(data = df_lite, dodge.width = 0.5, alpha = 0.7,
                  shape = 1, size = 1) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ f.layer,labeller = labeller(f.layer =layer.label),
               strip.position = "bottom") +
    labs(x = "", y= "", color = "group") +
    ggtitle(title)+
    theme(strip.text = element_text(size = 10), 
          strip.background = element_blank(),
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position, legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'cm'))+
    scale_color_d3() ->p
  
  return(p)
  
}

## custom contrasts --------------------------------------------------------

pw_emm = emmeans:::revpairwise.emmc(1:5)
as.matrix(names(pw_emm))
pw_emm_contr = pw_emm[, c(1, 3, 6, 10, 2, 9)]
names(pw_emm_contr) = c("age13 - age11", "age15 - age13", "age18 - age15", 
                        "age21 - age18", "age15 - age11", "age21 - age15")

## report main of group, layer, age to compare ANOVA table before boots -------------
#x = "synchrony_freq"
#var_name = x

report_main_effect = function(x, var_name = x) {
  
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  # anova of selected model-before bootstrapping
  anova = anova_list_select[[x]]
  
  # emmeans: age
  emm_age = emmeans(mod_boot, "f.age")
  #emm_contr_age = contrast(emm_age, method = pw_emm_contr)
  emm_contr_age = contrast(emm_age, method = "pairwise")
  emm_param_age = as.data.frame(model_parameters(
    emm_contr_age, centrality = c("median", "mean"), test = "pd"))
  emm_param_age$pval = pd_to_p(emm_param_age$pd)
  
  # emmeans: group
  emm_group = emmeans(mod_boot, "group")
  emm_contr_group = contrast(emm_group, method = "pairwise")
  emm_param_group = as.data.frame(model_parameters(
    emm_contr_group, centrality = c("median", "mean"), test = "pd"))
  emm_param_group$pval = pd_to_p(emm_param_group$pd)
  
  # emmeans: layer
  emm_layer = emmeans(mod_boot, "f.layer")
  emm_contr_layer = contrast(emm_layer, method = "pairwise")
  emm_param_layer = as.data.frame(model_parameters(
    emm_contr_layer, centrality = c("median", "mean"), test = "pd"))
  emm_param_layer$pval = pd_to_p(emm_param_layer$pd)
  
  table1 = kable(anova,caption = "ANOVA table before bootstrapping", digits = 3)
  print(table1)
  
  table2 = kable(list(emm_param_group, emm_param_age, emm_param_layer), digits = 3)
  print(table2)
}




## report group inter layer ---------
report_group_int_layer = function(x, var_name = x) {
  
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  # simple comparisons
  emm = emmeans(mod_boot, ~ group* f.layer)
  
  #compare
  emm_contr_group = contrast(emm, method = "pairwise")
  emm_param_group = as.data.frame(model_parameters(
    emm_contr_group, centrality = c("median", "mean"), test = "pd"))
  emm_param_group$pval = pd_to_p(emm_param_group$pd)
  emm_param_group$Median.1 = NULL
  
  # print formatted tables in sequence
  table = kable (emm_param_group, digits = 3)
  print(table)
}


## report age-group-layer -----the most comprehensive comparison--------------
report_inter_age_group_layer = function(x, var_name = x) {
  
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  # simple comparisons
  emm = emmeans(mod_boot, ~ group* f.layer*f.age)
  
  # age
  emm_contr_age = contrast(emm, method = "pairwise", simple = "f.age")
  #emm_contr_age = contrast(emm, method = pw_emm_contr, simple = "f.age")
  emm_param_age = as.data.frame(model_parameters(
    emm_contr_age, centrality = c("median", "mean"), test = "pd"))
  emm_param_age$pval = pd_to_p(emm_param_age$pd)
  emm_param_age$Median.1 = NULL
  
  # group
  emm_contr_group = contrast(emm, method = "consec", simple = "group")
  emm_param_group = as.data.frame(model_parameters(
    emm_contr_group, centrality = c("median", "mean"), test = "pd"))
  emm_param_group$pval = pd_to_p(emm_param_group$pd)
  emm_param_group$Median.1 = NULL
  
  # layer
  emm_contr_layer = contrast(emm, method = "consec", simple = "f.layer")
  emm_param_layer = as.data.frame(model_parameters(
    emm_contr_layer, centrality = c("median", "mean"), test = "pd"))
  emm_param_layer$pval = pd_to_p(emm_param_layer$pd)
  emm_param_layer$Median.1 = NULL
  
  # print formatted tables in sequence
  table = kable (list(emm_param_group,emm_param_age,emm_param_layer), digits = 3)
  print(table)
}


## report age-group-layer -----the most comprehensive comparison-simplified--------------
report_inter_age_group_layer_V2 = function(x, var_name = x) {
  
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  # simple comparisons
  emm = emmeans(mod_boot, ~ group* f.layer*f.age)
  
  # age
  emm_contr_age = contrast(emm, method = "pairwise", simple = "f.age")
  #emm_contr_age = contrast(emm, method = pw_emm_contr, simple = "f.age")
  emm_param_age = as.data.frame(model_parameters(
    emm_contr_age, centrality = c("median", "mean"), test = "pd"))
  emm_param_age$pval = pd_to_p(emm_param_age$pd)
  emm_param_age$Median.1 = NULL
  
  # group
  #emm_contr_group = contrast(emm, method = "consec", simple = "group")
  #emm_param_group = as.data.frame(model_parameters(
  #  emm_contr_group, centrality = c("median", "mean"), test = "pd"))
  #emm_param_group$pval = pd_to_p(emm_param_group$pd)
  #emm_param_group$Median.1 = NULL
  
  # layer
  #emm_contr_layer = contrast(emm, method = "consec", simple = "f.layer")
  #emm_param_layer = as.data.frame(model_parameters(
  #  emm_contr_layer, centrality = c("median", "mean"), test = "pd"))
  #emm_param_layer$pval = pd_to_p(emm_param_layer$pd)
  #emm_param_layer$Median.1 = NULL
  
  # print formatted tables in sequence
  #table = kable (list(emm_param_group,emm_param_age,emm_param_layer), digits = 3)
  table = kable (list(emm_param_group), digits = 3)
  print(table)
}





##merge layer then plot raw data and bootstrapping result------------------

merge_layer_plot_raw_boot = function(x, var_name = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1), 
                                     legend_position = "none") {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)
  
  # use original data set instead
  # https://tibble.tidyverse.org/reference/subsetting.html
  # the subsetting behavior is somewhat different in tibbles lol
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  # marginal effect point estimates and CI 
  emms = emmip(mod_boot, group ~ f.age, CIs = TRUE, plotit = FALSE)
  emm_plot = ggplot(emms, aes(x = f.age, y = yvar, color = group, shape = 1, size = 1)) +
    geom_pointrange(aes(ymin = LCL, ymax = UCL), size = 1,
                    position = position_dodge (width =0.5),
                    shape = 16, fatten = 1) +
    geom_beeswarm(data = df_lite, dodge.width = 0.5, alpha = 0.7,
                  shape = 1, size = 1) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    labs(x = "", y= "", color = "group") +
    theme(strip.text = element_text(size = 10), 
          strip.background = element_blank(),
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position, legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'cm'))+
    scale_color_d3() ->p
  
  return(p)
  
}



## Merge layer, then report interaction between age and group---------

merge_layer_report_age_int_group = function(x, var_name = x) {
  
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  # simple comparisons
  emm = emmeans(mod_boot, ~ f.age * group)
  
  # age
  emm_contr_age = contrast(emm, method = "pairwise", simple = "f.age")
  #emm_contr_age = contrast(emm, method = pw_emm_contr, simple = "f.age")
  emm_param_age = as.data.frame(model_parameters(
    emm_contr_age, centrality = c("median", "mean"), test = "pd"))
  emm_param_age$pval = pd_to_p(emm_param_age$pd)
  emm_param_age$Median.1 = NULL
  
  # group
  emm_contr_group = contrast(emm, "consec", simple = "group")
  emm_param_group = as.data.frame(model_parameters(
    emm_contr_group, centrality = c("median", "mean"), test = "pd"))
  emm_param_group$pval = pd_to_p(emm_param_group$pd)
  emm_param_group$Median.1 = NULL
  
  # print formatted tables in sequence
  table = kable (list(emm_param_group, emm_param_age), digits = 3)
  print(table)
}





