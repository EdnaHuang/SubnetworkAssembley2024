###for report----------------------------------------------------
library(dplyr)
library(emmeans)
library(parameters)
library(bayestestR)
library(ggplot2)
library(knitr)

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



##plot raw data
# var = "silhs_mean_before_stat"

plot_raw = function(var, title = NULL, data = df_lite,
                    y_breaks = c(1:6 * 0.2),
                    y_limits = c(0, 1),
                    legend_position = c(0.3, 0.8)) {
  # define label
  layer.label = c("L2/3", "L4")
  names(layer.label) = c("2", "4")
  sex.label = c("F","M")
  names(sex.label)= c("F","M")
  # define color....
  #color_sex = MetBrewer::met.brewer("Derain", 2)
  df_lite %>%
    # create a variable to denote the one-on-one correspondence 
    # between observations based on subject id, layer, and location
    unite(col = "set", c("subject_id", "f.layer", "location"), remove = FALSE) %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "+")), y = get(var), 
               group = set, color = f.layer)) +
    geom_point(aes(shape = p15_cno), lwd = 0.5, size = 1.5, alpha = 0.7) +
    geom_line(linewidth = 0.5, alpha = 1) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~sex*f.layer, labeller = labeller(f.layer =layer.label, sex = sex.label),
               strip.position = "bottom", nrow = 1)+
    labs(x = "", y = "", title = title) +
    scale_x_discrete(labels = c("-", "+")) +
    scale_shape_manual(values = c("neg" = 1, "CNO" = 19)) +
    #scale_color_manual(values = color_sex)+
    theme_classic() +
    theme(legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size  = "10"),
          strip.background = element_blank(),
          #strip.background = element_rect(color = "black", linewidth = 1),
          strip.placement = ("outside"),
          strip.text = element_text(size = 10),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2),
          axis.ticks.x = element_blank()) ->p
  return(p)
}


##plot raw data
# var = "silhs_mean_before_stat"
plot_raw_facet_layer = function(var, title = NULL, data = df_lite,
                                y_breaks = c(1:6 * 0.2),
                                y_limits = c(0, 1),
                                legend_position = c(0.3, 0.8)) {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  sex.label = c("Female","Male")
  names(sex.label)= c("F","M")
  
  #color_sex = MetBrewer::met.brewer("Hiroshige", 2)
  #color_sex = MetBrewer::met.brewer("Derain", 2)
  df_lite %>%
    # create a variable to denote the one-on-one correspondence 
    # between observations based on subject id, layer, and location
    unite(col = "set", c("subject_id", "f.layer", "location"), remove = FALSE) %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "CNO")), 
               y = get(var), 
               group = set, color = sex)) +
    geom_point(aes(shape = p15_cno, lwd = 0.5, fill = f.layer), size = 1.5, alpha = 1) +
    scale_x_discrete(labels = c("-", "+")) +
    scale_shape_manual(values = c("neg" = 1, "CNO" = 19)) +
    geom_line(linewidth = 0.5, alpha = 1) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~f.layer, labeller = labeller(f.layer =layer.label),
               strip.position = "bottom") +
    labs(x = "", y = "", title = title) +
    #scale_color_manual(values = color_sex)+
    theme_classic() +
    theme(legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size  = "10"),
          panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          #strip.background = element_rect(color = "black", linewidth = 1),
          strip.placement = ("outside"),
          strip.text = element_text(size = 10),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2),
          axis.ticks.x = element_blank())+
    scale_color_nejm()-> p
  return(p)
}


# merge layer plot row data-------------------------------------------------------------
#' ref. [Grouping on the x axis in ggplot2](https://stackoverflow.com/questions/44247239/grouping-on-the-x-axis-in-ggplot2)

# var = "silhs_mean_before_stat"
plot_raw_facet_sex = function(var, title = NULL, data = df_lite,
                    y_breaks = c(1:6 * 0.2),
                    y_limits = c(0, 1),
                    legend_position = c(0.3, 0.8)) {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  sex.label = c("Female","Male")
  names(sex.label)= c("F","M")
  #color_layer = MetBrewer::met.brewer("Hiroshige", 2)
  df_lite %>%
    # create a variable to denote the one-on-one correspondence 
    # between observations based on subject id, layer, and location
    unite(col = "set", c("subject_id", "layer", "location"), remove = FALSE) %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "CNO")), 
               y = get(var), 
               group = set, color = f.layer)) +
    geom_point(aes(shape = p15_cno, lwd = 0.5, fill = sex), size = 1.5, alpha = 1) +
    scale_x_discrete(labels = c("-", "+")) +
    scale_shape_manual(values = c("neg" = 1, "CNO" = 19)) +
    geom_line(linewidth = 0.5, alpha = 1) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~sex, labeller = labeller(sex = sex.label),
               strip.position = "bottom") +
    labs(x = "", y = "", title = title) +
    #scale_color_manual(values = color_layer) +
    theme_classic() +
    theme(legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size  = "10"),
          panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          #strip.background = element_rect(color = "black", linewidth = 1),
          strip.placement = ("outside"),
          strip.text = element_text(size = 10),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2),
          axis.ticks.x = element_blank())-> p
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

# x = "area"

report_main = function(x, var_name = x) {
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  mod_boot <<- mod_boot
  
  emm_list = lapply(cv_list[c(1:2, 4)], emm_var_main)
  
  anova = anova_list_select[[x]]
  
  table1 = kable(anova,caption = "ANOVA table before bootstrapping", digits = 3)
  print(table1)
  
  table2 = kable(emm_list, digits = 3)
  print(table2)
}


# report all interaction effect -----------------------------------------------
# x = "duration"
report_all_interaction = function(x, var_int1 = "p15_cno", var_int2 = "sex", var_int3 = "f.layer") {
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  f = paste(var_int1, "*", var_int2, "*", var_int3)
  emm = emmeans(mod_boot, eval(bquote(~ .(f))))
  emm_contr = contrast(emm, "pairwise")
  emm_param = as.data.frame(model_parameters(
    emm_contr, centrality = c("median", "mean"), test = "pd"))
  emm_param$pval = pd_to_p(emm_param$pd)
  
  table = kable(emm_param, caption = "complete comparisons", digits = 3)
  print(table)
}

report_all_interaction_simple = function(x, var_int1 = "p15_cno", var_int2 = "sex", var_int3 = "f.layer") {
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  f = paste(var_int1, "*", var_int2, "*", var_int3)
  emm = emmeans(mod_boot, eval(bquote(~ .(f))))
  emm_contr = contrast(emm, method = "consec")
  emm_param = as.data.frame(model_parameters(
    emm_contr, centrality = c("median", "mean"), test = "pd"))
  emm_param$pval = pd_to_p(emm_param$pd)
  
  table = kable(emm_param, caption = "complete comparisons-simple version", digits = 3)
  print(table)
}


# report interaction effect -----------------------------------------------
# x = "duration"
report_interaction = function(x, var_int1 = "p15_cno", var_int2 = "sex") {
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
  
  table = kable(emm_param, caption = "merge one cv-report others", digits = 3)
  print(table)
}


