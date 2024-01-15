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



##merge sex then plot raw data color by sex-----------
merge_sex_plot_raw_color_sex = function(var, data = df_lite, title = var,
                                        y_breaks = c(1:6 * 0.2), 
                                        y_limits = c(0, 1), 
                                        legend_position = "none") {
  
  df_lite %>%
    select(all_of(var), subject_id, sex, f.age, group) %>%
    ggplot(aes(x = f.age, y = get(var), group = group, color = sex)) +
    ggbeeswarm::geom_beeswarm(shape = 19, size = 1,
                              dodge.width = 0.5, alpha = 0.3) +
    stat_summary(fun.data = "mean_cl_boot", size = 0.5,linewidth = 1,
                 position = position_dodge(width = 0.5),fatten = 1.5) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    #facet_wrap(~ layer, labeller = labeller(layer = layer.label)) +
    labs(x = "", y = "", title = title) + 
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
    scale_color_nejm() ->p
  
  return(p)
  
}



## report main of group, sex, age to compare ANOVA table before boots -------------
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
  
  # emmeans: sex
  emm_sex = emmeans(mod_boot, "sex")
  emm_contr_sex = contrast(emm_sex, method = "pairwise")
  emm_param_sex = as.data.frame(model_parameters(
    emm_contr_sex, centrality = c("median", "mean"), test = "pd"))
  emm_param_sex$pval = pd_to_p(emm_param_sex$pd)
  
  table1 = kable(anova,caption = "ANOVA table before bootstrapping", digits = 3)
  print(table1)
  
  table2 = kable(emm_param_group, caption ="Main effect after bootstrapping", digits = 3)
  print(table2)
  
  table3 = kable(emm_param_age ,digits = 3)
  print(table3)
  
  table4 = kable(emm_param_sex, digits = 3)
  print(table4)
  
  
}


## report age-group-sex-----the most comprehensive comparison--------------
report_inter_age_group_sex = function(x, var_name = x) {
  
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  # simple comparisons
  emm = emmeans(mod_boot, ~ group* sex *f.age)
  
  # age
  emm_contr_age = contrast(emm, method = "pairwise", simple = "f.age")
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
  
  # sex
  emm_contr_sex = contrast(emm, method = "consec", simple = "sex")
  emm_param_sex = as.data.frame(model_parameters(
    emm_contr_sex, centrality = c("median", "mean"), test = "pd"))
  emm_param_sex$pval = pd_to_p(emm_param_sex$pd)
  emm_param_sex$Median.1 = NULL
  
  # print formatted tables in sequence
  table1 = kable (emm_param_group,
                  caption = "Post-hoc comparison with bootstrapping output", digits = 3)
  print(table1)
  
  table2 = kable (emm_param_age,digits = 3)
  print(table2)
  
  table3 = kable (emm_param_sex,digits = 3)
  print(table3)
  
}



## Merge sex, then report interaction between age and group---------

merge_sex_report_age_int_group = function(x, var_name = x) {
  
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
  table1 = kable (emm_param_group,
                  caption = "Merge-sex-Post-hoc comparison with bootstrapping output", digits = 3)
  print(table1)
  
  table2 = kable (emm_param_age,digits = 3)
  print(table2)
}





## report group inter layer ---------
report_group_int_sex = function(x, var_name = x) {
  
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  # simple comparisons
  emm = emmeans(mod_boot, ~ group* sex)
  
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






