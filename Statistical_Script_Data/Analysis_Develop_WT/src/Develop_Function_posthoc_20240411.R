library(pbkrtest)
library(lme4)
library(lmeresampler)

##extract the replicates and assign attributes and class -------

boot_model_extract = function(x) {
  # extract lmeresamp objects
  mod = case_boot_list[[x]]
  #mod = residual_boot_list[[x]]
 
   # save replicates
  mod_rep = mod$replicates
  
  # keep complete cases
  mod_rep = mod_rep[complete.cases(mod_rep), ]
  
  # assign attributes and class to the replicates
  attr(mod_rep, "original_model") = mod$model
  class(mod_rep) = append(class(mod_rep), "bootstrap_model")
  class(mod_rep) = append(class(mod_rep), "see_bootstrap_model")
  
  # return replicates
  return(mod_rep)
}

##extract the replicates and assign attributes and class -------
 # boot model extract raw function returns two data frames for plotting figure: 
 # (1) boot sample and (2) original data

boot_model_extract_raw = function(x) {
  # extract lmeresamp objects
  mod = case_boot_list[[x]]
 
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

##plot raw data and bootstrapping result  ------------------


plot_raw_boot_facet_layer_3age = function(x, var_name = x, title = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1), 
                                     legend_position = c(0.4, 0.85)) {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)[[1]]
  # extract the estimation data
  emms_all = emmip(mod_boot, sex ~ f.age |f.layer, CIs = TRUE,plotit=FALSE)
  emms = emms_all %>%
    filter (f.age != 13)%>%
    filter(f.age !=18)
  
  # extract data
  mod_data_all = boot_model_extract_raw(x)[[2]]
  mod_data = mod_data_all %>%
    filter (f.age != 13)%>%
    filter(f.age !=18)
  
  
  mod_data$yvar = mod_data[, x]
  mod_data$xvar = mod_data$f.age
  mod_data$tvar = mod_data$sex
  
  emm_plot = ggplot(emms, aes(x = xvar, y = yvar, color = sex)) +
    geom_pointrange(aes(ymin = LCL, ymax = UCL), size = 1,
                    position = position_dodge (width =0.5),
                    shape = 16, fatten = 1, alpha = 1, linewidth = 1) +
    
    geom_beeswarm(data = mod_data, dodge.width = 0.5,alpha = 0.7,
                  shape = 1, size = 1) +
    scale_x_discrete(limits=c("11","15", "21"), labels = c("P11","P15" ,"P21"))+
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ f.layer, labeller = labeller(f.layer =layer.label),
               strip.position = "bottom") +
    ggtitle(title) +
    labs(x = "", y= "", color = "Sex")+
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
          axis.ticks.x = element_blank())
  
  print(emm_plot)
}


plot_raw_boot_facet_layer_2age = function(x, var_name = x, title = x,
                                          y_breaks = c(1:6 * 0.2), 
                                          y_limits = c(0, 1), 
                                          legend_position = c(0.4, 0.85)) {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)[[1]]
  # extract the estimation data
  emms_all = emmip(mod_boot, sex ~ f.age |f.layer, CIs = TRUE,plotit=FALSE)
  emms = emms_all %>%
    filter (f.age != 13)%>%
    filter (f.age != 15)%>%
    filter(f.age !=18)
  
  # extract data
  mod_data_all = boot_model_extract_raw(x)[[2]]
  mod_data = mod_data_all %>%
    filter (f.age != 13)%>%
    filter (f.age != 15)%>%
    filter(f.age !=18)
  
  
  mod_data$yvar = mod_data[, x]
  mod_data$xvar = mod_data$f.age
  mod_data$tvar = mod_data$sex
  
  emm_plot = ggplot(emms, aes(x = xvar, y = yvar, color = sex)) +
    geom_pointrange(aes(ymin = LCL, ymax = UCL), size =1,
                    position = position_dodge (width =0.5),
                    shape = 16, fatten = 1, alpha = 1, linewidth = 1) +
    
    geom_beeswarm(data = mod_data, dodge.width = 0.5,alpha = 0.7,
                  shape = 1, size = 1) +
    scale_x_discrete(limits=c("11", "21"), labels = c("P11", "P21"))+
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ f.layer, labeller = labeller(f.layer =layer.label),
               strip.position = "bottom") +
    ggtitle(title) +
    labs(x = "", y= "", color = "Sex")+
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
          axis.ticks.x = element_blank())
  
  print(emm_plot)
}

##plot raw data and bootstrapping result -Facet by layer ------------------

plot_raw_boot_facet_layer = function(x, var_name = x, title = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1), 
                                     legend_position = c(0.4, 0.85)) {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)[[1]]
  
  # extract data
  mod_data = boot_model_extract_raw(x)[[2]]
  
  mod_data$yvar = mod_data[, x]
  mod_data$xvar = mod_data$f.age
  mod_data$tvar = mod_data$sex
  emm_plot = emmip(mod_boot, sex ~ f.age |f.layer, CIs = TRUE, 
                   CIarg = list(lwd = 1, alpha = 1), 
                   dodge = 0.5,linewidth = 0.5,fatten = 1) + 
    geom_beeswarm(data = mod_data, dodge.width = 0.5,alpha = 0.7,
                  shape = 1, size = 1) +
    
    scale_x_discrete(limits=c("11","13","15", "18","21"), labels = c("P11","P13","P15", "P18","P21"))+
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ f.layer, labeller = labeller(f.layer =layer.label),
               strip.position = "bottom") +
    ggtitle(title) +
    labs(x = "", y= "", color = "Sex")+
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
          axis.ticks.x = element_blank())
  
  print(emm_plot)
}

plot_raw_boot_merge_sex= function(x, var_name = x, title = x,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1), 
                                     legend_position = c(0.4, 0.85)) {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)[[1]]
  
  # extract data
  mod_data = boot_model_extract_raw(x)[[2]]
  
  mod_data$yvar = mod_data[, x]
  mod_data$xvar = mod_data$f.age
  mod_data$tvar = mod_data$f.layer
  
  emm_plot = emmip(mod_boot, f.layer ~ f.age, CIs = TRUE, 
                   CIarg = list(lwd = 2, alpha = 1), 
                   dodge = 0.5) +
  
    geom_beeswarm(data = mod_data, dodge.width = 0.5,alpha = 0.7,
                  shape = 1, size = 1) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    
    ggtitle(title) +
    labs(x = "", y= "", color = "Sex")+
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
          axis.ticks.x = element_blank()) +
          scale_color_jama()
  
  print(emm_plot)
}

##plot raw data and bootstrapping result -Facet by sex ------------------

plot_raw_boot_facet_sex = function(x, var_name = x,title = x,
                                   y_breaks = c(1:6 * 0.2), 
                                   y_limits = c(0, 1), 
                                   legend_position = c(0.4, 0.85)){
  # define label
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)[[1]]
  
  # extract data
  mod_data = boot_model_extract_raw(x)[[2]]
  
  mod_data$yvar = mod_data[, x]
  mod_data$xvar = mod_data$f.age
  mod_data$tvar = mod_data$f.layer 
  emm_plot = emmip(mod_boot, f.layer ~ f.age | sex, CIs = TRUE, 
                   CIarg = list(lwd = 1, alpha = 1), 
                   dodge = 0.5, linewidth = 0.5) +
    geom_beeswarm(data = mod_data, dodge.width = 0.5, alpha = 0.7,
                  shape = 1, size = 1)+
    scale_y_continuous(breaks = y_breaks, limits = y_limits)+
    facet_wrap(~ sex, labeller = labeller(sex =sex.label),
               strip.position = "bottom") +
    ggtitle(title) +
    labs(x = "",y = "", color = "Layer") +
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
          axis.ticks.x = element_blank())+
    scale_color_jama()
  
  print(emm_plot)
}


## custom contrasts --------------------------------------------------------

pw_emm = emmeans:::revpairwise.emmc(1:5)
as.matrix(names(pw_emm))
pw_emm_contr = pw_emm[, c(1, 3, 6, 10, 2, 9)]
names(pw_emm_contr) = c("age13 - age11", "age15 - age13", "age18 - age15", 
                        "age21 - age18", "age15 - age11", "age21 - age15")


## report main of sex, layer, age to compare ANOVA table before boots -------------

report_main_effect = function(x, var_name = x) {
  
  # extract replicates
  mod_boot = boot_model_extract(x)
  
  # anova of selected model-before bootstrapping
  anova = anova_list_select[[x]]
 
  # emmeans: age
  emm_age = emmeans(mod_boot, "f.age")
  #emm_contr_age = contrast(emm_age, method = pw_emm_contr)
  emm_contr_age = contrast(emm_age, method = "pairwise")
  emm_param_age = as.data.frame(model_parameters(
    emm_contr_age, centrality = c("median", "mean"), test = "pd"))
  emm_param_age$pval = pd_to_p(emm_param_age$pd)
  
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
  
  table2 = kable(emm_param_age, caption ="Main effect after bootstrapping", digits = 3)
  print(table2)
  
  table3 = kable(emm_param_sex ,digits = 3)
  print(table3)
  
  table4 = kable(emm_param_layer, digits = 3)
  print(table4)
  
}



## report age-sex-layer -----the most comprehensive comparison--------------
report_inter_age_sex_layer = function(x, var_name = x) {
  
  # extract replicates
  mod_boot = boot_model_extract(x)
  
  # simple comparisons
  emm = emmeans(mod_boot, ~ sex* f.layer*f.age)
  
  # age
  emm_contr_age = contrast(emm, method = "pairwise", simple = "f.age")
  #emm_contr_age = contrast(emm, method = pw_emm_contr, simple = "f.age")
  emm_param_age = as.data.frame(model_parameters(
    emm_contr_age, centrality = c("median", "mean"), test = "pd"))
  emm_param_age$pval = pd_to_p(emm_param_age$pd)
  emm_param_age$Median.1 = NULL
  
  # sex
  emm_contr_sex = contrast(emm, method = "consec", simple = "sex")
  emm_param_sex = as.data.frame(model_parameters(
    emm_contr_sex, centrality = c("median", "mean"), test = "pd"))
  emm_param_sex$pval = pd_to_p(emm_param_sex$pd)
  emm_param_sex$Median.1 = NULL
  
  # layer
  emm_contr_layer = contrast(emm, method = "consec", simple = "f.layer")
  emm_param_layer = as.data.frame(model_parameters(
    emm_contr_layer, centrality = c("median", "mean"), test = "pd"))
  emm_param_layer$pval = pd_to_p(emm_param_layer$pd)
  emm_param_layer$Median.1 = NULL
  
  # print formatted tables in sequence
  table1 = kable (emm_param_age,
                 caption = "Post-hoc comparison with bootstrapping output", digits = 3)
  print(table1)
  table2 = kable (emm_param_sex,digits = 3)
  print(table2)
  table3 = kable (emm_param_layer,digits = 3)
  print(table3)
}


