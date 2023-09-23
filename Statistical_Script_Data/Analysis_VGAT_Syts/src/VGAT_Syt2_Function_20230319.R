library(lme4)
library(pbkrtest)
library(simpleboot)
library(lmeresampler)
library(lme4)
library(foreach)
library(doParallel)
library(dplyr)
library(emmeans)
library(parameters)
library(bayestestR)
library(ggplot2)
library(knitr)

#define MLM w/ two-way interaction --------------------------------------------
# x = rv_list
lmer_two_way = function(x, data = df_lite) {
  mod = lmer(paste(x, "~ layer * sex + (1 | animal_id) + (1 | animal_id:image_number)"), data = data, REML = TRUE)
  return(mod)
}


#solve the isSingular issue by removing layer or change to lm model
model_select = function(x, data = df_lite) {
  
  mod = lmer(paste(x, "~ layer * sex + (1 | animal_id) + (1 | animal_id:image_number)"), data = data, REML = TRUE)
  
  if (isSingular(mod) == TRUE) {
    mod = lmer(paste(x, "~ layer * sex + (1 | animal_id)"), data = data, REML = TRUE)
  }
  
  if (isSingular(mod) == TRUE) {
    f = as.formula(paste(x, "~ layer * sex"))
    mod = eval(bquote(lm(.(f), data = data)))
  }
  return(mod)
}

# simple bootstrap --------------------------------------------------------

lm_bootstrap = function(mod, R = 10000) {
  lm_boot = lm.boot(mod, R = R, rows = TRUE)
  return(lm_boot)
}



# case bootstrap ----------------------------------------------------------

case_bootstrap = function(mod, b1 = 625 , b2 = 16) {
  tryCatch(
    {
      if (length(attributes(mod@flist)$names) == 1) {
        case_boot = foreach(B = rep(b1, b2), 
                            .combine = combine_lmeresamp,
                            .packages = c("lmeresampler", "lme4")) %dopar% {
                              bootstrap(mod, .f = fixef, 
                                        resample = c(TRUE, FALSE), 
                                        type = "case", B = B)
                            }
      } else {
        case_boot = foreach(B = rep(b1, b2), 
                            .combine = combine_lmeresamp,
                            .packages = c("lmeresampler", "lme4")) %dopar% {
                              bootstrap(mod, .f = fixef, 
                                        resample = c(TRUE, FALSE, FALSE), 
                                        type = "case", B = B)
                            }
      }
      return(case_boot)
    }, 
    error = function(e) {
      message("An error occurred.")
      print(e)
    }
  )
}




## plot raw data------------------------------------------
#var = "vgat_and_syt2_punctate_number_on_soma"

plot_raw = function(var, data = df_lite, title = var,
                                y_breaks = c(1:6 * 0.2), 
                                y_limits = c(0, 1), 
                                legend_position = c(0.4, 0.85)) {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("L2/3", "L4")
  sex.label = c("Female","Male")
  names(sex.label) = c("F", "M")
  
  df_lite %>%
    select(all_of(var), animal_id, sex) %>%
    ggplot(aes(x = sex, y = get(var),
               group = sex, color = sex)) +
    ggbeeswarm::geom_beeswarm(shape = 19, size = 1,
                              dodge.width = 0.5, alpha = 0.3) +
    #stat_summary(fun = "mean", linewidth = 1,
    #             position = position_dodge(width = 0.5)) +
    
    stat_summary(fun= mean, geom = "point",position = position_dodge(width = 0.5))+
        stat_summary(fun.data = "mean_cl_boot", size = 0.5,linewidth = 1,
                 position = position_dodge(width = 0.5)) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    #facet_wrap(~ layer, labeller = labeller(layer = layer.label)) +
    labs(x = "", y = "", title = title) + 
    theme(strip.text = element_text(size = 10), 
          strip.background = element_rect(color = "black", linewidth = 1), 
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position, legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank()) -> p
  return(p)
}

