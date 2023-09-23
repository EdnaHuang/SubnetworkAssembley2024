library(lme4)
library(lmeresampler)
library(foreach)
library(doParallel)
library(pbkrtest)
library(simpleboot)

#define MLM w/ two-way interaction --------------------------------------------
# x = rv_list
lmer_two_way = function(x, data = df_lite) {
  mod = lmer(paste(x, "~ f.layer * sex + (1 | subject_id) + (1 | subject_id:location)"), data = data, REML = TRUE)
  return(mod)
}



# model selection ---------------------------------------------------------

# x = rv_list[5]
# data = df_lite

model_select = function(x, data = df_lite) {
  
  mod = lmer(paste(x, "~ sex * f.layer + (1 | subject_id) + (1 | subject_id:location)"), data = data, REML = TRUE)
  
  if (isSingular(mod) == TRUE) {
    mod = lmer(paste(x, "~ sex * f.layer + (1 | subject_id)"), 
                data = data, REML = TRUE)
  }
  
  if (isSingular(mod) == TRUE) {
    f = as.formula(paste(x, "~ sex * f.layer"))
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

case_bootstrap = function(mod, b1 = 1000 , b2 = 10) {
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

# residual bootstrap ----------------------------------------------------------

residual_bootstrap = function(mod, b1 = 625 , b2 = 16) {
  tryCatch(
    {
      if (length(attributes(mod@flist)$names) == 1) {
        case_boot = foreach(B = rep(b1, b2), 
                            .combine = combine_lmeresamp,
                            .packages = c("lmeresampler", "lme4")) %dopar% {
                              bootstrap(mod, .f = fixef, 
                                        type = "residual", B = B)
                            }
      } else {
        case_boot = foreach(B = rep(b1, b2), 
                            .combine = combine_lmeresamp,
                            .packages = c("lmeresampler", "lme4")) %dopar% {
                              bootstrap(mod, .f = fixef, 
                                        type = "residual", B = B)
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

