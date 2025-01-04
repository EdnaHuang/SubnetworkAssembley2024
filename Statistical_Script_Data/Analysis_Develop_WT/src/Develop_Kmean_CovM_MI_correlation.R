freshr::freshr()
library(readxl)
library(here)
library(janitor)
library(dplyr)
library(tidyr)
library (tidyverse)
library(car)
library(sjPlot)
library(tictoc)

library(ggplot2)
theme_set(theme_bw())
source(here("src", "Develop_Function_20230219.R"))

# import data
df = read_xlsx(
  here("data", "Kmean_CovM_MI_Corr.xlsx"),
  sheet = "Sheet1")

# check column names
df = clean_names(df)
column_names_all = names(df)
print(column_names_all)

#define cv list
cv_list = column_names_all[c(4:6,8:11)]
cv_list

#define rv list
rv_list = column_names_all[c(33:38)]
rv_list 


# check if data is numeric
is.numeric(df$location) 

# check unique values of location
table(df$location)
#' very few 4, should remove in the analysis
df = filter(df, location != 4)

# create a temp version including rv_list and cv_list
df_lite = df %>%
  filter(location != 4) %>%
  mutate (f.age = factor(age), f.layer = factor(layer))%>%
  select(subject_id, all_of(cv_list), all_of(rv_list))

# update cv_list
cv_list[3] = 'f.layer'
cv_list [4] = 'f.age'
cv_list

# overview of the variables
simple_test <-skimr::skim(df_lite, c(all_of(cv_list), all_of(rv_list)))

# check if all outcome variables of interest are included in the data
summary(rv_list %in% names(df_lite))
summary(cv_list %in% names(df_lite))

#remove rows without value
df_lite = df_lite%>%
   filter(!is.na(in_r_mean))

# density estimates of response variables
lapply(rv_list, function(x) with(df_lite, densityPlot(get(x), xlab = x)))

# check distribution of variable as need
lapply(rv_list, function(x) with(df_lite, is.numeric(get(x))))


#Three-Way analysis
# apply to all outcome variables of interest
mod_list = lapply(rv_list, function(x) lmer_three_way(x))
names(mod_list) = rv_list
tab_model(mod_list)

# check which variables throw the warning message
sapply(rv_list, function(x) mod_list[[x]]@optinfo$conv$lme4$messages)

# model selection ---------------------------------------------------------
mod_list_select = lapply(rv_list, model_select)
names(mod_list_select) = rv_list
tab_model(mod_list_select)

# Wald test for linear mixed-effects models 
# Type III Wald F tests with Kenward-Roger degrees of freedom
anova_list_select = lapply(mod_list_select, Anova, type = 3, 
                           test.statistic = "F")
anova_list_select


# After determining the correct model for fitting, employ bootstrap---------
mod_class = sapply(mod_list_select, function(x) class(x)[1])
mod_list_lm = mod_list_select[which(mod_class == "lm")]
mod_list_lmer = mod_list_select[which(mod_class == "lmerMod")]


# lm bootstrap
# default: R = 10000
lm_boot_list = lapply(mod_list_lm, lm_bootstrap, R = 10000)

# lmer bootstrap
Sys.time()
tic("case bootstrap")
set.seed(47408)
cl = makeCluster(10)
registerDoParallel(cl)
# default: b1 = 625, b2 = 16 --> B = 10000
case_boot_list = lapply(mod_list_lmer, case_bootstrap, b1 = 1000, b2 = 10)
stopCluster(cl)
toc()


# save bootstrap output ---------------------------------------------------
save.image(file = here("output","SummaryDevAll_Kmean_CovM_MI_Corr.RData"))

