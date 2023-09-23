freshr::freshr()
library(readxl)
library(here)
library(janitor)
library(dplyr)
library(ggplot2)
library(car)
library(sjPlot)
library(tictoc)

theme_set(theme_bw())
source(here("src", "Develop_Function_20230219.R"))

#+ import data
df = read_xlsx(
  here("data", "SummaryDevAll_Kmean_CovM_withCorr_Dist_20230113.xlsx"),
  sheet = "ForMultilevel")


# check column names
df = clean_names(df)
names(df)

# list of response variables
rv_list = c("silhs_mean_before_Stat","n_cls","No_assemblies", "S_assemblies", 
            "M_assemblies", "Cells_not_in_assembly", "Cells_in_one_assembly", 
            "Cells_in_many_assembly", "R_allpair_mean", "Cluster_R_all_mean", 
            "R_ZF_allpair_mean", "Cluster_R_ZF_all_mean", "Neuron_Centroid", 
            "PW_Neuron_Neuron", "Centroid_Centroid")
rv_list = tolower(rv_list)
rv_list
cv_list= c("sex", "age", "layer", "location")

# check if all outcome variables of interest are included in the data
summary(rv_list %in% names(df))
summary(cv_list %in% names(df))

# overview of the variables
skimr::skim(df, c(all_of(cv_list), all_of(rv_list)))

# check unique values of location
table(df$location)
#' very few 4, should remove in the analysis


# create a lite version and convert the list of response variables to numeric
df_lite = df %>%
  filter(location != 4) %>%
  select(subject_id, all_of(cv_list), all_of(rv_list)) %>%
  mutate(centroid_centroid.num = as.numeric(centroid_centroid), 
         f.age = factor(age), 
         f.layer = factor(layer))
cv_list[2] = "f.age"

summary(df_lite$centroid_centroid.num)
df_lite$centroid_centroid[is.na(df_lite$centroid_centroid.num)]

#check if NA still in data
summary(df_lite$centroid_centroid == "NA")
#' okay, so NA does exist in the original variable

df_lite$centroid_centroid = df_lite$centroid_centroid.num

# density estimates of response variables
lapply(rv_list, function(x) with(df_lite, densityPlot(get(x), xlab = x)))

# check distribution of variable as need
summary(df_lite$n_cls)

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

tab_model(mod_list_select[1:3])
tab_model(mod_list_select[4:6])
tab_model(mod_list_select[7:10])


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
save.image(file = here("output","2023-02-19_SummaryDevAll_Kmean_CovM_withCorr_Dist_20230113.RData"))
