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
source(here("src", "Develop_Function_subnetworkID_20230219.R"))


# import data
df = read_xlsx(
  here("data", "SummaryDevAll_Kmean_CovM_Corr_in_outV2_20231204.xlsx"),
  sheet = "ForMultilevel")

# check column names
df = clean_names(df)
names(df)

# list of response variables
rv_list = c("in_r_mean_cluster_1","in_r_mean_cluster_2","in_r_mean_cluster_3","in_r_mean_cluster_4","in_r_mean_cluster_5",
           "in_r_mean_cluster_6", "in_r_mean_cluster_7", "in_r_mean_cluster_8", "in_r_mean_cluster_9","in_r_mean_cluster_10",
           "in_r_mean_cluster_11", "in_r_mean_cluster_12","in_r_mean_cluster_13", "in_r_mean_cluster_14","in_r_mean_cluster_15", 
           "out_r_mean_cluster_1","out_r_mean_cluster_2","out_r_mean_cluster_3" ,"out_r_mean_cluster_4","out_r_mean_cluster_5",
           "out_r_mean_cluster_6","out_r_mean_cluster_7","out_r_mean_cluster_8","out_r_mean_cluster_9","out_r_mean_cluster_10",
           "out_r_mean_cluster_11","out_r_mean_cluster_12","out_r_mean_cluster_13","out_r_mean_cluster_14","out_r_mean_cluster_15",
           "zf_in_r_mean_cluster_1","zf_in_r_mean_cluster_2","zf_in_r_mean_cluster_3","zf_in_r_mean_cluster_4","zf_in_r_mean_cluster_5",
           "zf_in_r_mean_cluster_6","zf_in_r_mean_cluster_7" ,"zf_in_r_mean_cluster_8","zf_in_r_mean_cluster_9","zf_in_r_mean_cluster_10",
           "zf_in_r_mean_cluster_11" ,"zf_in_r_mean_cluster_12","zf_in_r_mean_cluster_13","zf_in_r_mean_cluster_14","zf_in_r_mean_cluster_15",
           "zf_out_r_mean_cluster_1","zf_out_r_mean_cluster_2","zf_out_r_mean_cluster_3","zf_out_r_mean_cluster_4","zf_out_r_mean_cluster_5",
           "zf_out_r_mean_cluster_6","zf_out_r_mean_cluster_7","zf_out_r_mean_cluster_8","zf_out_r_mean_cluster_9","zf_out_r_mean_cluster_10",
           "zf_out_r_mean_cluster_11","zf_out_r_mean_cluster_12","zf_out_r_mean_cluster_13","zf_out_r_mean_cluster_14","zf_out_r_mean_cluster_15")
rv_list = tolower(rv_list)

cv_list= c("sex","f.age","f.layer","location","n_cls")

# check if data is numeric
is.numeric(df$location) 

# check unique values of location
table(df$location)
#' very few 4, should remove in the analysis
df = filter(df, location != 4)


# create a temp version including rv_list and cv_list
df_temp = df %>%
  filter(location != 4) %>%
  mutate (f.age = factor(age), f.layer = factor(layer))%>%
  select(subject_id, all_of(cv_list), all_of(rv_list))

# overview of the variables
skimr::skim(df_temp, c(all_of(cv_list), all_of(rv_list)))


##the function to take wide_list and convert to long list
wideTolong = function(input_dataframe, wide_list){
  output_long_form = input_dataframe %>%
    select (subject_id, all_of(cv_list), all_of(wide_list))%>%      
    gather(condition , r_value, all_of(wide_list))
  return (output_long_form)
}

#convert wide to long form for in_r value
rv_list_in_r = rv_list[1:15]
df_lite = wideTolong(df_temp, rv_list_in_r)

names(df_lite)
names(df_lite)[8] <- "in_r"
names(df_lite)

#add index if need
#df_lite_in_r$index = seq.int(nrow(df_lite_in_r))

rv_list_out_r = rv_list[16:30]
output_out_r = wideTolong(df_temp, rv_list_out_r)
df_lite$out_r = paste(output_out_r$r_value)
names(df_lite)

is.numeric(df_lite$out_r)
df_lite$out_r =  as.numeric(df_lite$out_r)

is.numeric(df_lite$out_r)

rv_list_zf_in_r = rv_list[31:45]
output_zf_in_r = wideTolong(df_temp, rv_list_zf_in_r)
df_lite$zf_in_r = paste(output_zf_in_r$r_value)
names(df_lite)
is.numeric(df_lite$zf_in_r)
df_lite$zf_in_r = as.numeric(df_lite$zf_in_r)

rv_list_zf_out_r = rv_list[46:60]
output_zf_out_r = wideTolong(df_temp, rv_list_zf_out_r)
df_lite$zf_out_r = paste(output_zf_out_r$r_value)
names(df_lite)
df_lite$zf_out_r = as.numeric(df_lite$zf_out_r )

# split the character to get the subnetwork ID
temp_str = str_split_fixed(df_lite$condition, "_",5)
temp_str = data.frame(temp_str)
names(temp_str)

df_lite$subnetwork_id = paste(temp_str$X5)
names(df_lite)


is.numeric(df_lite$subnetwork_id)

#re-define rv_list and cv_list
rv_list = c("in_r", "out_r", "zf_in_r", "zf_out_r")
cv_list
cv_list[5] = "subnetwork_id"


# check if all outcome variables of interest are included in the data
summary(rv_list %in% names(df_lite))
summary(cv_list %in% names(df_lite))


#remove rows without value

df_lite = df_lite%>%
   filter(!is.na(in_r))

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

tab_model(mod_list_select[1:4])


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
save.image(file = here("output","2023-12-06_by_clusterV2_SummaryDevAll_Kmean_CovM_Corr_in_outV2_20231204.RData"))

