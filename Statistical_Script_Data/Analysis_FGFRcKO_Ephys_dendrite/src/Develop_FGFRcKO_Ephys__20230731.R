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
source(here("src", "Develop_FGFRcKO_Ephys_Function_20230731.R"))

# import data
df = read_xlsx(
  here("data", "FGFR-tri-Ncre_sEPSC_sIPSC_all_age_20230731.xlsx"),
  sheet = "Multilevel")

#
# check column names
df = clean_names(df)
names(df)

# list of response variables
rv_list = c("s_epsc_amp", "s_epsc_freq","s_ipsc_amp","s_ipsc_freq",
            "e_i_ratio","synaptic_drive_ratio_e_i")
rv_list = tolower(rv_list)

cv_list= c("group","section_id","sex","f_age","pipette_resistance_mo")
cv_list = tolower(cv_list)

# check if all outcome variables of interest are included in the data
summary(rv_list %in% names(df))
summary(cv_list %in% names(df))

# overview of the variables
skimr::skim(df, c(all_of(cv_list), all_of(rv_list)))


# create a lite version and convert the list of response variables to numeric
df_lite = df %>%
  select(subject_id, all_of(cv_list), all_of(rv_list))%>%
  rename ("f.age" = "f_age")

df_lite$pipette_resistance_mo = as.numeric(as.character(df_lite$pipette_resistance_mo))
df_lite["group"][df_lite ["group"] == "Neg"] <- "Ctrl"
df_lite["group"][df_lite ["group"] == "+"] <- "KO"

# density estimates of response variables
lapply(rv_list, function(x) with(df_lite, densityPlot(get(x), xlab = x)))

# check distribution of variable as need
summary(df_lite$pipette_resistance_mo)
summary(df_lite$s_epsc_amp)


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
save.image(file = here("output", "2023-07-31_FGFR-tri-Ncre_sEPSC_sIPSC_all_age_20230731.RData"))

