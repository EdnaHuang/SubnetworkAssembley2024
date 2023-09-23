#+ warning=FALSE
freshr::freshr()
library(readxl)
library(here)
library(janitor)
library(dplyr)
library(ggplot2)
library(car)
library(sjPlot)
library(tictoc)
library(stringr)
theme_set(theme_bw())
source(here("src", "FGFR_cKO_Ephys_Function_20230729.R"))

# import data
df = read_xlsx(
  here("data", "FGFR-tri-Ncre_sEPSC_sIPSC_P11-P13_20230729.xlsx"),
  sheet = "Multilevel")

# check column names
df = clean_names(df)
names(df)

# list of response variables
rv_list = c("pipette_resistance_mo",
            "s_epsc_amp", "s_epsc_freq","s_ipsc_amp","s_ipsc_freq",
            "e_i_ratio","synaptic_drive_ratio_e_i")
            
cv_list= c("genotype_cre","section_id","cell_id","sex",
           "pipette_resistance_mo")

# check if all outcome variables of interest are included in the data
summary(rv_list %in% names(df))
summary(cv_list %in% names(df))

# overview of the variables
skimr::skim(df, c(all_of(cv_list), all_of(rv_list)))

# check unique values of animal_id
table(df$subject_id)

# create a lite version and convert the list of response variables to numeric
df_lite = df %>%
  select(recording_date,id, subject_id,all_of(cv_list), all_of(rv_list))

df_lite["genotype_cre"][df_lite ["genotype_cre"] == "Neg"] <- "Ctrl"
df_lite["genotype_cre"][df_lite ["genotype_cre"] == "+"] <- "KO"


# density estimates of response variables
lapply(rv_list, function(x) with(df_lite, densityPlot(get(x), xlab = x)))


# apply to all outcome variables of interest
mod_list = lapply(rv_list, function(x) lmer_two_way(x))
names(mod_list) = rv_list
tab_model(mod_list)

# check which variables throw the warning message
sapply(rv_list, function(x) mod_list[[x]]@optinfo$conv$lme4$messages)


#find the best fit
mod_list_select = lapply(rv_list, model_select)
names(mod_list_select) = rv_list


# check which variables throw the warning message
sapply(rv_list, function(x) mod_list_select[[x]]@optinfo$conv$lme4$messages)


#' Error in mod_list_select[[x]]@optinfo : 
#no applicable method for `@` applied to an object of class "lm"

tab_model(mod_list_select)


# Wald test for linear mixed-effects models 
# Type III Wald F tests with Kenward-Roger degrees of freedom
anova_list_select = lapply(mod_list_select, Anova, type = 3, 
                           test.statistic = "F")
anova_list_select



# check if the ANOVA result match the regression table
tab_model(mod_list_select)

# After determining the correct model for fitting, employ bootstrap----------
mod_class = sapply(mod_list_select, function(x) class(x)[1])
mod_list_lm = mod_list_select[which(mod_class == "lm")]
mod_list_lmer = mod_list_select[which(mod_class == "lmerMod")]


# lm bootstrap
# default: R = 10000
lm_boot_list = lapply(mod_list_lm, lm_bootstrap, R = 10000)

summary(lm_boot_list[[1]])


# case bootstrap ---------------------------------------------------------------

Sys.time()
tic("case bootstrap")
set.seed(47408)
cl = makeCluster(10)
registerDoParallel(cl)
# default: b1 = 625, b2 = 16 --> B = 10000
case_boot_list = lapply(mod_list, case_bootstrap, b1 = 1000, b2 = 10)
stopCluster(cl)
toc()


# check bootstrap results -------------------------------------------------
summary(case_boot_list[[1]])

# save bootstrap output ---------------------------------------------------
save.image(file = here("output", "2023-07-29_FGFR-tri-Ncre_sEPSC_sIPSC_P11-P13_20230729.RData"))
#load(here("output", "2023-07-29_FGFR-tri-Ncre_sEPSC_sIPSC_P11-P13_20230729.RData"))
