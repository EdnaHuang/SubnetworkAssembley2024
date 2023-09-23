
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
library(ggsci)
library(dplyr)
library(emmeans)
library(parameters)
library(bayestestR)
library(ggplot2)
library(knitr)

theme_set(theme_classic())
source(here("src", "VGAT_Syt2_Function_20230319.R"))

# import data
df = read_xlsx(
  here("data", "P15_VGAT_SYT2_20230318.xlsx"),
  sheet = "Somatic_input_ForMultilevel")

# check column names
df = clean_names(df)
names(df)

# list of response variables
rv_list = c("soma_volum_mm3", 
            "vgat_punctate_number_on_soma", 
            "syt2_punctate_number_on_soma",
            "vgat_and_syt2_punctate_number_on_soma",
            "syt2_neg_vgat_punctate_number_on_soma",
            "vgat_punctate_number_on_soma_normalized",
            "syt2_punctate_number_on_soma_normalized",
            "vgat_and_syt2_punctate_number_on_soma_normalized",
            "syt2_neg_vgat_punctate_number_on_soma_normalized")

cv_list= c("image_number","sex","layer")


# check if all outcome variables of interest are included in the data
summary(rv_list %in% names(df))
summary(cv_list %in% names(df))

# check unique values of animal_id
table(df$animal_id)

# overview of the variables
skimr::skim(df, c(all_of(cv_list), all_of(rv_list)))

# create a lite version and convert the list of response variables to numeric
df_lite = df %>%
  select(animal_id, file_name, all_of(cv_list), all_of(rv_list))


# density estimates of response variables
lapply(rv_list, function(x) with(df_lite, densityPlot(get(x), xlab = x)))

# check distribution of RV if needed
summary(df_lite$soma_volum_mm3)
summary(df_lite$vgat_punctate_number_on_soma)
summary(df_lite$syt2_punctate_number_on_soma)
summary(df_lite$vgat_and_syt2_punctate_number_on_soma)


# apply to all outcome variables of interest
mod_list = lapply(rv_list, function(x) lmer_two_way(x))
names(mod_list) = rv_list

# check which variables throw the warning message
sapply(rv_list, function(x) mod_list[[x]]@optinfo$conv$lme4$messages)
tab_model(mod_list)

#find the correct model: lmer or lm
mod_list_select = lapply(rv_list, model_select)
names(mod_list_select) = rv_list

# check if the ANOVA result match the regression table
tab_model(mod_list_select)


# Wald test for linear mixed-effects models 
# Type III Wald F tests with Kenward-Roger degrees of freedom
anova_list_select= lapply(mod_list_select, Anova, type = 3, test.statistic = "F")
names(anova_list_select) = rv_list
anova_list_select


# After determining the correct model for fitting, employ bootstrap----------
mod_class = sapply(mod_list_select, function(x) class(x)[1])
mod_list_lm = mod_list_select[which(mod_class == "lm")]
mod_list_lmer = mod_list_select[which(mod_class == "lmerMod")]


# lm bootstrap
# default: R = 10000
lm_boot_list = lapply(mod_list_lm, lm_bootstrap, R = 10000)




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
save.image(file = here("output", "2023-08-04_P15_VGAT_SYT2_20230318_somatic_input.RData"))
#load(here("output", "2023-08-04_P15_VGAT_SYT2_20230318_somatic_input.RData"))
