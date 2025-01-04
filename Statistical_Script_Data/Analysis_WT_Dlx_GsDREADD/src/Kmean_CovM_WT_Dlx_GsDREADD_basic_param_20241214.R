#' ---
#' title: "P15 euronal subnetwork assembly"
#' author: "Pei-Ying Chen"
#' date: "`r format(Sys.time(), '%B %d, %Y')`"
#' output: 
#'  pdf_document:
#'    latex_engine: xelatex
#' ---

# warning=FALSE
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
source(here("src", "WT_Dlx_GsDREADD_Function_20230211_caseboot.R"))

# import data
df = read_xlsx(
  here("data", "Kmean_CovM_SummaryDlxGsDREADD_CorrDistanceArea.xlsx"),
  sheet = "Sheet1")

# check column names
df = clean_names(df)
# names(df)

df <- df %>% rename(p15_cno = p15cno)

# list of response variables
rv_list = c("event_freq", "iei", "amp", "duration", "area", "synchrony_freq", 
            "synchrony_peak")

cv_list= c("sex", "layer", "location","p15_cno")

# check if all outcome variables of interest are included in the data
summary(rv_list %in% names(df))
summary(cv_list %in% names(df))

# overview of the variables
skimr::skim(df, c(all_of(cv_list), all_of(rv_list)))

# check unique values of animal_id
as.matrix(table(df$subject_id))

# create a lite version and convert the list of response variables to numeric
df_lite = df %>%
  select(subject_id, p15_cno, all_of(cv_list), all_of(rv_list)) %>%
  mutate(f.layer = factor(layer))
  
cv_list[2] = "f.layer"


# density estimates of response variables
# can see the plot of each response variable
lapply(rv_list, function(x) with(df_lite, densityPlot(get(x), xlab = x)))


# check distribution of RV
summary(df_lite$event_freq)
summary(df_lite$duration)

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
cl = makeCluster(16)
registerDoParallel(cl)
# default: b1 = 625, b2 = 16 --> B = 10000
case_boot_list = lapply(mod_list_lmer, case_bootstrap, b1 = 625, b2 = 16)
stopCluster(cl)
toc()

# save bootstrap output ---------------------------------------------------

save.image(file = here("output","Kmean_CovM_SummaryWT_Dlx_GsDREADD_basicPara.RData"))

