#'
#' script to automatically calculate seroconversion rates (SCRs) by EU for trachoma in 1-5 year olds
#' 
#' 
#' written by zach reynolds, 
#' based off of original code by Ben Arnold's group at UCSF
#'
#'      
#'        Script Overview
#'
# 1. load packages/estimate scr glm function  
# 2. analysis prep
# 3. calculate SCRs/seroprev by age (1-5)
# 4. fit age-seroprevalence curves (1-5)
# 5. write data to an excel file
# 6. save analysis objects for plotting

#------------------------------------------------------------------------------#

#load packages ----

#------------------------------------------------------------------------------#

library(MASS)
library(readxl)
library(writexl)
library(tidyverse)
library(binom)
library(foreach)
library(sandwich)
library(lmtest)
library(mgcv)

#NOTE: the MASS Package will overwrite the dplyr select function if you already have dplyr loaded  
#       
#     if that happens:

#detach("package:tidyr", unload = TRUE)
detach("package:dplyr", unload = TRUE)

library(dplyr)
library(tidyr)

#------------------------------------------------------------------------------#

#load estimate_scr_glm() function ----

#------------------------------------------------------------------------------#

#   formalized data cleaning code is available for LFA, MBA, and demo datasets

#------------------------------------------------------------------------------#

#if you're running this code by itself here's where to set your wd 
#which is the file path to where your data files are 
#(use the same pattern as i have shown here with your user id)
#remove the "#" before it and the next portion 

#setwd("C:/Users/userid/OneDrive - CDC/trachoma/analysis workflow/test env")

#and heres where to import your data
#replace with your data file name
#working_dataset <- read_excel("test_data_2.xlsx")

#------------------------------------------------------------------------------#

wd_mba <- read_excel("mba_2012.xlsx")
wd_elisa <- read_excel("elisa_2022.xlsx")


#uncomment this to run the analysis on the mba 2012 data
filtered_data <- wd_mba %>%
  mutate(res_cluster = "9999")

#uncomment this to run the analysis on the elisa 2022 data
#filtered_data <- wd_elisa %>%
#  mutate(res_cluster = "9999")
#------------------------------------------------------------------------------#

#   analysis prep + analysis for 1–5 year olds below                            

#------------------------------------------------------------------------------#

#variables that need to be included in your final cleaned dataframe


#sample_id     character                     unique id for each sample
#res_eu	       factor or character/integer	 inique id for evaluation unit (EU); used for stratification/faceting
#res_cluster	 factor or integer             cluster id for sampling unit; used for robust SEs
#age           numeric	                     age in years (integer); must be positive
#seropos	     binary numeric (1/0)	         seropositivity indicator: 1 = positive, 0 = negative



#prep data for 1–5 year olds ----

sero_df_1to5 <- filtered_data %>% 
  filter(!is.na(age) & !is.na(seropos)) %>% #filter out where age is NA and serostatus is NA
  filter(age %in% 1:5) %>%                  #set age from 1 to 5
  mutate(res_eu = as.factor(res_eu)) #  %>%
#can add a filter step here to filter out , or specifically look at one EU
#filter(res_eu == "1") #eg, only look at eu 1
#filter(res_eu != "3") #eg, exclude eu 3

#calculate seroprevalence by EU  ----

#creates our 'seroprev_by_eu' dataframe, where you can view seropositivty by EU
#(in ages 1-5)

seroprev_by_eu_1to5 <- sero_df_1to5 %>%
  group_by(res_eu) %>%
  summarise(
    age = "1 to 5 year olds - combined",
    n = n(),
    seroprev = mean(seropos, na.rm = TRUE) * 100,
    lower_ci = binom.confint(x = round(seroprev / 100 * n), n = n, method = "wilson")$lower * 100,
    upper_ci = binom.confint(x = round(seroprev / 100 * n), n = n, method = "wilson")$upper * 100,
    .groups = "drop"
  )

print(seroprev_by_eu_1to5)


#calculate seroprevalence by age (1–5) ----

#creates our 'seroprev_by_age_1to5' dataframe, where you can view seropositivty by EU and age

seroprev_by_age_1to5 <- sero_df_1to5 %>%
  group_by(res_eu, age) %>% #add ",age" to this if you want seroprev by age and eu
  summarise(
    n = n(),
    seroprev = mean(seropos, na.rm = TRUE) * 100,
    lower_ci = binom.confint(x = round(seroprev / 100 * n), n = n, method = "wilson")$lower * 100,
    upper_ci = binom.confint(x = round(seroprev / 100 * n), n = n, method = "wilson")$upper * 100,
    .groups = "drop"
  )


print(seroprev_by_age_1to5)
#------------------------------------------------------------------------------#

#   estimate scr by EU for 1–5 ----

#------------------------------------------------------------------------------#


#creates our 'scr_summary_1to5' dataframe, 
#this is the meat of the script and calls our estimate_scr_glm function from above

#this block loops over each EU/district (using the filtered 1 to 5 year old dataframe from above)
#it calculcates the SCR for each EU and the confidence intervals and scales them to 100 (for easier interpretability)

#e.g. an SCR of 3.6 per 100 person years suggests that 3.6 kids out of 100 children are seroconverting per year
#also when using the GLM we report confidence intervals and when you use a bayesian MCMC report credible intervals


#model choice: 

#This estimate_scr_glm() function is a SIR model that assumes no seroreversion 

#Ben Arnold's group compared this to an reversible catalytic model - 
#a variation of a SIS model and the rank-order scr estimates did not change

#Using the SIR (assuming no seroreversion) model is appropriate 

#they also compared this GLM model to the bayesian MCMC approach that they used in the manuscript
# and found that SCR estimates aligned closely and the simplified glm approach was 
# appropriate for future use


#threshold interpretation: 

# Our thresholds vary - see the Kamau/Arnold 2024 paper for more information 

# in areas with typical epidemiology: 
# A threshold of >2.2 suggests additional monitoring should considered
# a threshold of >4.5 sugests addditional moniting is required

# in areas with unusual epidemiology, those thresholds become 1.6 and 3.8 respectively 


scr_summary_1to5 <- foreach(curr_eu = unique(sero_df_1to5$res_eu), .combine = rbind) %do% {
  
  #filters the data for the current EU
  df_temp <- sero_df_1to5 %>%
    filter(res_eu == curr_eu)
  
  #count number of unique clusters in the current eu/district
  n_clusters_temp <- df_temp %>%
    filter(!is.na(res_cluster)) %>%
    summarise(
      n_clusters = n_distinct(res_cluster)
    ) %>%
    pull(n_clusters)
  
  #choose variance method automatically
  variance_used_temp <- if_else(
    n_clusters_temp > 1,
    "robust",
    "ML"
  )
  
  #estimate the SCR using the glm function
  scr_out <- estimate_scr_glm(
    serop = df_temp$seropos,
    agey = df_temp$age,
    id = df_temp$res_cluster,
    variance = variance_used_temp
  )
  
  #store the results in a tibble
  tibble(
    res_eu = curr_eu,
    n = nrow(df_temp),
    n_clusters = n_clusters_temp,
    variance_used = variance_used_temp,
    scr = scr_out$scr,
    scr_min95 = scr_out$scr_min95,
    scr_max95 = scr_out$scr_max95
  )
  
} %>%
  
  #scale SCRs and create new columns for excel output 
  mutate(
    across(c(scr, scr_min95, scr_max95), ~ round(. * 100, 2), .names = "{.col}_scaled"),
    
    #create a formatted label for graphs, e.g.: "3.45 (1.23–5.67)"
    #rounding to 1 decimal place 
    scr_label = sprintf("%.1f (%.1f–%.1f)", scr_scaled, scr_min95_scaled, scr_max95_scaled)
    
  ) %>%
  
  #reorder for excel table
  dplyr::select(res_eu,  n, n_clusters, variance_used,
                scr_scaled, scr_min95_scaled, scr_max95_scaled, scr_label, scr, scr_min95, scr_max95)


print(scr_summary_1to5)
#------------------------------------------------------------------------------#

#fit age-seroprevalence curves for 1–5 ----

#------------------------------------------------------------------------------#
#this adapts ben arnold's code for estimating the age-seroprevalence curves:

#We estimated seroprevalence by age using semiparametric cubic splines in a generalized additive model
#to allow for potential non-linear relationships with age, specifying binomial errors for seroprevalence, 
#and random effects for clusters to account for repeated observations

#"lines represent mean seroprevalence by age estimated using semiparametric cubic splines"


#NOTE: generating the curve is seperate from how we estimate the SCR, however, buried in one of bens papers is a 
#is a comparison between a GAM estimating an SCR at each age and it was essentially equiavlent to a fixed SCR across all ages 


age_sero_curve_results_1to5 <- foreach(curr_eu = unique(sero_df_1to5$res_eu), .combine = rbind) %do% {
  
  #filters the data for the current EU
  df_temp <- sero_df_1to5 %>%
    filter(res_eu == curr_eu)
  
  #count total numbe rof positive children district
  n_pos_total <- sum(df_temp$seropos, na.rm = TRUE)
  
  #fits a generalized additive model with:
  # - binomial family for binary seropositivity outcome (1/0)
  # - logit link (log-odds of being seropositive)
  # - cubic regression spline ("cr") with 4 basis functions (k = 4)
  # - CIs generated via parametric boostrapping
  
  #- note: confidence interval calculations here are not absolutely necessary
  #        removing them is not unreasonable 
  #        Also adding a random effect per cluster to account for correlation is not a bad idea
  #        but it may drive low prevalence district age-seroprevalence curves flat 
  
  gam_fit <- mgcv::gam(
    seropos ~ s(age, k = 4, bs = "cr"),
    data   = df_temp,
    family = binomial(link = "logit"),
    method = "REML"
  )
  
  #distinct observed age grid for prediction OR
  #age_grid <- df_temp %>%
  #  dplyr::select(age) %>%
  #  distinct() %>%
  #  arrange(age)
  
  #use this for much smoother lines
  age_grid <- seq(min(df_temp$age), max(df_temp$age), by = 0.1)
  
  #build the GAM design matrix at those ages 
  Xp <- predict(
    gam_fit,
    newdata = data.frame(age = age_grid),
    #if you decide to add random intercept by cluster  
    
 #  newdata = data.frame(age = age_grid, res_cluster = df_temp$res_cluster[1]),  
    type = "lpmatrix" 
 #  ,exclude = "s(res_cluster)"
  )
  
  #extract coefficient estimates (mean & covariance)
  beta_hat <- coef(gam_fit)                             #fitted gam coefficents
  V        <- vcov(gam_fit, unconditional = TRUE)       #covariance matrix
  
  #draw coefficients and build CI for the mean curve
  #this is the parametric bootstrapping
  set.seed(555)
  B <- 1000    #                                        #number of bootstrap draws
  beta_draws <- MASS::mvrnorm(B, mu = beta_hat, Sigma = V)
  
  #calculcate predicted values 
  eta_hat    <- as.vector(Xp %*% beta_hat)              #linear predictor using mean coefficients
  eta_draws  <- Xp %*% t(beta_draws)                    #predictions using bootstrap sampled coefficients
  
  #if our districts have less than 5 positives we can choose to remove the CIs around the 
  #seroprevalence curve due to them being extremely wide
  if (n_pos_total < 5) {
    nd <- tibble(
      age      = age_grid,
      seroprev = plogis(eta_hat) * 100,
      lower_ci = NA_real_,
      upper_ci = NA_real_
    )
  } else {
    
    #convert our predictions from log-odds to probabilities which equals percent seropositive
    #and compute 95% Ci bands from our bootstrap distribution 
    nd <- tibble(
      age      = age_grid,
      seroprev = plogis(eta_hat) * 100,
      lower_ci = apply(eta_draws, 1, function(z) quantile(plogis(z), 0.025)) * 100,
      upper_ci = apply(eta_draws, 1, function(z) quantile(plogis(z), 0.975)) * 100
    )
  }
  
  #return results 
  tibble(
    res_eu   = curr_eu,
    age      = nd$age,
    seroprev = nd$seroprev,
    lower_ci = nd$lower_ci,
    upper_ci = nd$upper_ci
  )
}


#------------------------------------------------------------------------------#

#join scr summary to smoothed estimates for 1–5 ----

age_sero_curve_results_1to5 <- age_sero_curve_results_1to5 %>%
  left_join(scr_summary_1to5, by = "res_eu")

#------------------------------------------------------------------------------#





#                      STOP                      HERE






#------------------------------------------------------------------------------#



#--------------------------------------------------------------------------#

#save MBA 1-5 outputs ----

#one for the mba dataset run 
#and one for the ELISA run (below)
#------------------------------------------------------------------------------#
#pull scr summary, seroprev summary by EU and by age and EU
mba_output_list_1to5 <- list(
  "SCR Estimates (1-5)" = scr_summary_1to5,
  "Seroprev by EU (1-5)" = seroprev_by_eu_1to5,
  "Seroprev by Age (1-5)" = seroprev_by_age_1to5
)

mba_1to5 <- list(
  dataset_label = "MBA, 2012",
  sero_df = sero_df_1to5,
  seroprev_by_eu = seroprev_by_eu_1to5,
  seroprev_by_age = seroprev_by_age_1to5,
  scr_summary = scr_summary_1to5,
  age_sero_curve_results = age_sero_curve_results_1to5
)

#write_xlsx(mba_output_list_1to5,"mba_scr_and_seroprev_results_1to5.xlsx")
#------------------------------------------------------------------------------#

#save ELISA 1-5 outputs ----

#------------------------------------------------------------------------------#
#pull scr summary, seroprev summary by EU and by age and EU
elisa_output_list_1to5 <- list(
  "SCR Estimates (1-5)" = scr_summary_1to5,
  "Seroprev by EU (1-5)" = seroprev_by_eu_1to5,
  "Seroprev by Age (1-5)" = seroprev_by_age_1to5
)


elisa_1to5 <- list(
  dataset_label = "ELISA, 2022",
  sero_df = sero_df_1to5,
  seroprev_by_eu = seroprev_by_eu_1to5,
  seroprev_by_age = seroprev_by_age_1to5,
  scr_summary = scr_summary_1to5,
  age_sero_curve_results = age_sero_curve_results_1to5
)

#write_xlsx(elisa_output_list_1to5,"elisa_scr_and_seroprev_results_1to5.xlsx")
#------------------------------------------------------------------------------#

#save analysis objects for plotting script ----

#------------------------------------------------------------------------------#

analysis_output_1to5 <- list(
  mba_1to5 = mba_1to5,
  elisa_1to5 = elisa_1to5
)

#save directly to working directory
saveRDS(
  analysis_output_1to5,
  "analysis_output_1to5.rds"
)
