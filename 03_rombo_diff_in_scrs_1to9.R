#'
#' script written by zach reynolds for comparing two (independently estimated) SCRs
#' 
#------------------------------------------------------------------------------#

#load packages ----
library(tidyverse)

#load estimate_scr_glm() function
#------------------------------------------------------------------------------#

#this estimates the scr for the mba (2012 study) data
res_mba <- estimate_scr_glm(
  serop = wd_mba$seropos,
  agey  = wd_mba$age,
  id    = wd_mba$sample_id,
  variance = "ML"
)

#this estimates the scr for the elisa (2022 study) data
res_elisa <- estimate_scr_glm(
  serop = wd_elisa$seropos,
  agey  = wd_elisa$age,
  id    = wd_elisa$sample_id,
  variance = "ML"
)

#------------------------------------------------------------------------------#

#extract log SCR estimates and variances ----

#------------------------------------------------------------------------------#
#notes:
#for a single-rate catalytic model with no covariates:

#log(-log(1 - p_i)) = log(a_i) + log(lambda)

#where:
#p_i     = probability individual i is seropositive
#a_i     = age of individual i
#lambda  = seroconversion rate

#the glm uses log(age_i) as an offset
#the model intercept estimates log(lambda)

#following Neal's notation:
#beta = log(lambda)

#therefore:
#beta_mba   = log(SCR_MBA)
#beta_elisa = log(SCR_ELISA)

#there are a few assumptions we are making here:
# 1. this comparison assumes that the two SCR estimates are independent 
#     which in this case, is reasonable because surveys were ~10 years apart
#
# 2. assumes the esimate log scrs are approx. normally distributed
#     this is a large-sample wald approximation
#     and is more reasonable when the sample size/# of pos is not small
#     and why we compared across 1-9s instead of 1-5s
#------------------------------------------------------------------------------#

#we have the log of the scr for both estimators
#beta = log(scr)
beta_mba   <- as.numeric(log(res_mba$scr))
beta_elisa <- as.numeric(log(res_elisa$scr))

#and their variances
V_mba   <- as.numeric(res_mba$logscr_se)^2
V_elisa <- as.numeric(res_elisa$logscr_se)^2

#assuming the log scr estimates are approximately normally distributed 
#this is calculating the difference in the two SCRs on the log scale
#positive values mean elisa scr > mba scr and vice versa
d  <- beta_elisa - beta_mba


#then we can sum the variances and 
#estimate the standard error of the difference between the two log SCRs
#assuming covariance between surveys is zero (independence)
se <- sqrt(V_elisa + V_mba)

#how big is the observed difference compared to what we expect
#calculcate z score

#essentially            the difference we saw (d)
# divided by        ---------------------------------     = z
#                    the difference due to chance (se)
z <- d / se


#if these two scrs were the same how often would we see a difference this large by chance
p <- 2 * (1 - pnorm(abs(z)))

#95% CIs on log scale
ci_log <- c(d - 1.96 * se, d + 1.96 * se)
ci_log

#then we exponentiate our difference between two SCRs and exponentiate it for our rate ratio
rr <- exp(d)
ci_rr  <- exp(ci_log)

#------------------------------------------------------------------------------#

#output in a labelled list 
list(
  SCR_MBA   = as.numeric(res_mba$scr),
  SCR_ELISA = as.numeric(res_elisa$scr),
  rate_ratio_ELISA_vs_MBA = rr,
  rr_lower = ci_rr[1],
  rr_upper = ci_rr[2],
  p_value = p
)
