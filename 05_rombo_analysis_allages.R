#####################################################################
#####################################################################
##                                                                 ##
## CODE FOR FITTING MODELS TO CROSS-SECTIONAL ANTIBODY TITRE DATA  ##
##                                                                 ##
## Please feel free to share modify the code as you see fit        ##   
## (but please maintain appropriate accreditation)                 ##
##                                                                 ##   
## Michael White                                                   ##
## Institut Pasteur                                                ##
## michael.white@pasteur.fr                                        ##
##                                                                 ##
#####################################################################
#####################################################################

#this script was originally written by Michael White, 
#and has since been edited by zach reynolds

#------------------------------------------------------------------------------#

#load packages ----

#------------------------------------------------------------------------------#

library(MASS)
library(compiler)
library(binom)
library(coda)
library(plyr)

#------------------------------------------------------------------------------#

#set run id ----

#------------------------------------------------------------------------------#

#change this for each independent MCMC run, e.g. "V1", "V2", "V3", etc. 
#each district should have its own v1, v2, v3 etc

RUN_ID <- "V1"
#RUN_ID <- "V2"
#RUN_ID <- "V3"

###############################################
###############################################
##          ##                               ##
##   ####   ##  ####    ####  ######  ####   ##
##  ##  ##  ##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ##  ##  ##  ## ######   ##   ######  ##
##  ##  ##  ##  ## ##  ##  ##   ##   ##  ##  ##
##   ####   ##  ####   ##  ##   ##   ##  ##  ##
##          ##                               ##
###############################################
###############################################

#------------------------------------------------------------------------------#

#0.1 - read in data ----

#------------------------------------------------------------------------------#
#set working directory to the folder containing the input CSV files
#setwd("C:/Users/.../OneDrive - CDC/Trachoma/tanzania_mpwapwa_rombo/scenario 2")


#read in antibody data

#column 1 must be age in years
#column 3 must be binary serostatus/seropositivity coded 1/0 (1= pos / 0= neg)

AB_data <- read.csv(file = "MBA2012allages.csv", header = TRUE)
#AB_data <- read.csv(file = "ELISAallages.csv", header = TRUE)

#make sure titre is numeric
#if AB_data$titre is used, the csv must also include a column named titre
AB_data$titre <- as.numeric(AB_data$titre)

#------------------------------------------------------------------------------#

#0.2 - prepare binned seroprevalence data for plotting ----

#------------------------------------------------------------------------------#

#create 5-year age bins from 0 to 60 years
age_bins     <- seq(from = 0, to = 60, by = 5)
age_bins_mid <- seq(from = 2.5, to = 57.5, by = 5)

#number of age bins used for the binned seroprevalence plot
N_bins <- length(age_bins) - 1

#set up matrix to store binomial point estimate and 95% CI for each age bin
SP_bins <- matrix(NA, nrow = N_bins, ncol = 3)
colnames(SP_bins) <- c("med", "low_95", "high_95")

#loop over age bins and calculate Wilson confidence intervals for seroprevalence
for (i in 1:N_bins) {
  index   <- which(AB_data[,1] > age_bins[i] & AB_data[,1] <= age_bins[i+1])
  temp_AB <- AB_data[index, 3]
  SP_bins[i, ] <- as.numeric(as.vector(
    binom.confint(sum(temp_AB), length(temp_AB), method = "wilson")[1, 4:6]
  ))
}

#------------------------------------------------------------------------------#

#0.3 - quick data plot ----

#------------------------------------------------------------------------------#

par(mfrow = c(1,1))
plot(x = age_bins_mid, y = SP_bins[,1],
     pch = 15, cex = 2,
     xlim = c(0, 65), ylim = c(0, 1),
     xlab = "age (years)", ylab = "Proportion seropositive",
     main = "Cross-sectional serological data")

for (i in 1:N_bins) {
  arrows(x0 = age_bins_mid[i], y0 = SP_bins[i,2],
         x1 = age_bins_mid[i], y1 = SP_bins[i,3],
         length = 0.03, angle = 90, code = 3, col = "black", lwd = 1)
}

###################################################
###################################################
##        ##                                     ##
##   ##   ##  #     #  ####  ####   ##### ##     ##
##  ###   ##  ##   ## ##  ## ## ##  ##    ##     ##
##   ##   ##  ####### ##  ## ##  ## ####  ##     ##
##   ##   ##  ## # ## ##  ## ## ##  ##    ##     ##
##  ####  ##  ##   ##  ####  ####   ##### #####  ##
##        ##                                     ##
###################################################
###################################################

#------------------------------------------------------------------------------#

#1.1 - define SIS catalytic model ----

#------------------------------------------------------------------------------#

#this function defines the reversible catalytic model used to estimate
#age specific seroprevalence under a change-point model

#lambda_0 = baseline seroconversion rate before the change point
#gamma    = multiplier that scales lambda_0 after the change point
#rho      = seroreversion rate
#time_c   = estimated time of change in transmission

#lambda_c is the post-change seroconversion rate:
#lambda_c = gamma * lambda_0

model_M2 <- function(a, par_M2) {
  lambda_0 <- par_M2[1]
  gamma    <- par_M2[2]
  rho      <- par_M2[3]
  time_c   <- par_M2[4]
  
  lambda_c <- gamma * lambda_0
  age_xx   <- a - time_c
  
  if (age_xx <= 0) {
    SP_prop <- (lambda_c / (lambda_c + rho)) * (1 - exp(-(lambda_c + rho) * a))
  }
  
  if ((age_xx > 0) && (age_xx < a)) {
    SP_prop <- lambda_c / (lambda_c + rho) +
      (((lambda_0 - lambda_c) * rho) / ((lambda_0 + rho) * (lambda_c + rho))) *
      exp(-(lambda_c + rho) * (a - age_xx)) -
      (lambda_0 / (lambda_0 + rho)) *
      exp(-(lambda_c + rho) * a) * exp(-(lambda_0 - lambda_c) * age_xx)
  }
  
  if (a <= age_xx) {
    SP_prop <- (lambda_0 / (lambda_0 + rho)) * (1 - exp(-(lambda_0 + rho) * a))
  }
  
  SP_prop
}

model_M2 <- cmpfun(model_M2, options = list(optimize = 3))


#------------------------------------------------------------------------------#

#1.2 - define likelihood ----

#------------------------------------------------------------------------------#

#this function calculates the binomial log-likelihood for a proposed parameter set
loglike_M2 <- function(par_M2) {
  #predict modeled seroprevalence for each observed age
  SP_model <- sapply(AB_data[,1], model_M2, par_M2 = par_M2)
  loglike  <- AB_data[,3] * log(SP_model) + (1 - AB_data[,3]) * log(1 - SP_model)
  sum(loglike)
}
loglike_M2 <- cmpfun(loglike_M2, options = list(optimize = 3))


#------------------------------------------------------------------------------#

#1.3 - define priors ----

#------------------------------------------------------------------------------#

LARGE <- 1e10     ## Large value for rejecting parameters with prior

#this function returns the log prior density for a proposed parameter set
#parameter values outside the allowed prior ranges are rejected with a very large penalty


prior_M2 <- function(par_M2) {
  lambda_0 <- par_M2[1]
  gamma    <- par_M2[2]
  rho      <- par_M2[3]
  time_c   <- par_M2[4]
  
#------------------------------------------------------------------------------#
  
#1.3.1 - lambda_0 prior options ----
  
#------------------------------------------------------------------------------#
  
#lambda_0:
#baseline/pre-change seroconversion rate

  #lambda_0 ~ U(0, 10)
  #allows up to 1000 seroconversions per 100 person-years
  #very broad prior; may implausibly high SCR values
  
  #if (lambda_0 > 0 && lambda_0 < 10) {
  #  prior_lambda_0 <- log(1/10)
  #} else {
  #  prior_lambda_0 <- -LARGE
  #}
  
  
  #lambda_0 ~ U(0, 1)
  #allows up to 100 seroconversions per 100 person-years
  #moderately broad prior
  
  #if (lambda_0 > 0 && lambda_0 < 1) {
  #  prior_lambda_0 <- log(1/1)
  #} else {
  #  prior_lambda_0 <- -LARGE
  #}
  
  
  #lambda_0 ~ U(0, 0.5)
  #allows up to 50 seroconversions per 100 person-years
  #this was used when broader lambda_0 priors resulted in poor mixing
  
  if (lambda_0 > 0 && lambda_0 < 0.5) {
    prior_lambda_0 <- log(1/0.5)
  } else {
    prior_lambda_0 <- -LARGE
  }
  
#------------------------------------------------------------------------------#
  
#1.3.2 - gamma prior ----
  
#------------------------------------------------------------------------------#
  
#gamma ~ U(0, 1)
  
  if (gamma > 0 && gamma < 1) {
    prior_gamma <- log(1/1)
  } else {
    prior_gamma <- -LARGE
  }
  
#original informative gamma prior from Pinsent et al. (i think)
#this prior is suggesting that post change is expected to be about 5.3% of pre-change transmission
#or approximately a 94.7% reduction

#this may be too informative if the goal is to estimate whether transmission changed
#rather than assume a large reduction a priori

  
#prior_gamma <- dnorm( x=gamma, mean=0.053, sd = 0.03, log=TRUE) #Pgp3 
  
  
#------------------------------------------------------------------------------#
  
#1.3.3 - rho prior options ----
  
#------------------------------------------------------------------------------#
  
#rho ~ U(0, 10)
#broad uniform prior on seroreversion
#allows up to 1000 seroreversions per 100 person-years
  
  if (rho > 0 && rho < 10) {
    prior_rho <- log(1/10)
  } else {
    prior_rho <- -LARGE
  }
  
  
#rho ~ normal(0.026, 0.003), 
#informative prior from pinsent et al
#centers rho around 2.6 per 100 person-years
  
  #if (rho > 0 && rho < 10) {
  #  prior_rho <- dnorm(x = rho, mean = 0.026, sd = 0.003, log = TRUE)  ## case for pgp3
  #} else {
  #  prior_rho <- -LARGE
  #}
  
#------------------------------------------------------------------------------#
  
#1.3.4 - time_c prior options ----
  
#------------------------------------------------------------------------------#
  
#time_c ~ U(0, 60)
  #allows the change point to occur 0 to 60 years before sampling
  if (time_c > 0 && time_c < 60) {
    prior_time_c <- log(1/60)
  } else {
    prior_time_c <- -LARGE
  }
  
  

#time_c ~ U(0, 90)
#allows the change point to occur 0 to 90 years before sampling
  
 #if (time_c > 0 && time_c < 90) {
 #  prior_time_c <- log(1/90)
 #} else {
 #  prior_time_c <- -LARGE
 #}

#------------------------------------------------------------------------------#
#sum log priors
  
  prior <- prior_lambda_0 + prior_gamma + prior_rho + prior_time_c
  
  prior
}

#compile prior function for faster repeated calls during MCMC
prior_M2 <- cmpfun(prior_M2, options = list(optimize = 3))


#################################################
#################################################
##          ##                                 ##
##   ####   ##  #     #  ####  #     #  ####   ##
##  ##  ##  ##  ##   ## ##  ## ##   ## ##  ##  ##
##     ##   ##  ####### ##     ####### ##      ##
##    ##    ##  ## # ## ##  ## ## # ## ##  ##  ##
##   #####  ##  ##   ##  ####  ##   ##  ####   ##
##          ##                                 ##
#################################################
#################################################

#------------------------------------------------------------------------------#

#2.1 - MCMC setup + sampling + diagnostics ----

#------------------------------------------------------------------------------#

#edit as needed
#set MCMC length and tuning/adaptation periods
N_mcmc       <- 10000
N_tune_start <- 500
N_tune_end   <- 5000
N_adapt      <- 6000

#scaled up
#N_mcmc       <- 40000
#N_tune_start <- 800
#N_tune_end   <- 8000
#N_adapt      <- 10000

#initialize proposal step size and acceptance counter
step_scale  <- 1
MCMC_accept <- 0

#------------------------------------------------------------------------------#

#2.2 - Robbins-munro step scaler ----

#------------------------------------------------------------------------------#
#Robbins-munro adaptation updates the proposal step size during early iterations
rm_scale <- function(step_scale, mc, log_prob) {
  dd <- exp(log_prob)
  if (dd < -30) dd <- 0
  dd <- min(dd, 1)
  rm_temp <- (dd - 0.23) / ((mc + 1) / (0.01 * N_adapt + 1))
  out <- step_scale * exp(rm_temp)
  out <- max(out, 0.05)
  out <- min(out, 5)
  out
}

#------------------------------------------------------------------------------#

#2.3 - prepare object for MCMC fitting ----

#------------------------------------------------------------------------------#

N_par <- 4

#store parameters plus log-likelihood and prior for every MCMC iteration
MCMC_par           <- matrix(NA, nrow = N_mcmc, ncol = N_par + 2)
colnames(MCMC_par) <- c("lambda_0", "gamma", "rho", "time_c", "loglike", "prior")

#------------------------------------------------------------------------------#

#2.4 - implement MCMC iterations ----

#------------------------------------------------------------------------------#

#starting values for lambda_0, gamma, rho, and time_c
par_MC   <- c(0.01, 0.1, 0.01, 10)     #(lambda_0, gamma, rho, time_c)
par_MCp1 <- rep(NA, N_par)

Sigma_MC <- diag((0.25 * par_MC)^2)    #initial guess of covariance of MVN proposal dist

prior_MC   <- prior_M2(par_MC)
loglike_MC <- loglike_M2(par_MC) + prior_MC     ## posterior log (ll + log prior)

#------------------------------------------------------------------------------#

#main Metropolis-Hastings loop
for (mc in 1:N_mcmc) {
  #-----------------------------------------------
  #propose new parameter
  par_MCp1 <- mvrnorm(n = 1, mu = par_MC, Sigma = step_scale * Sigma_MC)

  #-----------------------------------------------
  #prior check
  prior_MCp1 <- prior_M2(par_MCp1)
  
  if (prior_MCp1 > -0.5 * LARGE) {
    
    #-----------------------------------------------  
    ## MH accept/reject on posterior (ll + log prior)
    loglike_MCp1 <- loglike_M2(par_MCp1) + prior_MCp1
    log_prob     <- min(loglike_MCp1 - loglike_MC, 0)
    
    if (log(runif(1)) < log_prob) {
      par_MC     <- par_MCp1
      loglike_MC <- loglike_MCp1
      prior_MC   <- prior_MCp1
      MCMC_accept <- MCMC_accept + 1
    }
    
    #----------------------------------------------- 
    #RM scaling of step size
    if (mc < N_adapt) {
      step_scale <- rm_scale(step_scale, mc, log_prob)
    }
    
    #----------------------------------------------- 
    #adaptive tuning of covariance
    if ((mc > N_tune_start) && (mc < N_tune_end)) {
      Sigma_MC <- cov(MCMC_par[1:(mc-1), 1:N_par], use = "complete.obs")
    }
  }
  
  #----------------------------------------------- 
  #store output in a matrix
  MCMC_par[mc, 1:N_par] <- par_MC
  MCMC_par[mc, N_par+1] <- loglike_MC
  MCMC_par[mc, N_par+2] <- prior_MC
}

paste("Acceptance rate = ", round(100 * MCMC_accept / N_mcmc, 2), "%", sep = "")

#------------------------------------------------------------------------------#

#2.5 - examine MCMC chains ----

#------------------------------------------------------------------------------#
#examine chains for lambda_0, gamma, rho, and time_c
#with ess in title

par(mfrow = c(3,5))
par(mar = c(4,4,2,2))
for (j in 1:N_par) {
  plot(x = 1:N_mcmc, y = MCMC_par[, j],
       pch = 19, col = "grey", cex = 0.25,
       xlab = "MCMC iteration", ylab = colnames(MCMC_par)[j],
       main = paste("MCMC chain: ", colnames(MCMC_par)[j],
                    " (ESS = ", round(effectiveSize(MCMC_par[, j]), 0), ")", sep = ""))
}

#examine the mcmc log likelihood plot
par(mar = c(1,1,1,1))
plot(x = 1:N_mcmc, y = MCMC_par[, 5],
     pch = 19, col = "grey", cex = 0.25,
     ylim = quantile(MCMC_par[, 5], prob = c(0.01, 1)),
     xlab = "MCMC iteration", ylab = "log-likelihood",
     main = "log-likelihood")

#------------------------------------------------------------------------------#

#2.6 - examine posterior distributions ----

#------------------------------------------------------------------------------#
#drop the first 20% of samples as burn-in
MCMC_burn <- MCMC_par[floor(0.2 * nrow(MCMC_par)):(nrow(MCMC_par) - 1), ]

#examine posterior density plots for lambda_0, gamma, rho, and time_c
for (j in 1:N_par) {
  DEN   <- density(MCMC_burn[, j])
  QUANT <- quantile(MCMC_burn[, j], prob = c(0.025, 0.5, 0.975))
  
  plot(x = DEN$x, y = DEN$y, type = 'l',
       xlim = c(0, max(DEN$x)),
       xlab = colnames(MCMC_par)[j], ylab = "",
       main = paste("Posterior profile: ", colnames(MCMC_par)[j], sep = ""))
  
  low_index  <- which(DEN$x < QUANT[1])
  mid_index  <- intersect(which(DEN$x >= QUANT[1]), which(DEN$x <= QUANT[3]))
  high_index <- which(DEN$x > QUANT[3])
  
  polygon(x = c(DEN$x[low_index], rev(DEN$x[low_index])),
          y = c(rep(0, length(low_index)), rev(DEN$y[low_index])),
          col = "pink")
  
  polygon(x = c(DEN$x[mid_index], rev(DEN$x[mid_index])),
          y = c(rep(0, length(mid_index)), rev(DEN$y[mid_index])),
          col = "grey")
  
  polygon(x = c(DEN$x[high_index], rev(DEN$x[high_index])),
          y = c(rep(0, length(high_index)), rev(DEN$y[high_index])),
          col = "pink")
  
  points(x = rep(QUANT[2], 2), y = c(0, max(DEN$y)), type = 'l', lty = "dashed", lwd = 2)
}
plot.new()

#------------------------------------------------------------------------------#

#2.7 - examine autocorrelation ----

#------------------------------------------------------------------------------#

#examine autocorrelation plots for lambda_0, gamma, rho, and time_c
for (j in 1:N_par) {
  autocorr.plot(MCMC_par[, j], auto.layout = FALSE,
                main = paste("Auto-correlation: ", colnames(MCMC_par)[j]))
}
plot.new()

#ercentage of proposed MCMC steps that were accepted
#very low acceptance can indicate proposal jumps are too large
#very high acceptance can indicate proposal jumps are too small

paste("Acceptance rate = ", round(100 * MCMC_accept / N_mcmc, 2), "%", sep = "")


#############################################
#############################################
##          ##                             ##
##   ####   ##  ###### #####  ###  ######  ##
##  ##  ##  ##    ##   ##    ##      ##    ##
##     ##   ##    ##   ####   ###    ##    ##
##  ##  ##  ##    ##   ##       ##   ##    ##
##   ####   ##    ##   #####  ###    ##    ##
##          ##                             ##
#############################################
#############################################

#------------------------------------------------------------------------------#

#3.1 - calculate posterior median curve ----

#------------------------------------------------------------------------------#

#prediction grid for plotting the fitted seroprevalence curve
age_seq <- seq(from = 0, to = 60, by = 0.1)

#calculate the posterior median parameter set
par_median <- apply(X = MCMC_burn[, 1:N_par], MARGIN = 2, FUN = median)

M2_predict_median <- sapply(age_seq, model_M2, par_M2 = par_median)

#------------------------------------------------------------------------------#

#3.2 - calculate posterior predictive credible intervals ----

#------------------------------------------------------------------------------#

#sample posterior draws to generate the posterior predictive credible interval
N_sam   <- 500
sam_seq <- round(seq(from = 1, to = nrow(MCMC_burn), length = N_sam))

M2_predict <- matrix(NA, nrow = N_sam, ncol = length(age_seq))
for (k in 1:N_sam) {
  M2_predict[k, ] <- sapply(age_seq, model_M2, par_M2 = MCMC_burn[sam_seq[k], 1:N_par])
}

M2_quant <- matrix(NA, nrow = 3, ncol = length(age_seq))
for (j in 1:length(age_seq)) {
  M2_quant[, j] <- quantile(M2_predict[, j], prob = c(0.025, 0.5, 0.975))
}

#------------------------------------------------------------------------------#

#3.3 - plot and record current run ----

#------------------------------------------------------------------------------#

#this function plots the observed age seroprevalence points and overlays the
#posterior model-predicted seroprevalence curve
#
#the plot includes:
# - observed seroprevalence by age bin
# - wilson 95% confidence intervals for observed age-bin seroprevalence
# - posterior median model-predicted seroprevalence curve
# - posterior 95% credible interval around the model-predicted curve
# - posterior median estimate of the change point, time_c
#
#recordPlot() is used at the end so the plot can be saved/replayed later


#function for plotting and recording the current model fit
plot_current_panel <- function(panel_title = "ELISA, 2022 (all ages)",
                               draw_top_strip = FALSE,
                               top_strip_title = NULL) {
  par(mfrow = c(1,1),
      mar   = c(5, 5, 3, 2),
      oma   = c(0, 0, 2.2, 0),
      mgp   = c(2.2, 0.8, 0))
  
  #base plot with bin medians
  plot(x = age_bins_mid, y = SP_bins[,1],
       pch = 16, cex = 1.5,
       xlim = c(0, 60),
       ylim = c(0, 1),
       xlab = "Age (years)",
       ylab = "Proportion seropositive",
       cex.axis = 1.2, cex.lab = 1.3,
       xaxt = "n")
  
  axis(side = 1, at = seq(0, 60, by = 10), labels = seq(0, 60, by = 10))
  
  #add wilson 95% confidence intervals around observed age-bin seroprevalence
  for (i in seq_len(N_bins)) {
    arrows(x0 = age_bins_mid[i], y0 = SP_bins[i,2],
           x1 = age_bins_mid[i], y1 = SP_bins[i,3],
           length = 0.03, angle = 90, code = 3, col = "black", lwd = 1)
  }
  
  #add posterior 95% credible interval ribbon for the model-predicted curve
  polygon(x = c(age_seq, rev(age_seq)),
          y = c(M2_quant[1,], rev(M2_quant[3,])),
          col = rgb(154/256,170/256,235/256,0.5), border = NA)
  lines(age_seq, M2_quant[2,], lwd = 3, col = "blue")
  
  #posterior median and 95% credible interval for the change point
  tc_q <- quantile(MCMC_burn[, "time_c"], c(0.025, 0.5, 0.975))
  
  usr <- par("usr")
  
  #add vertical dashed line at posterior median change point
  segments(x0 = tc_q[2],
           y0 = usr[3] + 0.01*(usr[4]-usr[3]),
           x1 = tc_q[2],
           y1 = usr[4] - 0.02*(usr[4]-usr[3]),
           col = "grey60", lty = 2, lwd = 2)
  
  op <- par(xpd = NA)
  
  #place labels above the plot to distinguish post-change and pre-change periods
  y_lab <- usr[4] + 0.025*(usr[4]-usr[3])
  text(x = 6, y = y_lab, labels = "time 2 (post-change)", adj = c(0.0, 0), cex = 1)
  text(x = (tc_q[2] + usr[2]) / 2.0, y = y_lab, labels = "time 1 (pre-change)", adj = c(0.5, 0), cex = 1)
  par(op)
  
  mtext(panel_title, side = 3, line = 0.4, outer = TRUE, cex = 1.3, font = 2)
  
  #add change-point estimate text to the upper-left corner of the plot
  round1 <- function(x) sprintf("%.1f", x)
  
  #edit change-point label text here
  tc_expr <- bquote(
    bold("Change-point (tc): ") * .(round1(tc_q[2])) * "y" *
      " (" * .(round1(tc_q[1])) * "-" * .(round1(tc_q[3])) * "y)"
  )
  op <- par(xpd = NA)
  inset_xy <- c(0.02, 0.06)
  lg <- legend("topleft", inset = inset_xy, legend = list(tc_expr),
               bty = "n", plot = FALSE, cex = 1.0)
  rect(lg$rect$left, lg$rect$top - lg$rect$h,
       lg$rect$left + lg$rect$w, lg$rect$top,
       col = rgb(1,1,1,0.95), border = NA)

  legend("topleft", inset = inset_xy, legend = list(tc_expr),
         bty = "n", cex = 1.0)
  par(op)
  
  if (isTRUE(draw_top_strip) && !is.null(top_strip_title)) {
    op <- par(xpd = NA)
    x <- grconvertX(c(0, 1), from = "ndc", to = "user")
    y <- grconvertY(c(0.965, 1.00), from = "ndc", to = "user")
    rect(xleft = x[1], ybottom = y[1], xright = x[2], ytop = y[2],
         col = "white", border = NA)
    text(x = mean(x), y = mean(y), labels = top_strip_title, cex = 1.3, font = 2)
    par(op)
  }
  
  #record the plot so it can be replayed/saved later
  recordPlot()
}

#------------------------------------------------------------------------------#

#STOP HERE TO ADD TITLE (and make sure to change it between plots)

#this is what actually generates the plot in the viewer 

#------------------------------------------------------------------------------#
#mba
plot_obj <- plot_current_panel(panel_title = "MBA, 2012 (all ages)",
                               draw_top_strip = FALSE)

#elisa
#plot_obj <- plot_current_panel(panel_title = "ELISA, 2022 (all ages)",
#                               draw_top_strip = TRUE, top_strip_title = "ELISA, 2022 (all ages)")


######################################
######################################
##          ##                      ##
##  ##      ##  ####   ####  ####   ##
##  ## ##   ##  ## ##   ##  ##  ##  ##
##  ######  ##  ##  ##  ##  ##      ##
##     ##   ##  ## ##   ##  ##  ##  ##
##     ##   ##  ####   ####  ####   ## 
##          ##                      ##
######################################
######################################

#------------------------------------------------------------------------------#

#4.1.1 - calculate DIC using posterior mean ----

#------------------------------------------------------------------------------#
#DIC = Deviance Information Criterion
#used here to compare model fit across candidate models/prior choices
#lower DIC indicates better relative fit

#calculate the center of the posterior distribution using the posterior mean
#parameters are lambda_0, gamma, rho, and time_c
theta_bar   <- apply(X = MCMC_burn[,1:4], FUN = mean, MARGIN = 2)

#calculate deviance at the posterior mean
D_theta_bar <- -2 * loglike_M2(theta_bar)

#MCMC_burn[,5] - MCMC_burn[,6] is the log-likelihood
D_bar       <- -2 * mean(MCMC_burn[,5] - MCMC_burn[,6])

#effective number of parameters
pD          <- D_bar - D_theta_bar

#mean DIC estimate
DIC_mean         <- pD + D_bar

#------------------------------------------------------------------------------#

#4.1.2 - calculate DIC using posterior median ----

#------------------------------------------------------------------------------#

#calculate the center of the posterior distribution using the posterior median
#parameters are lambda_0, gamma, rho, and time_c
theta_bar   <- apply(X = MCMC_burn[,1:4], FUN = median, MARGIN = 2)
D_theta_bar <- -2 * loglike_M2(theta_bar)

#calculate the deviance at the posterior median parameter values
D_bar       <- -2 * median(MCMC_burn[,5] - MCMC_burn[,6])
pD_median   <- D_bar - D_theta_bar
DIC_median  <- pD_median + D_bar

#median-centered DIC
DIC_median

#------------------------------------------------------------------------------#

#4.2 - summarize posterior medians and 95% credible intervals ----

#------------------------------------------------------------------------------

#credible interval (CI) summaries
#summarize posterior medians and 95% credible intervals
#scale to 100 person years
lm0 <- as.data.frame(t(100 * quantile(MCMC_burn[,1], prob = c(0.5, 0.025, 0.975)))); names(lm0) <- c("lm0","lm0_lower","lm0_upper"); lm0$id <- 1
lmc <- as.data.frame(t(100 * quantile(MCMC_burn[,1]*MCMC_burn[,2], prob = c(0.5, 0.025, 0.975)))); names(lmc) <- c("lmc","lmc_lower","lmc_upper"); lmc$id <- 1
g   <- as.data.frame(t(quantile(MCMC_burn[,2], prob = c(0.5, 0.025, 0.975)))); names(g) <- c("g","g_lower","g_upper"); g$id <- 1
rho <- as.data.frame(t(100 * quantile(MCMC_burn[,3], prob = c(0.5, 0.025, 0.975)))); names(rho) <- c("rho","rho_lower","rho_upper"); rho$id <- 1
tc  <- as.data.frame(t(quantile(MCMC_burn[,4], prob = c(0.5, 0.025, 0.975)))); names(tc) <- c("tc","tc_lower","tc_upper"); tc$id <- 1

result <- join_all(list(lm0,lmc,g,rho,tc), by = 'id', type = 'full')
#set column order
col_order <- c("id","lm0","lm0_lower","lm0_upper","lmc","lmc_lower","lmc_upper",
               "g","g_lower","g_upper","rho","rho_lower","rho_upper","tc","tc_lower","tc_upper")
result <- result[, col_order]
print(result)

#------------------------------------------------------------------------------#

#4.3 - save run-specific results to global environment ----

#------------------------------------------------------------------------------#

#save plots/objects with RUN_ID that we set at the beginning of the run 
#if RUN_ID = "V1", this creates result_V1
assign(paste0("result_", RUN_ID), result, envir = .GlobalEnv)

#save recorded plot object for this run
#example: plot_V1
assign(paste0("plot_", RUN_ID), plot_obj,        envir = .GlobalEnv)

#save median-centered DIC for this run
#example: m2_all_ages_DIC_V1
assign(paste0("m2_all_ages_DIC_", RUN_ID), DIC_median, envir = .GlobalEnv)
assign(paste0("m2_all_ages_par_", RUN_ID), MCMC_par,   envir = .GlobalEnv)

burn_tmp   <- MCMC_burn
lambda_c   <- burn_tmp[,1] * burn_tmp[,2]
burn_named <- cbind(burn_tmp[,1], lambda_c, burn_tmp[,2], burn_tmp[,3],
                    burn_tmp[,4], burn_tmp[,5], burn_tmp[,6])
colnames(burn_named) <- c("lambda_0","lambda_c","gamma","rho","time_c","loglike","prior")
assign(paste0("m2_all_ages_burn_", RUN_ID), burn_named, envir = .GlobalEnv)

assign(paste0("age_seq_",   RUN_ID), age_seq,      envir = .GlobalEnv)
assign(paste0("M2_quant_",  RUN_ID), M2_quant,     envir = .GlobalEnv)
assign(paste0("SP_bins_",   RUN_ID), SP_bins,      envir = .GlobalEnv)
assign(paste0("age_mids_",  RUN_ID), age_bins_mid, envir = .GlobalEnv)


#

#

#


#------------------------------------------------------------------------------#

#4.4 - comparing diagnostics across multiple runs ----

#------------------------------------------------------------------------------#

#diagnostics

#the v1 / v2/ v3 here WILL need to change with your run IDs
#if you have a v4 / v5 / v6 you will to update each thing respectively 

list<-mcmc.list(as.mcmc(m2_all_ages_burn_V1), as.mcmc(m2_all_ages_burn_V2), as.mcmc(m2_all_ages_burn_V3))
gelman.diag(list[,1:5], confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)

coda::effectiveSize(list)
par(mfrow=c(1,5))
coda::traceplot(list[,1:5])

replayPlot(plot_V1)
replayPlot(plot_V2)
replayPlot(plot_V3)

#------------------------------------------------------------------------------#

#4.5 - save recorded plots ----

#------------------------------------------------------------------------------#

#helper function to save recorded plots as TIFF files
save_run_plot <- function(x, file,
                          width_in = 8, height_in = 6,  # bigger figure
                          res = 600,                    # print quality
                          pointsize = 14,               # bigger text
                          bg = "white",
                          compression = "lzw") {
  #accept recordedplot or list with $plot
  plt <- if (inherits(x, "recordedplot")) x else if (!is.null(x$plot)) x$plot
  else stop("Pass a recordedplot or a list with $plot")
  
  tiff(filename = file,
       width = width_in, height = height_in, units = "in",
       res = res, pointsize = pointsize,
       bg = bg, compression = compression)
  
  replayPlot(plt)
  dev.off()
}

#change with run ID, as before
replayPlot(plot_V1)
replayPlot(plot_V2)
replayPlot(plot_V3)

#uncomment to save, change run ID, same as before
#save_run_plot(plot_V1, "elisa2022_run1_plot_V1.tiff")
#save_run_plot(plot_V2, "elisa2022_run2_plot_V2.tiff")
#save_run_plot(plot_V3, "elisa2022_run3_plot_V3.tiff")


#replayPlot(plot_V4)
#replayPlot(plot_V5)
#replayPlot(plot_V6)

#------------------------------------------------------------------------------#

#4.6 - build diagnostics tables after three runs exist ----

#------------------------------------------------------------------------------#

#only build the diagnostic table once all three independent chains exist
if (exists("m2_all_ages_burn_V1") && exists("m2_all_ages_burn_V2") && exists("m2_all_ages_burn_V3")) {
  params <- c("lambda_0","lambda_c","gamma","rho","time_c")
  lst <- mcmc.list(
    as.mcmc(m2_all_ages_burn_V1[, params]),
    as.mcmc(m2_all_ages_burn_V2[, params]),
    as.mcmc(m2_all_ages_burn_V3[, params])
  )
  
  gr  <- gelman.diag(lst, confidence = 0.95, transform = FALSE,
                     autoburnin = TRUE, multivariate = FALSE)$psrf
  GR_point <- setNames(gr[,1], rownames(gr))
  GR_upper <- setNames(gr[,2], rownames(gr))
  
  ess_pooled_vec <- effectiveSize(lst)
  ESS_pooled <- setNames(as.numeric(ess_pooled_vec), names(ess_pooled_vec))
  
  fmt_cri <- function(x, digits = 3) {
    q <- quantile(x, c(0.5, 0.025, 0.975), names = TRUE)
    sprintf(paste0("%.",digits,"f (%.",digits,"f - %.",digits,"f)"),
            unname(q["50%"]), unname(q["2.5%"]), unname(q["97.5%"]))
  }
  ess_run <- function(M, p) as.numeric(effectiveSize(as.mcmc(M[, p])))
  
  one_run_row <- function(M, DIC_val) {
    l0  <- 100 * M[,"lambda_0"]
    lc  <- 100 * M[,"lambda_c"]
    rr  <- 100 * M[,"rho"]
    gam <-       M[,"gamma"]
    tc  <-       M[,"time_c"]
    
    data.frame(
      `lambda_0`      = fmt_cri(l0),
      ESS_per_l0      = round(ess_run(M, "lambda_0"), 1),
      GR_l0           = round(GR_point["lambda_0"], 2),
      `GR Upper_l0`   = round(GR_upper["lambda_0"], 2),
      ESS_pooled_l0   = round(ESS_pooled["lambda_0"], 1),
      
      rho             = fmt_cri(rr),
      ESS_per_rho     = round(ess_run(M, "rho"), 1),
      GR_rho          = round(GR_point["rho"], 2),
      `GR Upper_rho`  = round(GR_upper["rho"], 2),
      ESS_pooled_rho  = round(ESS_pooled["rho"], 1),
      
      DIC             = round(DIC_val, 3),
      
      lambda_c        = fmt_cri(lc),
      GR_lc           = round(GR_point["lambda_c"], 2),
      `GR Upper_lc`   = round(GR_upper["lambda_c"], 2),
      ESS_lc          = round(ESS_pooled["lambda_c"], 1),
      
      gamma           = fmt_cri(gam),
      GR_g            = round(GR_point["gamma"], 2),
      `GR Upper_g`    = round(GR_upper["gamma"], 2),
      ESS_g           = round(ESS_pooled["gamma"], 1),
      
      time_c          = fmt_cri(tc),
      GR_tc           = round(GR_point["time_c"], 2),
      `GR Upper_tc`   = round(GR_upper["time_c"], 2),
      ESS_tc          = round(ESS_pooled["time_c"], 1),
      check.names = FALSE
    )
  }
  
  #open tab_m2 in the global environment if you want to see the data not in the console
  tab_m2 <- rbind(
    one_run_row(m2_all_ages_burn_V1, m2_all_ages_DIC_V1),
    one_run_row(m2_all_ages_burn_V2, m2_all_ages_DIC_V2),
    one_run_row(m2_all_ages_burn_V3, m2_all_ages_DIC_V3)
  )
  rownames(tab_m2) <- c("Run 1","Run 2","Run 3")
  print(tab_m2)
  
  #diagnostics plots
  par(mfrow = c(1,5))
  traceplot(lst)
  
  ##check for time_c
  qc <- function(x) quantile(x, c(0.5, 0.025, 0.975))
  cat("\n(time_c):\n")
  print(list(
    V1_time_c = qc(m2_all_ages_burn_V1[,"time_c"]),
    V2_time_c = qc(m2_all_ages_burn_V2[,"time_c"]),
    V3_time_c = qc(m2_all_ages_burn_V3[,"time_c"])
  ))
  View(tab_m2)
}

