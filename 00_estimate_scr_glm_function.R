#code written by Ben Arnold's group at UCSF

#-------------------------------
# estimate_scr_glm()
#
# fit a catalytic model
# with a single seroconversion rate
# using a Generalized Linear Model (GLM)
# this is equivalent to an SIR model
#
# vectors serop, agey, and id (if specified) must be the same length!
#
# @serop : seropositivity indicator (1=seropos+, 0=seroneg-)
# @agey  : age in years
# @id    : ID variable for independent units (individuals or clusters)
#          if unspecified, assumes 1:length(serop) (individuals)
# @variance : Variance estimator: "ML", "robust", "bootstrap"
#             ML: maximum likelihood (the default), assumes i.i.d obs
#             robust: Huber-White robust SEs, clustered at the level of id
#             bootstrap: non-parametric bootstrap, resampling ids w/ replacement
# @nboots    : number of bootstrap replicates, 
#              default is 100 (use 1000+ for publication-quality CIs)
# 
# returns a list with the following objects
# scr       : ML estimate of the seroconversion rate (SCR)
# scr_se    : ML estimate of the SCR standard error (NA if using a bootstrap)
# scr_min95        : lower 95% confidence interval for the SCR
# scr_max95        : upper 95% confidence interval for the SCR
# scr_ci_method    : method used to estimate the 95% CI (character) 
# glm_err_warn_msg : warning message (in the case of a failed fit) (character)
#-------------------------------

estimate_scr_glm <- function(serop,agey,id=1:length(serop),
                             variance="ML",nboots=100,seed=NULL) {
  
  # set a seed (if specified)
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # confirm that serop, agey, and id are all the same length
  if((length(serop) != length(agey)) | (length(serop) != length(id))) {
    e <- simpleError(paste0("\n\nThe length of serop (",length(serop),"), agey (",length(agey), "), and id (",length(id),") vectors are not the same.\nCheck to ensure they are the same length."))
    stop(e)
  }
  
  # create a data frame and restrict
  # to non-missing observations
  df  <- data.frame(id,serop,agey) %>%
    filter(!is.na(serop) & !is.na(agey))
  
  # identify the number of unique IDs
  ids <- unique(df$id)
  
  # fit the GLM model
  glm_fit <- tryCatch(glm(serop~1, offset = log(agey),
                          data = df,
                          family = binomial(link = "cloglog")),
                      error = function(cond){return(cond$message)},
                      warning = function(cond){return(cond$message)})
  
  if(class(glm_fit)[1] == "character") {
    scr_hat <- NA
    logscr_se  <- NA
    scr_min95 <- NA
    scr_max95 <- NA
    scr_ci_method <- NA
    glm_err_warn_msg <- glm_fit
  } else {
    # retrieve the SCR and the SE of the log SCR
    scr_hat <- exp(glm_fit$coefficients)
    logscr_se <- sqrt(summary(glm_fit)$cov.unscaled)
    scr_ci_method = "Maximum likelihood SE"
    glm_err_warn_msg <- NA
    cat("\n   -------------------------------",
        "\n   GLM model fit",
        "\n   -------------------------------")
    print(summary(glm_fit))
    
    # if variance=="robust"
    # adjust the SEs using A clustered, Huber-White robust variance.
    if(variance == "robust") {
      glm_fit_rb <- coeftest(glm_fit,vcov. = vcovCL(glm_fit,cluster=df$id))
      logscr_se <- glm_fit_rb[1,2]
      scr_ci_method = "Robust, Huber-White SE"
      cat("\n   -------------------------------",
          "\n   GLM model fit, with Robust, SEs",
          "\n   -------------------------------")
      print(glm_fit_rb)
    } 
    
  }
  
  # estimate the 95% CIs
  scr_min95 = exp(log(scr_hat) - 1.96*logscr_se)
  scr_max95 = exp(log(scr_hat) + 1.96*logscr_se)
  
  # if variance = "bootstrap"
  # resample the ID variable with replacement and re-estimate the SCR
  if(variance=="bootstrap") {
    
    cat(paste0("\n\nFitting a catalytic model (SIR) using GLM\n with inference using a nonparametric bootstrap (",nboots," iterations)\n"))
    
    bsamp <- matrix(sample(ids,size=length(ids)*nboots,replace=TRUE),
                    nrow=length(ids),ncol=nboots)
    
    bootests <- foreach(brep=1:nboots,.combine=rbind) %do% {
      if(brep %% 50 == 0) {
        cat(".",brep,"\n")
      }else{
        cat(".")
      }
      
      di <- df %>%
        left_join(data.frame(id=bsamp[,brep]),di,by=c("id")) 
      
      # fit the GLM model
      glm_fiti <- tryCatch(glm(serop~1, offset = log(agey),
                               data = di,
                               family = binomial(link = "cloglog")),
                           error = function(cond){return(cond$message)},
                           warning = function(cond){return(cond$message)})
      
      if(class(glm_fiti)[1] == "character") {
        scr_hati <- NA
      } else {
        scr_hati <- exp(glm_fiti$coefficients)
      }
      res <- data.frame(brep,scr_hati)
      return(res)
      
    }
    # from bootstrap, estimate the 95% CI
    logscr_se <- NA
    scr_min95 <- quantile(bootests$scr_hati,prob=0.025)
    scr_max95 <- quantile(bootests$scr_hati,prob=0.975)
    scr_ci_method <- paste0("Bootstrap (",nboots," reps)")
    
  }
  
  cat("\n\n-----------------------------------------------\n",
      "GLM model fit of the seroconversion rate (SCR)\n",
      "SCR units are seroconversions per child-year \n",
      "----------------------------------------------\n",
      "Non-missing observations:   ", nrow(df),"\n",
      "Number of independent units:",length(ids),"\n",
      "95% CI estimated using", scr_ci_method,"\n",
      "Estimates:\n",
      "SCR estimate:", sprintf("%1.4f",scr_hat),"\n",
      "SCR min 95  :", sprintf("%1.4f",scr_min95),"\n",
      "SCR max 95  :", sprintf("%1.4f",scr_max95),"\n",
      "----------------------------------------------\n")
  
  return(list(scr=scr_hat,logscr_se=logscr_se,scr_min95=scr_min95,scr_max95=scr_max95,scr_ci_method=scr_ci_method,glm_err_warn_msg=glm_err_warn_msg))
  
}

