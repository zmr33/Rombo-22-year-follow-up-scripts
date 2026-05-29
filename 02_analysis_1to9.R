################################################################################
# SCR and seroprevalence analysis: ages 1-9
#
# Purpose:
#   Calculate seroprevalence, seroconversion rates, and smoothed
#   age-seroprevalence estimates for children ages 1-9.
#
# Author:
#   Zach Reynolds
#
# Inputs:
#   - elisa_2022.xlsx
#
# Required columns:
#   - sample_id: unique sample or participant ID
#   - res_eu: evaluation unit identifier
#   - age: age in years
#   - seropos: binary serostatus indicator, where 1 = positive and 0 = negative
#
# Outputs:
#   - outputs/analysis_1to9.rds
#   - outputs/scr_and_seroprev_results_1to9.xlsx
################################################################################


# ---- packages ----------------------------------------------------------------
library(readxl)
library(writexl)
library(tidyverse)
library(binom)
library(foreach)
library(sandwich)
library(lmtest)
library(mgcv)


# ---- estimate SCR helper -----------------------------------------------------
# This function was adapted from Ben Arnold's group at UCSF.
#
# It fits a catalytic model with a single seroconversion rate using a GLM.
# The model uses a complementary log-log link and an age offset.
#
# Inputs:
#   serop    Binary seropositivity indicator: 1 = seropositive, 0 = negative
#   agey     Age in years
#   id       Independent unit ID; defaults to individual-level IDs
#   variance Variance estimator: "ML", "robust", or "bootstrap"
#   nboots   Number of bootstrap replicates if variance = "bootstrap"
#   seed     Optional random seed for bootstrap reproducibility
#
# Output:
#   A list containing SCR estimate, uncertainty interval, CI method, and any
#   GLM warning/error message.
estimate_scr_glm <- function(serop, agey, id = 1:length(serop),
                             variance = "ML", nboots = 100, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if ((length(serop) != length(agey)) | (length(serop) != length(id))) {
    stop(
      "\nThe length of serop, agey, and id vectors are not the same.",
      "\nCheck to ensure they are the same length."
    )
  }

  df <- data.frame(id, serop, agey) %>%
    filter(!is.na(serop), !is.na(agey))

  ids <- unique(df$id)

  glm_fit <- tryCatch(
    glm(
      serop ~ 1,
      offset = log(agey),
      data = df,
      family = binomial(link = "cloglog")
    ),
    error = function(cond) cond$message,
    warning = function(cond) cond$message
  )

  if (class(glm_fit)[1] == "character") {
    scr_hat <- NA_real_
    logscr_se <- NA_real_
    scr_min95 <- NA_real_
    scr_max95 <- NA_real_
    scr_ci_method <- NA_character_
    glm_err_warn_msg <- glm_fit
  } else {
    scr_hat <- exp(glm_fit$coefficients)
    logscr_se <- sqrt(summary(glm_fit)$cov.unscaled)
    scr_ci_method <- "Maximum likelihood SE"
    glm_err_warn_msg <- NA_character_

    if (variance == "robust") {
      glm_fit_rb <- lmtest::coeftest(
        glm_fit,
        vcov. = sandwich::vcovCL(glm_fit, cluster = df$id)
      )
      logscr_se <- glm_fit_rb[1, 2]
      scr_ci_method <- "Robust, Huber-White SE"
    }
  }

  scr_min95 <- exp(log(scr_hat) - 1.96 * logscr_se)
  scr_max95 <- exp(log(scr_hat) + 1.96 * logscr_se)

  if (variance == "bootstrap") {
    bsamp <- matrix(
      sample(ids, size = length(ids) * nboots, replace = TRUE),
      nrow = length(ids),
      ncol = nboots
    )

    bootests <- foreach(brep = 1:nboots, .combine = rbind) %do% {
      di <- df %>%
        left_join(data.frame(id = bsamp[, brep]), by = "id")

      glm_fiti <- tryCatch(
        glm(
          serop ~ 1,
          offset = log(agey),
          data = di,
          family = binomial(link = "cloglog")
        ),
        error = function(cond) cond$message,
        warning = function(cond) cond$message
      )

      if (class(glm_fiti)[1] == "character") {
        scr_hati <- NA_real_
      } else {
        scr_hati <- exp(glm_fiti$coefficients)
      }

      data.frame(brep, scr_hati)
    }

    logscr_se <- NA_real_
    scr_min95 <- quantile(bootests$scr_hati, prob = 0.025, na.rm = TRUE)
    scr_max95 <- quantile(bootests$scr_hati, prob = 0.975, na.rm = TRUE)
    scr_ci_method <- paste0("Bootstrap (", nboots, " reps)")
  }

  list(
    scr = scr_hat,
    logscr_se = logscr_se,
    scr_min95 = scr_min95,
    scr_max95 = scr_max95,
    scr_ci_method = scr_ci_method,
    glm_err_warn_msg = glm_err_warn_msg
  )
}


# ---- user settings -----------------------------------------------------------
input_file <- "elisa_2022.xlsx"
dataset_label <- "ELISA, 2022"

age_min <- 1
age_max <- 9
age_band_label <- paste0(age_min, " to ", age_max)

output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE)


# ---- read and validate input data -------------------------------------------
wd <- read_excel(input_file)

required_columns <- c("sample_id", "res_eu", "age", "seropos")
missing_columns <- setdiff(required_columns, names(wd))

if (length(missing_columns) > 0) {
  stop(
    "The input dataset is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

filtered_data <- wd


# ---- prepare analysis dataset ------------------------------------------------
sero_df <- filtered_data %>%
  filter(!is.na(age), !is.na(seropos)) %>%
  filter(age %in% age_min:age_max) %>%
  mutate(res_eu = as.factor(res_eu))

# Optional EU filters:
# sero_df <- sero_df %>% filter(res_eu == "1")   # keep only EU 1
# sero_df <- sero_df %>% filter(res_eu != "3")   # exclude EU 3


# ---- seroprevalence by EU ----------------------------------------------------
seroprev_by_eu <- sero_df %>%
  group_by(res_eu) %>%
  summarise(
    age = paste0(age_band_label, " year olds - combined"),
    n = n(),
    n_pos = sum(seropos, na.rm = TRUE),
    seroprev = 100 * n_pos / n,
    lower_ci = binom.confint(x = n_pos, n = n, method = "wilson")$lower * 100,
    upper_ci = binom.confint(x = n_pos, n = n, method = "wilson")$upper * 100,
    .groups = "drop"
  )


# ---- seroprevalence by EU and age -------------------------------------------
seroprev_by_age <- sero_df %>%
  group_by(res_eu, age) %>%
  summarise(
    n = n(),
    n_pos = sum(seropos, na.rm = TRUE),
    seroprev = 100 * n_pos / n,
    lower_ci = binom.confint(x = n_pos, n = n, method = "wilson")$lower * 100,
    upper_ci = binom.confint(x = n_pos, n = n, method = "wilson")$upper * 100,
    .groups = "drop"
  )


# ---- SCR by EU ---------------------------------------------------------------
# SCR is reported per 100 child-years in the scaled columns.
scr_summary <- foreach(curr_eu = unique(sero_df$res_eu), .combine = rbind) %do% {

  df_temp <- sero_df %>%
    filter(res_eu == curr_eu)

  scr_out <- estimate_scr_glm(
    serop = df_temp$seropos,
    agey = df_temp$age,
    id = df_temp$sample_id,
    variance = "ML"
  )

  tibble(
    res_eu = curr_eu,
    scr = scr_out$scr,
    scr_min95 = scr_out$scr_min95,
    scr_max95 = scr_out$scr_max95,
    scr_ci_method = scr_out$scr_ci_method,
    glm_err_warn_msg = scr_out$glm_err_warn_msg
  )

} %>%
  mutate(
    across(
      c(scr, scr_min95, scr_max95),
      ~ round(. * 100, 2),
      .names = "{.col}_scaled"
    ),
    scr_label = sprintf(
      "<b>SCR:</b> %.2f (%.2f-%.2f)",
      scr_scaled,
      scr_min95_scaled,
      scr_max95_scaled
    )
  ) %>%
  select(
    res_eu,
    scr_scaled,
    scr_min95_scaled,
    scr_max95_scaled,
    scr_label,
    scr,
    scr_min95,
    scr_max95,
    scr_ci_method,
    glm_err_warn_msg
  )


# ---- smoothed age-seroprevalence curves by EU -------------------------------
age_sero_curve_results <- foreach(curr_eu = unique(sero_df$res_eu), .combine = rbind) %do% {

  df_temp <- sero_df %>%
    filter(res_eu == curr_eu)

  gam_fit <- mgcv::gam(
    seropos ~ s(age, k = 4, bs = "cr"),
    data = df_temp,
    family = binomial(link = "logit"),
    method = "REML"
  )

  age_grid <- tibble(age = seq(age_min, age_max, by = 0.1))

  preds <- predict(
    gam_fit,
    newdata = age_grid,
    se.fit = TRUE
  )

  tibble(
    res_eu = curr_eu,
    age = age_grid$age,
    seroprev = plogis(preds$fit) * 100,
    lower_ci = plogis(preds$fit - 1.96 * preds$se.fit) * 100,
    upper_ci = plogis(preds$fit + 1.96 * preds$se.fit) * 100
  )
}

age_sero_curve_results <- age_sero_curve_results %>%
  left_join(scr_summary, by = "res_eu")


# ---- save outputs ------------------------------------------------------------
analysis_output <- list(
  dataset_label = dataset_label,
  age_min = age_min,
  age_max = age_max,
  age_band_label = age_band_label,
  sero_df = sero_df,
  seroprev_by_eu = seroprev_by_eu,
  seroprev_by_age = seroprev_by_age,
  scr_summary = scr_summary,
  age_sero_curve_results = age_sero_curve_results
)

saveRDS(
  analysis_output,
  file.path(output_dir, paste0("analysis_", age_min, "to", age_max, ".rds"))
)

write_xlsx(
  list(
    paste0("SCR Estimates (", age_min, "-", age_max, ")") = scr_summary,
    paste0("Seroprev by EU (", age_min, "-", age_max, ")") = seroprev_by_eu,
    paste0("Seroprev by Age (", age_min, "-", age_max, ")") = seroprev_by_age
  ),
  file.path(output_dir, paste0("scr_and_seroprev_results_", age_min, "to", age_max, ".xlsx"))
)

message(
  "Done. Saved outputs for ages ",
  age_min,
  "-",
  age_max,
  " to: ",
  output_dir
)
