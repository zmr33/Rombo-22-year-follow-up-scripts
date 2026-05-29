################################################################################
# Plot SCR and age-specific seroprevalence figures
#
# Purpose:
#   Build and save figures using the analysis outputs created by:
#     - 01_analysis_1to5.R
#     - 02_analysis_1to9.R
#
# This script does not rerun the analysis. It only reads saved .rds files,
# creates plots, and saves figures.
################################################################################


# ---- packages ----------------------------------------------------------------
library(tidyverse)
library(ggtext)
library(patchwork)
library(grid)
library(scales)
library(writexl)


# ---- user settings -----------------------------------------------------------
output_dir <- "outputs"
figure_dir <- "figures"

dir.create(figure_dir, showWarnings = FALSE)


# ---- read analysis outputs ---------------------------------------------------
analysis_1to5 <- readRDS(file.path(output_dir, "analysis_1to5.rds"))
analysis_1to9 <- readRDS(file.path(output_dir, "analysis_1to9.rds"))


# ---- plotting helper ---------------------------------------------------------
make_age_sero_plot <- function(analysis_output, tag = NULL) {

  age_min <- analysis_output$age_min
  age_max <- analysis_output$age_max
  dataset_label <- analysis_output$dataset_label

  x_breaks <- age_min:age_max

  age_sero_curve_results <- analysis_output$age_sero_curve_results
  seroprev_by_age <- analysis_output$seroprev_by_age
  scr_summary <- analysis_output$scr_summary

  p <- age_sero_curve_results %>%
    ggplot(aes(x = age, group = res_eu)) +

    # Observed seroprevalence estimates and Wilson confidence intervals.
    geom_point(
      data = seroprev_by_age,
      aes(x = age, y = seroprev / 100),
      inherit.aes = FALSE,
      shape = 15,
      size = 4,
      color = "black",
      na.rm = TRUE
    ) +
    geom_errorbar(
      data = seroprev_by_age,
      aes(x = age, ymin = lower_ci / 100, ymax = upper_ci / 100),
      inherit.aes = FALSE,
      width = 0.2,
      linewidth = 0.8,
      color = "black",
      na.rm = TRUE
    ) +

    # Smoothed GAM seroprevalence curve and model-based uncertainty interval.
    geom_ribbon(
      aes(ymin = lower_ci / 100, ymax = upper_ci / 100),
      fill = "blue",
      alpha = 0.30,
      color = NA,
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = seroprev / 100),
      color = "blue",
      linewidth = 1.2,
      na.rm = TRUE
    ) +

    # SCR label.
    geom_richtext(
      data = scr_summary,
      aes(x = age_min, y = 0.93, label = scr_label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 5,
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt")
    ) +

    scale_x_continuous(
      breaks = x_breaks,
      limits = c(age_min - 0.1, age_max + 0.1)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    labs(
      x = "Age (years)",
      y = "Proportion seropositive",
      title = paste0(dataset_label, " (ages ", age_min, "-", age_max, ")"),
      tag = tag
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.tag = element_text(face = "bold", size = 12),
      plot.tag.position = c(0.02, 0.98),
      legend.position = "none"
    ) +
    coord_cartesian(clip = "off")

  # If your dataset contains multiple EUs and you want one panel per EU,
  # uncomment the next line:
  # p <- p + facet_wrap(~ res_eu)

  p
}


# ---- create figures ----------------------------------------------------------
age_sero_plot_1to5 <- make_age_sero_plot(analysis_1to5, tag = "A.")
age_sero_plot_1to9 <- make_age_sero_plot(analysis_1to9, tag = "B.")

combined_age_sero_plot <- age_sero_plot_1to5 / age_sero_plot_1to9


# ---- print figures to viewer -------------------------------------------------
print(age_sero_plot_1to5)
print(age_sero_plot_1to9)
print(combined_age_sero_plot)


# ---- save figures ------------------------------------------------------------
ggsave(
  file.path(figure_dir, "age_seroprev_1to5.tiff"),
  age_sero_plot_1to5,
  device = "tiff",
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "age_seroprev_1to9.tiff"),
  age_sero_plot_1to9,
  device = "tiff",
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "age_seroprev_combined_1to5_1to9.tiff"),
  combined_age_sero_plot,
  device = "tiff",
  width = 7,
  height = 10,
  dpi = 300
)


# ---- optional combined workbook ---------------------------------------------
# This keeps the old single-workbook behavior without rerunning analysis.
combined_output_list <- list(
  "SCR Estimates (1-5)"   = analysis_1to5$scr_summary,
  "Seroprev by EU (1-5)"  = analysis_1to5$seroprev_by_eu,
  "Seroprev by Age (1-5)" = analysis_1to5$seroprev_by_age,

  "SCR Estimates (1-9)"   = analysis_1to9$scr_summary,
  "Seroprev by EU (1-9)"  = analysis_1to9$seroprev_by_eu,
  "Seroprev by Age (1-9)" = analysis_1to9$seroprev_by_age
)

write_xlsx(
  combined_output_list,
  file.path(output_dir, "scr_and_seroprev_combined.xlsx")
)

message("Done. Figures saved to: ", figure_dir)
