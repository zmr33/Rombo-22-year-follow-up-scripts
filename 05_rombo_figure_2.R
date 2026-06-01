#'
#' script to generate district level age-seroprevalence curves for trachoma
#' 
#' see: '01_rombo_analysis_1to5.R' and '02_rombo_analysis_1to9.R'
#' for proper set up 
#'
#' 
#' written by zach reynolds, 
#' based off of original code by Ben Arnold's group at UCSF
#'

#------------------------------------------------------------------------------#

#load packages ----

#------------------------------------------------------------------------------#

library(tidyverse)
library(ggtext)
library(patchwork)
library(grid)
library(scales)
library(writexl)

#------------------------------------------------------------------------------#

#create folder for saving figures ----

#------------------------------------------------------------------------------#

figure_dir <- "figures"
dir.create(figure_dir, showWarnings = FALSE)

#------------------------------------------------------------------------------#

#read analysis outputs ----

#------------------------------------------------------------------------------#

analysis_output_1to5 <- readRDS("analysis_output_1to5.rds")
analysis_output_1to9 <- readRDS("analysis_output_1to9.rds")

mba_1to5   <- analysis_output_1to5$mba_1to5
elisa_1to5 <- analysis_output_1to5$elisa_1to5

mba_1to9   <- analysis_output_1to9$mba_1to9
elisa_1to9 <- analysis_output_1to9$elisa_1to9

#------------------------------------------------------------------------------#

#plotting function ----

#------------------------------------------------------------------------------#

make_age_sero_plot <- function(analysis_output, age_min, age_max, tag = NULL) {
  
  #pull dataset label from the analysis object
  #older versions used dataset; newer versions may use dataset_label
  dataset_label <- analysis_output$dataset_label
  
  if (is.null(dataset_label)) {
    dataset_label <- analysis_output$dataset
  }
  
  x_breaks <- age_min:age_max
  
  age_sero_curve_results <- analysis_output$age_sero_curve_results
  seroprev_by_age        <- analysis_output$seroprev_by_age
  scr_summary            <- analysis_output$scr_summary
  
  #make SCR label robust to whether scr_label already includes "SCR:"
  scr_summary <- scr_summary %>%
    mutate(
      scr_label_plot = if_else(
        str_detect(scr_label, "^SCR:"),
        str_replace(scr_label, "^SCR:\\s*", "<b>SCR:</b> "),
        paste0("<b>SCR:</b> ", scr_label)
      )
    )
  
  age_sero_curve_results %>%
    ggplot(aes(x = age, group = res_eu)) +
    
    #observed points
    geom_point(
      data = seroprev_by_age,
      aes(x = age, y = seroprev / 100),
      inherit.aes = FALSE,
      shape = 16,
      size = 4,
      color = "black"
    ) +
    
    #wilson 95% confidence intervals around observed age-specific seroprevalence
    geom_errorbar(
      data = seroprev_by_age,
      aes(x = age, ymin = lower_ci / 100, ymax = upper_ci / 100),
      inherit.aes = FALSE,
      width = 0.2,
      linewidth = 0.8,
      color = "black"
    ) +
    
    #estimated age-seroprevalence curve
    geom_ribbon(
      aes(ymin = lower_ci / 100, ymax = upper_ci / 100),
      fill = "blue",
      alpha = 0.30,
      color = NA,
      na.rm = TRUE
    ) +
    
    #smoothed line on top
    geom_line(
      aes(y = seroprev / 100),
      color = "blue",
      linewidth = 1.2
    ) +
    
    #SCR label
    geom_richtext(
      data = scr_summary,
      aes(x = age_min, y = 0.93, label = scr_label_plot),
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
}

#------------------------------------------------------------------------------#

#create individual panels ----

#------------------------------------------------------------------------------#

mba_1to5_plot   <- make_age_sero_plot(mba_1to5, age_min = 1, age_max = 5, tag = "2A.")
elisa_1to5_plot <- make_age_sero_plot(elisa_1to5, age_min = 1, age_max = 5, tag = "2B.")

mba_1to9_plot   <- make_age_sero_plot(mba_1to9, age_min = 1, age_max = 9, tag = "2C.")
elisa_1to9_plot <- make_age_sero_plot(elisa_1to9, age_min = 1, age_max = 9, tag = "2D.")

#------------------------------------------------------------------------------#

#create 2 x 1 figures ----

#------------------------------------------------------------------------------#

fig_1to5_2x1 <- mba_1to5_plot | elisa_1to5_plot
fig_1to9_2x1 <- mba_1to9_plot | elisa_1to9_plot

print(fig_1to5_2x1)
print(fig_1to9_2x1)

#------------------------------------------------------------------------------#

#create 2 x 2 figure ----

#------------------------------------------------------------------------------#

fig_2x2 <- (mba_1to5_plot | elisa_1to5_plot) /
  (mba_1to9_plot | elisa_1to9_plot)

print(fig_2x2)

#------------------------------------------------------------------------------#

#save figures ----

#------------------------------------------------------------------------------#

ggsave(
  filename = file.path(figure_dir, "age_seroprev_1to5_2x1.tiff"),
  plot = fig_1to5_2x1,
  device = "tiff",
  width = 10,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = file.path(figure_dir, "age_seroprev_1to9_2x1.tiff"),
  plot = fig_1to9_2x1,
  device = "tiff",
  width = 10,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = file.path(figure_dir, "age_seroprev_2x2.tiff"),
  plot = fig_2x2,
  device = "tiff",
  width = 10,
  height = 8,
  dpi = 300,
  compression = "lzw"
)

#------------------------------------------------------------------------------#

#save plotting tables ----

#------------------------------------------------------------------------------#

plot_output_tables <- list(
  "MBA 2012 SCR 1-5" = mba_1to5$scr_summary,
  "MBA 2012 Age Prev 1-5" = mba_1to5$seroprev_by_age,
  "MBA 2012 Curve 1-5" = mba_1to5$age_sero_curve_results,
  
  "ELISA 2022 SCR 1-5" = elisa_1to5$scr_summary,
  "ELISA 2022 Age Prev 1-5" = elisa_1to5$seroprev_by_age,
  "ELISA 2022 Curve 1-5" = elisa_1to5$age_sero_curve_results,
  
  "MBA 2012 SCR 1-9" = mba_1to9$scr_summary,
  "MBA 2012 Age Prev 1-9" = mba_1to9$seroprev_by_age,
  "MBA 2012 Curve 1-9" = mba_1to9$age_sero_curve_results,
  
  "ELISA 2022 SCR 1-9" = elisa_1to9$scr_summary,
  "ELISA 2022 Age Prev 1-9" = elisa_1to9$seroprev_by_age,
  "ELISA 2022 Curve 1-9" = elisa_1to9$age_sero_curve_results
)

write_xlsx(
  x = plot_output_tables,
  path = file.path(figure_dir, "age_seroprev_plotting_tables.xlsx")
)
