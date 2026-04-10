# 01_EDA_environmental_plots.R
# Author: Rezvan Hatami
# Original date: Oct 2016
# Revised date: 01 Mar 2026
# Objective: This script reads the raw Wang biological and environmental data and produces the exploratory environmental figures used in early workflow checks.
# It is fully self-sufficient, writes figure outputs to figures/01, and saves the data objects needed by the next script to output/.
# Set the working directory to the downloaded project folder before running this script.

# ---- 1. Load modern packages used in this script -------------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(download.file.method = "wininet")

required_pkgs <- c("ggplot2", "dplyr", "readr", "fs", "here")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE, method = "wininet")
}

library(ggplot2)
library(dplyr)
library(readr)
library(fs)
library(here)
library(patchwork)

# ---- 2. Define project paths and create output folders -------------------------------------------------------

# Set the working directory to the downloaded project folder, For example:
# project_dir <- "C:/Users/yourname/Documents/Causality_Wang"  
project_dir <- here::here() # Replace this line with your own project path
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "output")
figure_dir <- file.path(project_dir, "figures", "01")
dir_create(output_dir)
dir_create(figure_dir)

# ---- 3. Read the raw input files required by this script -----------------------------------------------------
wangbug <- read_csv(file.path(data_dir, "wangbug.csv"), show_col_types = FALSE)
wangenv <- read_csv(file.path(data_dir, "wangenv.csv"), show_col_types = FALSE)
diversity <- read_csv(file.path(data_dir, "diversity.csv"), show_col_types = FALSE)

# ---- 4. Prepare date fields used later in the rainfall and flow figures -------------------------------------
if ("date1" %in% names(wangenv)) wangenv <- wangenv %>% mutate(date = as.Date(date1, format = "%d/%m/%Y"))
if (!"date" %in% names(wangenv) && "date1" %in% names(wangenv)) wangenv$date <- as.Date(wangenv$date1, format = "%d/%m/%Y")
if ("date" %in% names(wangenv) && all(is.na(wangenv$date)) && "date1" %in% names(wangenv)) wangenv$date <- as.Date(wangenv$date1, format = "%d-%m-%Y")

# ---- 5. Define common plotting settings used across all environmental figures --------------------------------
sample_days <- c(1, 126, 260, 336, 518)
day_palette <- c("1" = "black", "126" = "gold2", "260" = "blue", "336" = "green3", "518" = "red")
day_labels <- c("1" = "Dec 13", "126" = "April 14", "260" = "Aug 2014", "336" = "Nov 2014", "518" = "May 2015")

make_spatial_plot <- function(data, var_name, y_label, panel_title, y_limits) {
  ggplot(filter(data, day %in% sample_days), aes(x = dist, y = .data[[var_name]], color = factor(day), group = factor(day))) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.4) +
    geom_vline(xintercept = 4, linetype = 2) +
    scale_color_manual(values = day_palette, labels = day_labels, name = "months") +
    coord_cartesian(ylim = y_limits) +
    labs(x = "Distance (km)", y = y_label, title = panel_title) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(size = 11))
}

save_panel_figure <- function(plot_list, filename, width = 14, height = 14) {
  wrapped_plot <- wrap_plots(plotlist = plot_list, ncol = 2, guides = "collect") & theme(legend.position = "top")
  ggsave(filename = file.path(figure_dir, filename), plot = wrapped_plot, width = width, height = height, dpi = 800, bg = "transparent")
}

# ---- 6. Generate and save nutrient and phosphate figure ------------------------------------------------------
nutrient_plots <- list(
  make_spatial_plot(wangenv, "tn", "TN (mg/L)", "a) Total Nitrogen plotted against spatial position", c(0, 3.5)),
  make_spatial_plot(wangenv, "no2", "NO2 (mg/L)", "b) Nitrite plotted against spatial position", c(0, 0.7)),
  make_spatial_plot(wangenv, "no3", "NO3 (mg/L)", "c) Nitrate plotted against spatial position", c(0, 0.7)),
  make_spatial_plot(wangenv, "nh3", "NH3 (mg/L)", "d) Ammonia plotted against spatial position", c(0, 0.7)),
  make_spatial_plot(wangenv, "op", "Orthophosphate (mg/L)", "e) Orthophosphate plotted against spatial position", c(0, 0.08)),
  make_spatial_plot(wangenv, "tp", "TP", "f) Total Phosphorus plotted against spatial position", c(0, 0.7))
)
save_panel_figure(nutrient_plots, "nutrients.tif")

# ---- 7. Generate and save heavy metals, COD, and antimony figure ------------------------------------------------
metal_plots <- list(
  make_spatial_plot(wangenv, "cd", "Cd (mg/kg)", "a) Cadmium plotted against spatial position", c(0, 1.5)),
  make_spatial_plot(wangenv, "cr", "Cr (mg/kg)", "b) Chromium plotted against spatial position", c(0, 10)),
  make_spatial_plot(wangenv, "cu", "Cu (mg/kg)", "c) Copper plotted against spatial position", c(0, 16)),
  make_spatial_plot(wangenv, "zn", "Zn (mg/kg)", "d) Zinc plotted against spatial position", c(0, 120)),
  make_spatial_plot(wangenv, "sb", "Antimony (mg/L)", "e) Antimony plotted against spatial position", c(0, 0.2)),
  make_spatial_plot(wangenv, "cod", "COD", "f) COD (mg/L) plotted against spatial position", c(0, 100))
)
save_panel_figure(metal_plots, "heavy metals, cod, and antimony.tif")

# ---- 8. Generate and save pH, temperature, conductivity, alkalinity, TOC, and chlorophyll figure -------------
physchem_plots <- list(
  make_spatial_plot(wangenv, "temp", "Temperature (-C)", "a) Temperature plotted against spatial position", c(5, 25)),
  make_spatial_plot(wangenv, "ph2", "pH", "b) pH plotted against spatial position", c(6.5, 8.5)),
  make_spatial_plot(wangenv, "alk", "Alkalinity (mg/L)", "c) Alkalinity plotted against spatial position", c(0, 200)),
  make_spatial_plot(wangenv, "chla", "Chlorophyll A (mg/L)", "d) Chlorophyll A plotted against spatial position", c(0, 40)),
  make_spatial_plot(wangenv, "toc", "TOC (mg/L)", "e) Total Organic Carbon plotted against spatial position", c(0, 12)),
  make_spatial_plot(wangenv, "cond", "Conductivity (mg/L)", "f) Conductivity plotted against spatial position", c(0, 700))
)
save_panel_figure(physchem_plots, "ph, temperature, conductivity, alkalinity, toc, and chla.jpg")

# ---- 9. Generate and save dissolved oxygen, canopy, velocity, turbidity, and habitat figure -------------------
habitat_plots <- list(
  make_spatial_plot(wangenv, "vel", "Velocity (m/s)", "a) Velocity plotted against spatial position", c(0, 0.8)),
  make_spatial_plot(wangenv, "turb", "Turbidity (NTU)", "b) Turbidity plotted against spatial position", c(0, 100)),
  make_spatial_plot(wangenv, "sed", "Sediment size (mm)", "c) Sediment size plotted against spatial position", c(0.5, 2.5)),
  make_spatial_plot(wangenv, "canop", "Canopy cover (%)", "d) Canopy cover plotted against spatial position", c(0, 100)),
  make_spatial_plot(wangenv, "do", "Dissolved Oxygen (mg/L)", "e) DO plotted against spatial position", c(4, 20)),
  make_spatial_plot(wangenv, "cfpom", "CPOM/ FPOM", "f) CPOM/FPOM plotted against spatial position", c(0, 50))
)
save_panel_figure(habitat_plots, "DO, canop, velocity and turbidity.jpg")

# ---- 10. Select a valid reference site and generate rainfall and flow rate figure ------------------------------
reference_dist <- wangenv %>% filter(!is.na(date), !is.na(dayflow), !is.na(dflow), !is.na(rain3)) %>% arrange(dist, date) %>% slice(1) %>% pull(dist)
reference_data <- wangenv %>% filter(dist == reference_dist) %>% arrange(date)
flow_rain_plots <- list(
  ggplot(reference_data, aes(date, dayflow)) + geom_point(size = 2.4) + labs(x = "day", y = "Creek flow rate (ML/day)", title = "Average daily Creek flow rate plotted against time") + theme_bw(base_size = 11),
  ggplot(reference_data, aes(date, dflow)) + geom_point(size = 2.4, color = "blue") + labs(x = "day", y = "Discharge flow rate (ML/day)", title = "Daily effluent discharge flow rate plotted against time") + theme_bw(base_size = 11),
  ggplot(reference_data, aes(date, rain3)) + geom_point(size = 2.4, color = "red") + labs(x = "day", y = "Rainfall(mm)", title = "Three-month average rainfall plotted against time") + theme_bw(base_size = 11)
)
ggsave(filename = file.path(figure_dir, "rain and flow rate.tif"), plot = wrap_plots(plotlist = flow_rain_plots, ncol = 1), width = 12, height = 14, dpi = 800, bg = "transparent")

# ---- 11. Generate and save one-month rainfall figure ----------------------------------------------------------
rain_month_plot <- ggplot(reference_data, aes(date, rain2)) + geom_point(size = 2.4, color = "blue") + labs(x = "day", y = "Rainfall(mm)", title = "Average rainfall (a month) plotted against time") + theme_bw(base_size = 11)
ggsave(filename = file.path(figure_dir, "rainfall (a month).tif"), plot = rain_month_plot, width = 6, height = 3.5, dpi = 800, bg = "transparent")

# ---- 12. Save the data objects required by the next script ----------------------------------------------------
write_rds(wangbug, file.path(output_dir, "01_wangbug_raw.rds"))
write_rds(wangenv, file.path(output_dir, "01_wangenv_prepped.rds"))
write_rds(diversity, file.path(output_dir, "01_diversity_raw.rds"))
