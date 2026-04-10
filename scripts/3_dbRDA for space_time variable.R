# 03_dbRDA_space_time_variable.R
# Author: Rezvan Hatami
# Date: Oct 2016
# Revised date: 03 March 2026
# This script reads the PCO outputs from Script 02, fits the main dbRDA spatiotemporal model with effluent terms, and decomposes observed ordination structure into fitted and residual components.
# It also generates the main spatial, temporal, prediction, and residual figures, then saves the model objects and derived data required by the next script.

# Workflow in this script:
# 1. Load the packages required for modelling, plotting, and file handling.
# 2. Define the project, output, and figure folders used by this script.
# 3. Read the prepared inputs created by Scripts 01 and 02.
# 4. Combine the environmental data and PCO scores into a single analysis table.
# 5. Generate Figure 4 showing PCO1 to PCO6 across spatial position.
# 6. Generate Figure 5 showing PCO1 to PCO6 across time, plus the separate PCO7 time plot.
# 7. Fit the main space-time-effluent dbRDA model and the reduced spatiotemporal model.
# 8. Test model terms and compare the reduced and full dbRDA models.
# 9. Derive fitted and residual PCO components from the dbRDA design matrix.
# 10. Generate prediction grids and produce Figure 6 with observed values and fitted model lines.
# 11. Produce residual diagnostics and residual figures after modelling.
# 12. Save the model objects, fitted values, residual values, and derived tables required by the next script.

# ---- 1. Load packages required for modelling, plotting, and file handling -------------------------------------
rm(list = ls())
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(download.file.method = "wininet")
required_pkgs <- c("vegan", "permute", "dplyr", "readr", "tibble", "fs", "here", "lattice")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, dependencies = TRUE, method = "wininet")
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---- 2. Define project paths and create output folders for this script ----------------------------------------
project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
figure_dir <- fs::path(project_dir, "figures", "03")
fs::dir_create(output_dir)
fs::dir_create(figure_dir)

# ---- 3. Read the prepared inputs created by Scripts 01 and 02 -------------------------------------------------
wangbug_file <- fs::path(output_dir, "01_wangbug_raw.rds")
wangenv_file <- fs::path(output_dir, "01_wangenv_prepped.rds")
wangbug_bc_file <- fs::path(output_dir, "02_wangbug_BC.rds")
pwangbug_file <- fs::path(output_dir, "02_pwangbug_scores.rds")
diagnostics_file <- fs::path(output_dir, "02_pco_diagnostics.rds")
if (!file.exists(wangbug_file)) stop("Missing input file: ", wangbug_file)
if (!file.exists(wangenv_file)) stop("Missing input file: ", wangenv_file)
if (!file.exists(wangbug_bc_file)) stop("Missing input file: ", wangbug_bc_file)
if (!file.exists(pwangbug_file)) stop("Missing input file: ", pwangbug_file)
wangbug <- readRDS(wangbug_file)
wangenv <- readRDS(wangenv_file)
wangbug_BC <- readRDS(wangbug_bc_file)
pwangbug <- readRDS(pwangbug_file)
if (file.exists(diagnostics_file)) pco_diagnostics <- readRDS(diagnostics_file)
p <- ncol(pwangbug)
if (!"date" %in% names(wangenv) && "date1" %in% names(wangenv)) wangenv$date <- as.Date(as.character(wangenv$date1), format = "%d/%m/%Y")
if ("date" %in% names(wangenv) && inherits(wangenv$date, "character")) wangenv$date <- as.Date(as.character(wangenv$date), format = "%d-%m-%Y")
bugenv <- cbind(wangenv, pwangbug)
pco_varpercent <- if (exists("pco_diagnostics") && !is.null(pco_diagnostics$pco_varpercent)) pco_diagnostics$pco_varpercent else round(100 * apply(pwangbug, 2, var) / sum(apply(pwangbug, 2, var)), 1)

# ---- 4. Create helper vectors used in repeated spatial and temporal panels ------------------------------------
day_values <- c(1, 126, 260, 336, 518)
day_cols <- c("black", "gold2", "blue", "green3", "red")
day_pch <- c(19, 15, 8, 17, 18)
day_lty <- 1:5
day_labels <- c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015")
dist_values <- sort(unique(bugenv$dist))
dist_cols <- c("black", "deep pink", "blue", "green", "red", "gold2", "lightseagreen", "sienna2")[seq_along(dist_values)]
dist_pch <- c(19, 15, 8, 17, 18, 18, 18, 18)[seq_along(dist_values)]
dist_lty <- c(2, 3, 4, 5, 6, 7, 6, 6)[seq_along(dist_values)]
if ("date" %in% names(bugenv)) {
  bugenv$date <- as.Date(as.character(bugenv$date), format = "%Y-%m-%d")
  if (all(is.na(bugenv$date)) && "date1" %in% names(bugenv)) bugenv$date <- as.Date(as.character(bugenv$date1), format = "%d-%m-%Y")
  mydates <- bugenv[bugenv$dist == min(dist_values), "date"]
  mydaterange <- c(as.POSIXct(min(mydates, na.rm = TRUE)), as.POSIXct(max(mydates, na.rm = TRUE)))
}

# ---- 5. Generate Figure 4 showing PCO1 to PCO6 across spatial position ---------------------------------------
tiff(file = fs::path(figure_dir, "03_Figure_4.tif"), width = 14, height = 14, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow = c(3, 2), mar = c(5, 5, 4, 2), cex = 0.9, cex.axis = 0.8, las = 1)
for (axis_id in 1:6) {
  y_name <- paste0("pco", axis_id)
  y_lab <- if (axis_id == 1) "PCO1 score" else paste0("PCO", axis_id)
  ttl <- if (axis_id == 1) "PCO1 plotted against spatial position" else paste0(letters[axis_id], ") PCO", axis_id, " plotted against spatial position")
  plot(bugenv[bugenv$day == day_values[1], "dist"], bugenv[bugenv$day == day_values[1], y_name], type = "b", ylim = c(-0.6, 1), pch = day_pch[1], cex = 1.5, xlab = "Distance(km)", ylab = y_lab, main = ttl, lwd = 3, lty = day_lty[1])
  for (j in 2:length(day_values)) points(bugenv[bugenv$day == day_values[j], "dist"], bugenv[bugenv$day == day_values[j], y_name], type = "b", col = day_cols[j], pch = day_pch[j], lwd = 3, cex = if (j == 5) 2 else 1.5, lty = day_lty[j])
  legend("topright", inset = c(0, 0), legend = day_labels, lty = c(1, 5), pch = day_pch, lwd = 2, col = day_cols, ncol = 2, horiz = FALSE, cex = 0.9, title = "months")
  abline(v = 4, lty = 2)
}
dev.off()

# ---- 6. Generate Figure 5 showing PCO1 to PCO6 through time, plus a separate PCO7 panel ----------------------
tiff(file = fs::path(figure_dir, "03_Figure_5.tif"), width = 14, height = 14, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow = c(3, 2), mar = c(2, 4.5, 2, 0.5), cex = 1.5, cex.axis = 0.7, las = 1, cex.main = 0.8)
for (axis_id in 1:6) {
  y_name <- paste0("pco", axis_id)
  ttl <- paste0(letters[axis_id], ") PCO", axis_id, " plotted against time")
  plot(mydates, bugenv[bugenv$dist == dist_values[1], y_name], type = "b", ylim = c(-0.8, 0.8), pch = 19, cex = 1, xlab = "day", ylab = paste0("PCO", axis_id), main = ttl, lwd = 3, lty = 2, cex.main = 1, xaxt = "n")
  if (length(dist_values) > 1) for (j in 2:length(dist_values)) points(mydates, bugenv[bugenv$dist == dist_values[j], y_name], type = "b", col = dist_cols[j], pch = dist_pch[j], lwd = 3, cex = if (j >= 5) 1.5 else 1, lty = dist_lty[j])
  axis.Date(side = 1, mydates, at = seq(mydaterange[1], mydaterange[2], by = "month"), format = "%b-%y")
  legend("bottomright", inset = c(0, 0), legend = levels(as.factor(bugenv$dist)), lty = 1, lwd = 2, col = dist_cols, ncol = 4, horiz = FALSE, cex = 0.6, title = "distance(km)")
}
dev.off()
if (ncol(pwangbug) >= 7) {
  tiff(file = fs::path(figure_dir, "03_PCO7_time.tif"), width = 6, height = 3.5, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
  par(mar = c(4, 4, 2.5, 2) + 0.1, cex = 0.6, cex.axis = 0.8, las = 1)
  plot(mydates, bugenv[bugenv$dist == dist_values[1], "pco7"], type = "b", ylim = c(-0.7, 0.7), pch = 19, cex = 1.5, xlab = "day", ylab = "PCO7", main = "PCO7 plotted against time", lwd = 3, lty = 2, xaxt = "n")
  if (length(dist_values) > 1) for (j in 2:length(dist_values)) points(mydates, bugenv[bugenv$dist == dist_values[j], "pco7"], type = "b", col = dist_cols[j], pch = dist_pch[j], lwd = 3, cex = if (j >= 5) 2 else 1.5, lty = dist_lty[j])
  axis.Date(side = 1, mydates, at = seq(mydaterange[1], mydaterange[2], by = "month"), format = "%b-%y")
  legend("bottomright", inset = c(0, 0), legend = levels(as.factor(bugenv$dist)), lty = 1, lwd = 2, col = dist_cols, ncol = 4, horiz = FALSE, cex = 0.8, title = "distance(km)")
  dev.off()
}

# ---- 7. Fit the main dbRDA space-time-effluent model and the reduced spatiotemporal model ---------------------
wangenv$eff <- ifelse(wangenv$dist < 4, 0, 1)
wangenv$time <- as.factor(wangenv$time)
wang_cap1 <- vegan::capscale(formula = wangbug_BC ~ dist + time + dist:time + eff + eff:time + eff:dist:time, data = wangenv, comm = wangbug, add = TRUE, na.action = na.omit)
wang_anova1 <- anova(wang_cap1, by = "term", permutations = permute::how(nperm = 9999))
wang_cap1b <- vegan::capscale(formula = wangbug_BC ~ dist + time + dist:time, data = wangenv, comm = wangbug, add = TRUE, na.action = na.omit)
wang_anova1b <- anova(wang_cap1b, by = "term")
wang_model_compare <- anova(wang_cap1b, wang_cap1, permutations = permute::how(nperm = 9999))

# ---- 8. Construct the design matrix and derive fitted and residual PCO values ---------------------------------
design <- scale(model.matrix(~ -1 + dist + time + dist:time + eff + eff:time + eff:dist:time, data = wangenv), center = FALSE, scale = FALSE)
pco_beta <- qr.coef(qr(design), as.matrix(pwangbug))
pco_predict <- qr.fitted(qr(design), as.matrix(pwangbug))
pco_resid <- as.matrix(pwangbug) - pco_predict
colnames(pco_predict) <- paste0("pred", 1:p)
colnames(pco_resid) <- paste0("res", 1:p)
envpcopred <- cbind(wangenv, as.data.frame(pco_predict))
envpcores <- cbind(wangenv, as.data.frame(pco_resid))
alldata <- cbind(envpcopred, pwangbug)
ss_total <- sum(diag(var(as.matrix(pwangbug)))) * p
ss_predict <- sum(diag(var(as.matrix(pco_predict)))) * p
ss_res <- sum(diag(var(as.matrix(pco_resid)))) * p
aic_value <- nrow(wangenv) * log(ss_res / nrow(wangenv)) + 2 * ncol(design)
bic_value <- nrow(wangenv) * log(ss_res / nrow(wangenv)) + log(nrow(wangenv)) * ncol(design)

# ---- 9. Generate a regular space-time grid and compute model predictions for plotting -------------------------
xygrid <- expand.grid(seq(1, 5, 1), seq(0, 7.1, 0.1))
colnames(xygrid) <- c("time", "dist")
for (i in 1:nrow(xygrid)) {
  if (xygrid$time[i] == 1) xygrid$day[i] <- 1
  if (xygrid$time[i] == 2) xygrid$day[i] <- 125
  if (xygrid$time[i] == 3) xygrid$day[i] <- 260
  if (xygrid$time[i] == 4) xygrid$day[i] <- 336
  if (xygrid$time[i] == 5) xygrid$day[i] <- 518
}
xygrid$eff <- ifelse(xygrid$dist < 4, 0, 1)
colnames(xygrid) <- c("time", "dist", "day", "eff")
xygrid$time <- as.factor(xygrid$time)
design2 <- scale(model.matrix(~ -1 + dist + time + dist:time + as.factor(eff) + as.factor(eff):time + as.factor(eff):dist:time, data = xygrid), center = FALSE, scale = FALSE)
pco_predict2 <- as.matrix(design2) %*% pco_beta
colnames(pco_predict2) <- paste0("pred", 1:p)
pco_predict2 <- cbind(xygrid, as.data.frame(pco_predict2))
wang_table2_lm <- lm(pco1 ~ dist + time + dist:time + eff + eff:time + eff:dist:time, data = alldata)

# ---- 10. Generate Figure 6 showing observed PCOs with spatiotemporal model predictions ------------------------
tiff(file = fs::path(figure_dir, "03_Figure_6.tif"), width = 14, height = 14, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2, 0.5), cex = 1.5, cex.axis = 0.7, las = 1, cex.main = 0.8)

for (axis_id in 1:6) {
  pco_no <- paste0("pco", axis_id)
  pred_no <- paste0("pred", axis_id)
  pco_lab <- paste0("PCO", axis_id, " (", pco_varpercent[axis_id], "% of total variation)")
  
  plot(wangenv$dist[wangenv$time == levels(wangenv$time)[1]], pwangbug[wangenv$time == levels(wangenv$time)[1], pco_no],
       type = "p", pch = 19, col = day_cols[1], xlab = "Distance (km)", ylab = pco_lab, ylim = c(-0.6, 1),
       las = 1, cex.lab = 0.7, main = paste0("PCO", axis_id, " with spatiotemporal model predictions"))
  
  for (j in 2:5) {
    points(wangenv$dist[wangenv$time == levels(wangenv$time)[j]], pwangbug[wangenv$time == levels(wangenv$time)[j], pco_no],
           col = day_cols[j], pch = day_pch[j])
  }
  
  for (j in 1:5) {
    idx <- xygrid$time == levels(xygrid$time)[j]
    lines(xygrid$dist[idx], pco_predict2[idx, pred_no], type = "l", lty = day_lty[j], col = day_cols[j], lwd = 2)
  }
  
  legend("topright", inset = c(0, 0), legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
         lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2, col = day_cols, ncol = 2, horiz = FALSE, cex = 0.5, title = "months")
}

dev.off()

# ---- 11. Produce residual diagnostics and save the residual figure for the first PCO axis ---------------------
resid_d <- vegan::vegdist(as.data.frame(pco_resid), method = "euclidean")
dist_d <- stats::dist(scale(wangenv[, c("day")], scale = FALSE), method = "euclidean")
if (requireNamespace("ecodist", quietly = TRUE)) {
  tiff(file = fs::path(figure_dir, "03_mantel_correlograms.tif"), width = 8, height = 8, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
  par(mfrow = c(2, 1))
  plot(ecodist::mgram(resid_d, dist_d), xlab = "Lag")
  plot(vegan::mantel.correlog(resid_d, dist_d, nperm = 99, n.class = 16, cutoff = FALSE))
  dev.off()
}
bugenvres <- cbind(bugenv, as.data.frame(pco_resid))
tiff(file = fs::path(figure_dir, "03_residuals_after_modeling_pco.tif"), width = 8, height = 10, units = "in", pointsize = 16, bg = "transparent", res = 800, compression = "lzw")
par(mfrow = c(2, 1), cex.axis = 0.7, las = 1, cex.main = 1)
plot(bugenvres[bugenvres$day == day_values[1], "dist"], bugenvres[bugenvres$day == day_values[1], "res1"], type = "b", ylim = c(-0.8, 0.88), pch = 19, cex = 1.5, xlab = "Distance(km)", ylab = "Residuals", main = "Residual plotted against spatial position", lwd = 3, lty = 1)
for (j in 2:length(day_values)) points(bugenvres[bugenvres$day == day_values[j], "dist"], bugenvres[bugenvres$day == day_values[j], "res1"], type = "b", col = day_cols[j], pch = day_pch[j], lwd = 3, cex = if (j == 5) 2 else 1.5, lty = day_lty[j])
legend("topright", inset = c(0, 0), legend = day_labels, lty = c(1, 5), pch = day_pch, lwd = 2, col = day_cols, ncol = 2, horiz = FALSE, cex = 0.8, title = "months")
abline(v = 4, lty = 2)
plot(mydates, bugenvres[bugenvres$dist == dist_values[1], "res1"], type = "b", ylim = c(-0.8, 0.8), pch = 19, cex = 1.5, xlab = "day", ylab = "Residuals", main = "Residuals plotted against time", lwd = 3, lty = 2, cex.main = 1, xaxt = "n")
if (length(dist_values) > 1) for (j in 2:length(dist_values)) points(mydates, bugenvres[bugenvres$dist == dist_values[j], "res1"], type = "b", col = dist_cols[j], pch = dist_pch[j], lwd = 3, cex = if (j >= 5) 1.5 else 1, lty = dist_lty[j])
axis.Date(side = 1, mydates, at = seq(mydaterange[1], mydaterange[2], by = "month"), format = "%b-%y")
legend("bottomright", inset = c(0, 0), legend = levels(as.factor(bugenv$dist)), lty = 1, lwd = 2, col = dist_cols, ncol = 4, horiz = FALSE, cex = 0.8, title = "distance(km)")
dev.off()

# ---- 12. Save the model objects, predictions, residuals, and plotting tables needed downstream ----------------
saveRDS(
  list(
    wangbug = wangbug,
    wangenv = wangenv,
    wangbug.BC = wangbug_BC,
    pwangbug = pwangbug,
    alldata = alldata,
    envpcores = envpcores,
    bugenvres = bugenvres
  ),
  fs::path(output_dir, "03_environmental_script_inputs.rds")
)