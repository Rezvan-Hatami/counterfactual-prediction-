# 11_nitrate_models_and_figures.R
# Author: Rezvan Hatami
# Date: 16 March 26
#
# Purpose:
# This script reproduces the nitrate modelling workflow in a self-contained, and reproducible form. 

rm(list = ls())

# ---- Setup: package management ------------------------------------------------
required_pkgs <- c("fs", "here", "lattice", "latticeExtra", "car")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "These packages must be installed before running the script: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(scipen = 1)
set.seed(123)

# ---- Setup: path configuration ------------------------------------------------
project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
figure_dir <- fs::path(project_dir, "figures", "11")

fs::dir_create(output_dir)
fs::dir_create(figure_dir)

# ---- Setup: read required inputs ----------------------------------------------
wangenv_file <- fs::path(output_dir, "01_wangenv_prepped.rds")
bugenv_file <- fs::path(output_dir, "02_bugenv_with_pco.rds")
script3_inputs_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")

if (!fs::file_exists(wangenv_file)) {
  stop("Missing input file: ", wangenv_file)
}
if (!fs::file_exists(bugenv_file)) {
  stop("Missing input file: ", bugenv_file)
}
if (!fs::file_exists(script3_inputs_file)) {
  stop("Missing input file: ", script3_inputs_file)
}

wangenv <- readRDS(wangenv_file)
bugenv <- readRDS(bugenv_file)

script3_inputs <- readRDS(script3_inputs_file)
if (!is.list(script3_inputs) || !"alldata" %in% names(script3_inputs)) {
  stop("03_environmental_script_inputs.rds does not contain alldata.")
}
alldata <- as.data.frame(script3_inputs$alldata)

# ---- Setup: helper functions --------------------------------------------------
panel_cor <- function(x, y, digits = 2, ...) {
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(invisible(NULL))
  }
  
  if (sd(x) == 0 || sd(y) == 0) {
    text(0.5, 0.5, "sd = 0")
    return(invisible(NULL))
  }
  
  r <- suppressWarnings(cor(x, y, method = "spearman"))
  txt <- paste0("r= ", formatC(r, digits = digits, format = "f"))
  text(0.5, 0.8, txt)
  
  p <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  txt2 <- if (is.na(p)) {
    "p= NA"
  } else if (p < 0.01) {
    "p= <0.01"
  } else {
    paste0("p= ", formatC(p, digits = digits, format = "f"))
  }
  text(0.5, 0.3, txt2)
  
  invisible(NULL)
}

num_or_na <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

log_pos <- function(x) {
  x <- num_or_na(x)
  x[!is.finite(x) | x <= 0] <- NA_real_
  log(x)
}

# -------------------------------------------------------------------
# 1. Scatterplot checks for nitrate
# -------------------------------------------------------------------
pairs_df <- data.frame(
  ph2 = num_or_na(alldata$ph2),
  log_alk = log_pos(alldata$alk),
  log_dno3_dayflow = log_pos(num_or_na(alldata$dno3) / num_or_na(alldata$dayflow)),
  eff = num_or_na(alldata$eff),
  dayflow = num_or_na(alldata$dayflow),
  dist = num_or_na(alldata$dist),
  log_temp = log_pos(alldata$temp),
  log_toc = log_pos(alldata$toc),
  log_no3 = log_pos(alldata$no3),
  no3 = num_or_na(alldata$no3)
)

pairs_df <- pairs_df[, colSums(is.finite(as.matrix(pairs_df))) > 0, drop = FALSE]
pairs_df <- pairs_df[stats::complete.cases(pairs_df), , drop = FALSE]

if (ncol(pairs_df) < 2) {
  stop("Not enough finite columns available for nitrate scatterplot checks.")
}
if (nrow(pairs_df) < 3) {
  stop("Not enough complete rows available for nitrate scatterplot checks.")
}

par(mar = c(5, 5, 4, 2), cex = 0.9)

pairs(
  pairs_df,
  upper.panel = panel_cor,
  pch = 20
)
# -------------------------------------------------------------------
# 2. Nitrate regression with other variables
# -------------------------------------------------------------------
nitrate.model1 <- lm(no3 ~ eff, data = alldata)
nitrate.model2 <- lm(no3 ~ dno3, data = alldata)
nitrate.model3 <- lm(no3 ~ dayflow, data = alldata)
nitrate.model4 <- lm(no3 ~ dist, data = alldata)
nitrate.model5 <- lm(no3 ~ temp, data = alldata)

print(summary(nitrate.model1))
print(summary(nitrate.model2))
print(summary(nitrate.model3))
print(summary(nitrate.model4))
print(summary(nitrate.model5))

# -------------------------------------------------------------------
# 3. Nitrate plotted against spatial position
# -------------------------------------------------------------------
plot(
  bugenv[bugenv$day == 1, "dist"],
  bugenv[bugenv$day == 1, "no3"],
  type = "b", ylim = c(0, 0.7), pch = 19,
  xlab = "Distance(km)", ylab = "Nitrate (mg/L)",
  main = "Nitrate plotted against spatial position",
  lwd = 3, lty = 1
)
points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "no3"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "no3"], type = "b", col = "blue", pch = 8, lwd = 3, cex = 1.2, lty = 3)
points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "no3"], type = "b", col = "green3", pch = 17, lwd = 3, cex = 1.2, lty = 4)
points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "no3"], type = "b", col = "red", pch = 18, lwd = 3, cex = 1.4, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

# -------------------------------------------------------------------
# 4. Multiple regression analysis for nitrate vs time, distance and eff
# -------------------------------------------------------------------
tempo <- lm(log(no3) ~ time + eff:time + eff:dist + eff:dist:time + eff:I(dist^2):time, data = alldata)
print(anova(tempo))
print(summary(tempo))

tempo1 <- cbind(cbind(exp(fitted(tempo)), exp(residuals(tempo))), alldata)
print(attributes(tempo1))
names(tempo1)[1:2] <- c("fitted", "residual")

plot(
  tempo1$dist[tempo1$time == 1],
  tempo1[tempo1$time == 1, "no3"],
  type = "p", pch = 19, col = "black",
  xlab = "Distance (km)", ylab = "Nitrate (mg/L)",
  main = "Nitrate with spatiotemporal prediction",
  ylim = c(0, 0.7), cex = 1, lwd = 3, lty = 1
)
points(tempo1$dist[tempo1$time == 2], tempo1[tempo1$time == 2, "no3"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(tempo1$dist[tempo1$time == 3], tempo1[tempo1$time == 3, "no3"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1$dist[tempo1$time == 4], tempo1[tempo1$time == 4, "no3"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1$dist[tempo1$time == 5], tempo1[tempo1$time == 5, "no3"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
points(tempo1[tempo1$day == 1,   "dist"], tempo1[tempo1$day == 1,   "fitted"], type = "l", col = "black", pch = 19, cex = 1.5, lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "fitted"], type = "l", col = "gold2", pch = 15, cex = 1.5, lwd = 3, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "fitted"], type = "l", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "fitted"], type = "l", col = "green3", pch = 17, cex = 1.5, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "fitted"], type = "l", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

plot(
  tempo1[tempo1$day == 1, "dist"],
  tempo1[tempo1$day == 1, "residual"],
  type = "b", ylim = c(0, 2), pch = 19, col = "black",
  cex = 1.5, xlab = "Distance(km)", ylab = "Nitrate",
  main = "residuals", lwd = 3, lty = 1
)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "residual"], type = "l", col = "gold2", pch = 15, cex = 1.5, lwd = 3, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "residual"], type = "l", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "residual"], type = "l", col = "green3", pch = 17, cex = 1.5, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "residual"], type = "l", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

plot(tempo)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
hist(tempo$residuals)

# -------------------------------------------------------------------
# 5. Multiple regression analysis for no3 and variables affecting it
# -------------------------------------------------------------------
#tempo <- lm(no3 ~ I(time == 3) + eff:dist + eff:I((dno3)/dayflow) + eff:I((dno3)/dayflow):dist + eff:temp2:dist, data = alldata)
#tempo <- lm(log(no3) ~ time + eff:dist + eff:I(log(dno3)/dayflow)), data = alldata)
#tempo <- lm(log(no3) ~ time + eff:I(log(dno3)/dayflow):time + eff:I(log(dno3)/dayflow):time:dist, data = alldata)
#tempo <- lm(log(no3) ~ temp + toc + time + eff:I(log(dno3)/dayflow):time, data = alldata)
#tempo <- lm(log(no3) ~ temp + time + eff:dist + eff:dist:time + eff:I(dist^2):time + eff:I(log(dno3/dayflow)) + eff:I(log(dno3/dayflow)):time, data = alldata)
tempo <- lm(log(no3) ~ temp + time + alk:time + eff:I(log(dno3)/dayflow) + eff:I(log(dno3)/dayflow):time, data = alldata)
tempo <- lm(log(no3) ~ temp + rain2 + rain3 + alk:time + eff:I(log(dno3)/dayflow) + eff:I(log(dno3)/dayflow):time, data = alldata)
tempo <- lm(log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3)/dayflow):time + eff:I(log(dno3)/dayflow):I(dist == 4.08):time, data = alldata)
print(anova(tempo))
print(summary(tempo))

tempo1 <- cbind(cbind(exp(fitted(tempo)), exp(residuals(tempo))), alldata)
names(tempo1)[1:2] <- c("fitted", "residual")

plot(
  tempo1[tempo1$day == 1, "dist"],
  tempo1[tempo1$day == 1, "fitted"],
  type = "l", ylim = c(0, 0.7), pch = 19, col = "black",
  xlab = "Distance(km)", ylab = "Nitrate (mg/L)",
  main = "Prediction with environmental variables",
  lwd = 3, lty = 1
)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "fitted"], type = "l", col = "gold2", pch = 15, lwd = 3, cex = 1.5, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "fitted"], type = "l", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "fitted"], type = "l", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "fitted"], type = "l", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

avPlots(tempo)
plot(tempo)
qqnorm(tempo$residuals)
qqline(tempo$residuals)
hist(tempo$residuals)

nitrate.check <- lm(residuals(tempo) ~ time + eff:time + eff:dist + eff:dist:time + eff:I(dist^2):time, data = alldata)
print(summary(nitrate.check))
print(anova(nitrate.check))

# -------------------------------------------------------------------
# 6. Some more graphs
# -------------------------------------------------------------------
a <- xyplot(no3 ~ dist, group = time, data = alldata, ylim = c(0, 0.8), type = "p", pch = 19, lty = 2, col = c("black", "gold2", "blue", "green3", "red"))
b <- xyplot(exp(fitted(tempo)) ~ dist, group = time, data = alldata, type = "l", lwd = 2, col = c("black", "gold2", "blue", "green3", "red"))
c <- xyplot(toc ~ dist, data = alldata, panel = function(x, y) { panel.abline(v = 4, lty = 2) })
print(a + as.layer(b) + as.layer(c))

xyplot(residuals(tempo) ~ dist, groups = time, data = alldata, type = "b", pch = 19, lwd = 2, col = c("black", "gold2", "blue", "green3", "red"))
xyplot(no3 ~ dist, group = time, data = alldata, type = "b", pch = 19, lty = 2, col = c("black", "gold2", "blue", "green3", "red"))
xyplot(nh3 ~ dist, group = time, data = alldata, type = "b", pch = 19, lty = 2, col = c("black", "gold2", "blue", "green3", "red"))
xyplot(do ~ dist, group = time, data = alldata, type = "b", pch = 19, lty = 2, col = c("black", "gold2", "blue", "green3", "red"))
xyplot(log(no3) ~ toc, groups = time, data = alldata, type = "p", pch = 19:24, lwd = 2, col = c("black", "gold2", "blue", "green3", "red"))
xyplot(do ~ toc, groups = time, data = alldata, type = "p", pch = 19:24, lwd = 2, col = c("black", "gold2", "blue", "green3", "red"))
xyplot(log(no3) ~ log(nh3), groups = time, data = alldata, type = "p", pch = 19:24, lwd = 2, col = c("black", "gold2", "blue", "green3", "red"))
xyplot(log(nh3) ~ log(toc), groups = time, data = alldata, type = "p", pch = 19:24, lwd = 2, col = c("black", "gold2", "blue", "green3", "red"))

print(cor.test(alldata$do, alldata$toc))

# -------------------------------------------------------------------
# 7. Final nitrate figure
# -------------------------------------------------------------------
tiff(file = fs::path(figure_dir, "11_nitrate_plot.tif"),
     width = 10, height = 14, units = "in",
     pointsize = 12, bg = "transparent",
     res = 800, compression = "lzw")

par(mfrow = c(3, 1), mar = c(3, 4.5, 2.5, 3.5), cex = 1.5, cex.axis = 0.9, las = 1, cex.main = 1, cex.lab = 0.8)

plot(bugenv[bugenv$day == 1, "dist"], bugenv[bugenv$day == 1, "no3"], type = "b", ylim = c(0, 0.7), pch = 19,
     xlab = "Distance(km)", ylab = "Nitrate (mg/L)", main = "a) Nitrate plotted against spatial position", lwd = 3, lty = 1)
points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "no3"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.2, lty = 2)
points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "no3"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "no3"], type = "b", col = "green3", pch = 17, lwd = 3, cex = 1.2, lty = 4)
points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "no3"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

tempo <- lm(log(no3) ~ time + eff:time + eff:dist + eff:dist:time + eff:I(dist^2):time, data = alldata)
tempo1 <- cbind(cbind(exp(fitted(tempo)), exp(residuals(tempo))), alldata)
names(tempo1)[1:2] <- c("fitted", "residual")

plot(tempo1$dist[tempo1$time == 1], tempo1[tempo1$time == 1, "no3"], type = "p", pch = 19, col = "black",
     xlab = "Distance (km)", ylab = "Nitrate (mg/L)", main = "b) Nitrate with spatiotemporal prediction",
     ylim = c(0, 0.7), cex = 1, lwd = 3, lty = 1)
points(tempo1$dist[tempo1$time == 2], tempo1[tempo1$time == 2, "no3"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(tempo1$dist[tempo1$time == 3], tempo1[tempo1$time == 3, "no3"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1$dist[tempo1$time == 4], tempo1[tempo1$time == 4, "no3"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1$dist[tempo1$time == 5], tempo1[tempo1$time == 5, "no3"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
points(tempo1[tempo1$day == 1,   "dist"], tempo1[tempo1$day == 1,   "fitted"], type = "l", col = "black", pch = 19, cex = 1.5, lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "fitted"], type = "l", col = "gold2", pch = 15, cex = 1.5, lwd = 3, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "fitted"], type = "l", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "fitted"], type = "l", col = "green3", pch = 17, cex = 1.5, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "fitted"], type = "l", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

tempo <- lm(log(no3) ~ temp + rain2 + rain3 + eff:I(log(dno3)/dayflow):time + eff:I(log(dno3)/dayflow):I(dist == 4.08):time, data = alldata)
tempo1 <- cbind(cbind(exp(fitted(tempo)), exp(residuals(tempo))), alldata)
names(tempo1)[1:2] <- c("fitted", "residual")

plot(tempo1[tempo1$day == 1, "dist"], tempo1[tempo1$day == 1, "fitted"], type = "b", ylim = c(0, 0.7), pch = 19,
     col = "black", xlab = "Distance(km)", ylab = "Nitrate (mg/L)",
     main = "c) Prediction with environmental variables", lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "fitted"], type = "b", col = "gold2", pch = 15, lwd = 3, cex = 1.5, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "fitted"], type = "b", col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "fitted"], type = "b", col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "fitted"], type = "b", col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

dev.off()

# -------------------------------------------------------------------
# 8. Nonlinear modelling for no3
# -------------------------------------------------------------------
temp3 <- alldata$temp2 - 20
print(temp3)

alldata2 <- alldata
alldata2$eff_num <- suppressWarnings(as.numeric(as.character(alldata2$eff)))
if (any(is.na(alldata2$eff_num))) alldata2$eff_num <- as.numeric(alldata2$eff) - 1

no3nlm <- 0.01 + 0.33 * (alldata2$time == 3) +
  1.8 * alldata2$eff_num * (alldata2$dno3 / alldata2$dayflow) *
  exp(-0.01 * 1.1^alldata2$temp2 * alldata2$dist / alldata2$vel2)

alldata2 <- cbind(alldata2, no3nlm, temp3)

par(mar = c(5, 5, 4, 2), cex = 0.9)
plot(alldata2$dist[alldata2$time == 1], alldata2[alldata2$time == 1, "no3"], type = "p", pch = 19, col = "black",
     xlab = "Distance (km)", ylab = "Nitrate (mg/L)", main = "Nitrate with spatiotemporal prediction",
     ylim = c(0, 0.7), cex = 1, lwd = 3, lty = 1)
points(alldata2$dist[alldata2$time == 2], alldata2[alldata2$time == 2, "no3"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(alldata2$dist[alldata2$time == 3], alldata2[alldata2$time == 3, "no3"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(alldata2$dist[alldata2$time == 4], alldata2[alldata2$time == 4, "no3"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(alldata2$dist[alldata2$time == 5], alldata2[alldata2$time == 5, "no3"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
points(alldata2[alldata2$day == 1,   "dist"], alldata2[alldata2$day == 1,   "no3nlm"], type = "l", col = "black", pch = 15, lwd = 3, cex = 1.5, lty = 1)
points(alldata2[alldata2$day == 126, "dist"], alldata2[alldata2$day == 126, "no3nlm"], type = "l", col = "gold2", pch = 15, lwd = 3, cex = 1.5, lty = 2)
points(alldata2[alldata2$day == 260, "dist"], alldata2[alldata2$day == 260, "no3nlm"], type = "l", col = "blue", pch = 8, lwd = 3, cex = 1.5, lty = 3)
points(alldata2[alldata2$day == 336, "dist"], alldata2[alldata2$day == 336, "no3nlm"], type = "l", col = "green3", pch = 17, lwd = 3, cex = 1.5, lty = 4)
points(alldata2[alldata2$day == 518, "dist"], alldata2[alldata2$day == 518, "no3nlm"], type = "l", col = "red", pch = 18, lwd = 3, cex = 2, lty = 5)
legend("topright", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.9, title = "months")
abline(v = 4, lty = 2)

no3nls <- c(b1 = 0.01, b2 = 0.33, b3 = 1.8, b4 = -0.02, b5 = 1.05)

dat_nls <- subset(
  alldata2,
  is.finite(no3) &
    is.finite(time) &
    is.finite(eff_num) &
    is.finite(dno3) &
    is.finite(dayflow) &
    is.finite(temp3) &
    is.finite(dist) &
    is.finite(vel2) &
    dayflow != 0 &
    vel2 != 0
)

no3nls.fm1 <- nls(
  no3 ~ b1 + b2 * (time == 3) +
    b3 * eff_num * (dno3/dayflow) * exp(b4 * b5^temp3 * dist/vel2),
  data = dat_nls,
  start = no3nls,
  trace = TRUE
)

print(summary(no3nls.fm1))
plot(no3nls.fm1)

fm2 <- as.matrix(fitted(no3nls.fm1))
tempo1 <- cbind(fm2, dat_nls)
names(tempo1)[1] <- "fm2"

plot(tempo1$dist[tempo1$time == 1], tempo1[tempo1$time == 1, "no3"], type = "p", pch = 19, col = "black",
     xlab = "Distance (km)", ylab = "Nitrate (mg/L)", main = "Nitrate with spatiotemporal prediction",
     ylim = c(0, 0.7), cex = 1, lwd = 3, lty = 1)
points(tempo1$dist[tempo1$time == 2], tempo1[tempo1$time == 2, "no3"], col = "gold2", pch = 15, cex = 1.2, lwd = 3, lty = 2)
points(tempo1$dist[tempo1$time == 3], tempo1[tempo1$time == 3, "no3"], col = "blue", pch = 8, cex = 1.2, lwd = 3, lty = 3)
points(tempo1$dist[tempo1$time == 4], tempo1[tempo1$time == 4, "no3"], col = "green3", pch = 17, cex = 1.2, lwd = 3, lty = 4)
points(tempo1$dist[tempo1$time == 5], tempo1[tempo1$time == 5, "no3"], col = "red", pch = 18, cex = 1.4, lwd = 3, lty = 5)
points(tempo1[tempo1$day == 1,   "dist"], tempo1[tempo1$day == 1,   "fm2"], type = "l", col = "black", pch = 19, cex = 1.5, lwd = 3, lty = 1)
points(tempo1[tempo1$day == 126, "dist"], tempo1[tempo1$day == 126, "fm2"], type = "l", col = "gold2", pch = 15, cex = 1.5, lwd = 3, lty = 2)
points(tempo1[tempo1$day == 260, "dist"], tempo1[tempo1$day == 260, "fm2"], type = "l", col = "blue", pch = 8, cex = 1.5, lwd = 3, lty = 3)
points(tempo1[tempo1$day == 336, "dist"], tempo1[tempo1$day == 336, "fm2"], type = "l", col = "green3", pch = 17, cex = 1.5, lwd = 3, lty = 4)
points(tempo1[tempo1$day == 518, "dist"], tempo1[tempo1$day == 518, "fm2"], type = "l", col = "red", pch = 18, cex = 2, lwd = 3, lty = 5)
legend("topleft", inset = c(0, 0), legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
       lty = c(1, 5), pch = c(19, 15, 8, 17, 18), lwd = 2,
       col = c("black", "gold2", "blue", "green3", "red"),
       ncol = 2, horiz = FALSE, cex = 0.6, title = "months")
abline(v = 4, lty = 2)

if (requireNamespace("caret", quietly = TRUE)) {
  lmfit <- caret::train(
    log(no3) ~ temp + rain2 + rain3 +
      eff:I(log(dno3)/dayflow):time +
      eff:I(log(dno3)/dayflow):I(dist == 4.08):time,
    data = alldata
  )
  print(lmfit)
} else {
  cat("\nSkipping caret::train() because package 'caret' is not installed.\n")
}

cat("Script 11 completed successfully.\n")