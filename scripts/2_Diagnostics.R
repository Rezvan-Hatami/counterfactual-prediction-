# 02_pco_diagnostics.R
# Author: Rezvan Hatami
# Date: Oct 2016
# Revised date: 02 March 2026

# ======================================================================
# Bray-Curtis dissimilarity, principal coordinates analysis,
# and ordination diagnostics
# Project root: Causality_Wang
# Data folder: Causality_Wang/data
# Output folders: Causality_Wang/output and Causality_Wang/figures/02
# ======================================================================

# ======================================================================
# Workflow in this script:
# 1. Load the packages required for ordination, diagnostics, and export.
# 2. Define output and figure folders for this script.
# 3. Read the prepared datasets created in 01_data_prep.R.
# 4. Compute Bray-Curtis dissimilarity from square-root transformed
#    abundance data.
# 5. Run principal coordinates analysis with correction for negative
#    eigenvalues.
# 6. Extract eigenvalues, variance explained, cumulative variance, and
#    principal coordinate scores.
# 7. Evaluate ordination structure using broken-stick expectations,
#    bootstrap eigenvalue intervals, and permutation-based diagnostics.
# 8. Combine environmental data and PCO scores for downstream scripts.
# 9. Save ordination objects, diagnostic summaries, and Figure 3.
# ======================================================================

# ======================================================================
# 1. Setup
rm(list = ls())
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(download.file.method = "wininet")
pkgs <- c("vegan", "readr", "dplyr", "tibble", "fs", "here")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) install.packages(miss, dependencies = TRUE, method = "wininet")
invisible(lapply(pkgs, library, character.only = TRUE))

project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
figure_dir <- fs::path(project_dir, "figures", "02")
fs::dir_create(output_dir)
fs::dir_create(figure_dir)

# 2. Read inputs from Script 01
wangbug_file <- fs::path(output_dir, "01_wangbug_raw.rds")
wangenv_file <- fs::path(output_dir, "01_wangenv_prepped.rds")
if (!file.exists(wangbug_file)) stop("Missing input file: ", wangbug_file)
if (!file.exists(wangenv_file)) stop("Missing input file: ", wangenv_file)
wangbug <- readRDS(wangbug_file)
wangenv <- readRDS(wangenv_file)

# 3. Compute Bray-Curtis dissimilarity and PCO
Y <- sqrt(as.matrix(wangbug))
wangbug_BC <- vegan::vegdist(Y, method = "bray")
n <- nrow(wangbug)
p <- n - 1
nmax <- 20
wangbug_mds <- cmdscale(wangbug_BC, k = p, eig = TRUE, add = TRUE, x.ret = FALSE)
wangeigens <- wangbug_mds$eig
pwangbug <- as.data.frame(wangbug_mds$points)
colnames(pwangbug) <- paste0("pco", seq_len(ncol(pwangbug)))
pco_varpercent <- round(wangeigens / sum(wangeigens) * 100, digits = 1)
cumulative_varpercent <- round(cumsum(wangeigens) / sum(wangeigens) * 100, digits = 2)

# 4. Compute diagnostics for Figure 3
D <- as.matrix(vegan::vegdist(Y, method = "bray", diag = TRUE, upper = TRUE))
low <- function(x) quantile(x, probs = 0.025)
high <- function(x) quantile(x, probs = 0.975)

broken_stick <- rep(0, n)
for (k in 1:p) broken_stick[k] <- sum(1 / (k:p))
broken_stick_percent <- 100 * broken_stick / sum(broken_stick)

nboot <- 999
wangeigens_boot <- matrix(0, nrow = nboot, ncol = n)
wangeigens_perc_boot <- matrix(0, nrow = nboot, ncol = n)
for (iboot in 1:nboot) {
  index <- sample(seq_len(n), replace = TRUE)
  D_boot <- D[index, index]
  eig_boot <- cmdscale(D_boot, k = p, eig = TRUE, add = TRUE, x.ret = FALSE)$eig
  wangeigens_boot[iboot, ] <- eig_boot
  wangeigens_perc_boot[iboot, ] <- 100 * eig_boot / sum(eig_boot)
}
lower <- apply(wangeigens_perc_boot, 2, low)
upper <- apply(wangeigens_perc_boot, 2, high)
centre <- apply(wangeigens_perc_boot, 2, mean)

nperm <- 999
wangeigens_perm <- matrix(0, nrow = nperm, ncol = n)
wangeigens_perc_perm1 <- matrix(0, nrow = nperm, ncol = n)
wangeigens_perc_perm2 <- matrix(0, nrow = nperm, ncol = n)
nvars <- ncol(Y)
for (iperm in 1:nperm) {
  Y_perm <- matrix(0, nrow = n, ncol = nvars)
  for (k in 1:nvars) Y_perm[, k] <- Y[sample(seq_len(n), replace = FALSE), k]
  D_perm <- as.matrix(vegan::vegdist(Y_perm, method = "bray", diag = TRUE, upper = FALSE))
  eig_perm <- cmdscale(D_perm, k = p, eig = TRUE, add = TRUE, x.ret = FALSE)$eig
  wangeigens_perm[iperm, ] <- eig_perm
  wangeigens_perc_perm1[iperm, ] <- 100 * eig_perm / sum(eig_perm)
  for (k in 1:n) wangeigens_perc_perm2[iperm, k] <- 100 * eig_perm[k] / sum(eig_perm[k:n])
}
lower_p1 <- apply(wangeigens_perc_perm1, 2, low)
upper_p1 <- apply(wangeigens_perc_perm1, 2, high)
middle_p1 <- apply(wangeigens_perc_perm1, 2, median)
lower_p2 <- apply(wangeigens_perc_perm2, 2, low)
upper_p2 <- apply(wangeigens_perc_perm2, 2, high)
middle_p2 <- apply(wangeigens_perc_perm2, 2, median)

real <- rep(0, n)
for (k in 1:n) real[k] <- 100 * wangeigens[k] / sum(wangeigens[k:n])

# 5. Combine environmental data and PCO scores
bugenv <- cbind(wangenv, pwangbug)

# 6. Save required outputs only
saveRDS(wangbug_BC, fs::path(output_dir, "02_wangbug_BC.rds"))
saveRDS(pwangbug, fs::path(output_dir, "02_pwangbug_scores.rds"))
saveRDS(bugenv, fs::path(output_dir, "02_bugenv_with_pco.rds"))
saveRDS(
  list(
    eigenvalues = wangeigens,
    pco_varpercent = pco_varpercent,
    cumulative_varpercent = cumulative_varpercent,
    broken_stick_percent = broken_stick_percent,
    bootstrap_mean = centre,
    bootstrap_lower = lower,
    bootstrap_upper = upper,
    permutation1_median = middle_p1,
    permutation1_lower = lower_p1,
    permutation1_upper = upper_p1,
    permutation2_median = middle_p2,
    permutation2_lower = lower_p2,
    permutation2_upper = upper_p2,
    real = real,
    nmax = nmax
  ),
  fs::path(output_dir, "02_pco_diagnostics.rds")
)

# 7. Export Figure 3
jpeg(fs::path(figure_dir, "02_Figure_3.jpg"), width = 6, height = 6, units = "in", pointsize = 12, bg = "transparent", res = 800)
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 2) + 0.1, cex = 0.6, cex.axis = 0.8, las = 1)

plot(1:n, 100 * wangeigens / sum(wangeigens), type = "b", xlab = "PCO axis number", ylab = "Percent variation explained", pch = 19, xlim = c(0, nmax))
abline(h = 0)
points(1:n, broken_stick_percent, type = "b", lty = "dotted", cex = 1.3)
mtext("a", 3, -5, cex = 1.8)

plot(1:n, centre, type = "b", xlab = "PCO axis number", ylab = "Percent variation explained", pch = 19, ylim = c(0, max(upper, na.rm = TRUE)), xlim = c(0, nmax))
abline(h = 0)
arrows(1:n, lower, 1:n, upper, code = 3, length = 0.1, angle = 90)
mtext("b", 3, -5, cex = 1.8)

plot(1:n, 100 * wangeigens / sum(wangeigens), type = "b", xlab = "PCO axis number", ylab = "Percent variation explained", pch = 19, ylim = c(0, max(110 * wangeigens / sum(wangeigens), na.rm = TRUE)), xlim = c(0, nmax))
abline(h = 0)
points(1:n, middle_p1, type = "b", lty = "dotted", cex = 1.3)
arrows(1:n, lower_p1, 1:n, upper_p1, code = 3, length = 0.1, angle = 90)
mtext("c", 3, -5, cex = 1.8)

plot(1:(n - 2), real[1:(n - 2)], type = "b", xlab = "PCO axis number", ylab = "Percent variation explained", pch = 19, ylim = c(0, max(110 * wangeigens / sum(wangeigens), na.rm = TRUE)), xlim = c(0, nmax))
abline(h = 0)
points(1:(n - 2), middle_p2[1:(n - 2)], type = "b", lty = "dotted", cex = 1.3)
idx1 <- which((upper_p2 - lower_p2) != 0 & is.finite(lower_p2) & is.finite(upper_p2))
idx1 <- idx1[idx1 <= (n - 2)]
idx0 <- which((upper_p2 - lower_p2) == 0 & is.finite(lower_p2) & is.finite(upper_p2))
idx0 <- idx0[idx0 <= (n - 2)]
arrows(idx1, lower_p2[idx1], idx1, upper_p2[idx1], code = 3, length = 0.1, angle = 90)
points(idx0, lower_p2[idx0], pch = 19)
mtext("d", 3, -5, cex = 1.8)

dev.off()