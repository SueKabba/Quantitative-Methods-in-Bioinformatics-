############################################################
# Part 4 — Drug Effectiveness
# 1) Pre-treatment cancer cell counts follow Poisson(λ=52.6)
# Plot PMF of X = # cells per 1 mm^2 image
# 95% central interval for X
# After treatment, observed X_obs = 42 → interpret
#
# 2) Clinical trial (n=50 subjects):
# 'Tumor_Data.csv' contains per-subject reductions (after - before).
# Compute 95% CI for the mean reduction
# Comment on effectiveness (does CI exclude 0?)

# -----------------------------
# 1) Poisson model for pre-treatment counts
# -----------------------------
lambda <- 52.6

# Range for PMF plot: mean ± 4*sd is a good visual window
rng <- c(max(0, floor(lambda - 4*sqrt(lambda))),ceiling(lambda + 4*sqrt(lambda)))
k <- rng[1]:rng[2]
pmf <- dpois(k, lambda)

# PMF plot (base R)
plot(k, pmf, type = "h",main = "PMF of Pre-Treatment Cell Counts: X ~ Poisson(λ = 52.6)",
     xlab = "Cells per 1 mm^2 image", ylab = "Probability")
points(k, pmf, pch = 16)

# 95% central interval for X (as a random variable)
ci_lower <- qpois(0.025, lambda)
ci_upper <- qpois(0.975, lambda)

cat("\n--- PART 4.1: Pre-treatment Poisson model ---\n")
cat(sprintf("95%% central interval for X: [%d, %d]\n", ci_lower, ci_upper))

# After treatment, a new image shows X_obs = 42.
X_obs <- 42
# Deductive statement: Compute tail probability under pre-treatment model.
# P(X <= 42 | λ = 52.6)
p_le_42 <- ppois(q = X_obs, lambda = lambda)

cat(sprintf("P(X <= %d | λ=52.6) = %.4f\n", X_obs, p_le_42))

# -----------------------------
# Interpretation guidance:
# - If P <= 0.05, strong evidence of reduction vs baseline distribution.
# - Here we will print a friendly note for the report text.
if (p_le_42 <= 0.05) {
  cat("Interpretation: This is unlikely under the baseline (<=5%), strong evidence of reduction.\n")
} else if (p_le_42 <= 0.10) {
  cat("Interpretation: Moderately unlikely under baseline (~10%), suggestive evidence of reduction.\n")
} else {
  cat("Interpretation: Plausible under baseline; not strong evidence of reduction by itself.\n")
}

# -----------------------------
# 2) Clinical trial: 50-subject reductions
# -----------------------------
# Expect a CSV with a single numeric column named, e.g., 'Reduction'
# Positive numbers = reduction in cell count; negative = increase.

file_path <- "Tumor_Data.csv"

if (file.exists(file_path)) {
  dat <- read.csv(file_path)
  # Try common column names
  if ("Reduction" %in% names(dat)) {
    y <- dat$Reduction
  } else if ("reduction" %in% names(dat)) {
    y <- dat$reduction
  } else if (ncol(dat) == 1) {
    y <- dat[[1]]
    names(y) <- "Reduction"
  } else {
    stop("Could not find the reduction column. Ensure the CSV has a single 'Reduction' column.")
  }
  
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < 2) stop("Need at least 2 observations for a CI.")
  
  ybar <- mean(y)
  s <- sd(y)
  se <- s / sqrt(n)
  
  # t CI for mean reduction (unknown sigma)
  alpha <- 0.05
  tcrit <- qt(1 - alpha/2, df = n - 1)
  ci_mean <- c(ybar - tcrit * se, ybar + tcrit * se)
  
  cat("\n--- PART 4.2: Clinical trial (n=50) ---\n")
  cat(sprintf("Mean reduction = %.3f\n", ybar))
  cat(sprintf("95%% CI for mean reduction = [%.3f, %.3f]\n", ci_mean[1], ci_mean[2]))
  
  # Simple effectiveness decision rule:
  if (ci_mean[1] > 0) {
    cat("Effectiveness: CI is strictly > 0 → evidence the drug reduces tumor cell density on average.\n")
  } else if (ci_mean[2] < 0) {
    cat("Effectiveness: CI is strictly < 0 → evidence the drug INCREASES tumor cell density (not desired).\n")
  } else {
    cat("Effectiveness: CI includes 0 → inconclusive effect at 95% confidence.\n")
  }
  
} else {
  cat("\n--- PART 4.2: Clinical trial (n=50) ---\n")
  cat("Note: 'Tumor_Data.csv' not found in working directory.\n")
  cat("Place the CSV (one numeric column named 'Reduction') and re-run to compute the 95% CI.\n")
}
