## Pharmaceutics Reliability Calculation
## Goal is to analyze the probability that the blood concentration
## stays >= 70 ug/mL for the entire 24 hr after an IV bolus dose of 200 ug/mL

##  Over 24h, drug is REMOVED by three independent processes:
##  D1 (absorption from bloodstream), D2 (excretion), D3 (metabolism).
##  Units given as μg/mL/day; the period is exactly 1 day (24h).
##  D1 ~ Normal(μ1=37,  σ1^2=12), D2 ~ Normal(μ2=45, σ2^2=18), D3 ~ Normal(μ3=52, σ3^2=16)
##  Initial concentration C0 = 200 μg/mL
##  Effectiveness threshold = 70 μg/mL

##-----------------------
## setting up parameters
##-----------------------
C0 <-200         # initial concentration
threshold <- 70  #effectiveness threshold

mu1 <- 37; var1 <- 12  #abortion form bloodsteam removal
mu2 <- 45; var2 <- 18  #excretion removal 
mu3 <- 52; var3 <- 16  #metabolism removal

##---------------------------------------
##distribution of total removal, D_total
##---------------------------------------
## Sum of independent Nirmals is Normal with:

mu_total <- mu1 + mu2 + mu3
var_total <- var1 + var2 + var3
sd_total <- sqrt(var_total)

##------------------------
##Analytical probability
##------------------------
cutoff <- C0 - threshold
p_effective_analytical <- pnorm(q = cutoff, mean = mu_total, sd = sd_total)

##------------------------
##Monte Carlo Sanity Check
##-------------------------
set.seed(531)
N <- 100000
D1 <- rnorm(N, mean = mu1, sd = sqrt(var1))
D2 <- rnorm(N, mean = mu2, sd = sqrt(var2))
D3 <- rnorm(N, mean = mu3, sd = sqrt(var3))
D_total_sim <- D1 + D2 + D3

p_effective_sim <- mean(D_total_sim <= cutoff)


##------------------------------------------------------
##Sensitivity: What C0 is needed for target reliability?
##-------------------------------------------------------
# Solve for C0 in: P(D_total <= C0 - threshold) >= p_target
# => (C0 - threshold - mu_total) / sd_total >= z_{p_target}
# => C0 >= threshold + mu_total + z_{p_target} * sd_total
targets <- c(0.75, 0.80, 0.90, 0.95)
z_vals  <- qnorm(targets)
C0_needed <- threshold + mu_total + z_vals * sd_total
dose_table <- data.frame(Target_Reliability = targets, z_value = round(z_vals, 3),
                         Required_C0_ug_per_mL = round(C0_needed, 2))


# -----------------------------
# 6) Quick, friendly summaries for your write-up
# -----------------------------
cat("\n--- SUMMARY NUMBERS ---\n")
cat(sprintf("Mean total removal (μ): %.0f μg/mL/day\n", mu_total))
cat(sprintf("Total variance (σ^2): %.0f (μg/mL/day)^2\n", var_total))
cat(sprintf("Total SD (σ): %.3f μg/mL/day\n", sd_total))
cat(sprintf("Cutoff for total removal (C0 - threshold): %.0f μg/mL/day\n", cutoff))
cat(sprintf("Analytical probability of staying >= %.0f all day: %.3f\n", threshold, p_effective_analytical))
cat(sprintf("Simulation probability (N=%d): %.3f\n", N, p_effective_sim))
cat("\nDose needed for reliability targets:\n")
print(dose_table)
