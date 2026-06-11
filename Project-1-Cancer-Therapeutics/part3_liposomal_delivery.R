############################################################
# Part 3 — Liposomal Drug Delivery
# Goal: Compare two delivery methods for hitting ≥4 molecules
# delivered inside a *single* cancer cell.
#
# Method 1: 1 liposome, packed with 20 molecules
# Method 2: 20 liposomes, each with 1 molecule
#
# Each liposome independently reaches and is internalized
# by a cell with probability p = 0.25.
#
# X = # of drug molecules internalized by the cell
# What are E[X], Var[X], PMF, and P(X >= 4) for each method?
############################################################

# -----------------------------
# Parameters
# -----------------------------
p <- 0.25              # success (internalization) probability per liposome
cap <- 20              # molecules per liposome (method 1 capacity)
need <- 4              # clinical minimum needed in cell

# -----------------------------
# 1) Method 1: 1 liposome, 20 molecules
# -----------------------------
# If the liposome is internalized (Bernoulli(p)), all 20 arrive; else 0 arrive.
# So X1 = 20 * Bernoulli(p)
# PMF: P(X1=20)=p; P(X1=0)=1-p
# E[X1] = 20*p; Var[X1] = (20^2)*p*(1-p)
E_X1   <- cap * p
Var_X1 <- (cap^2) * p * (1 - p)

# Probability to achieve at least 4 molecules:
# Since it's either 0 or 20, P(X1 >= 4) = P(X1 = 20) = p
P_X1_ge_need <- p

pmf_method1 <- data.frame(
  molecules = c(0, cap),
  probability = c(1 - p, p)
)

# -----------------------------
# 2) Method 2: 20 liposomes, 1 molecule each
# -----------------------------
# Each liposome independently gets internalized with probability p.
# X2 ~ Binomial(n=20, p)
n2 <- cap   # 20 one-molecule liposomes
E_X2   <- n2 * p
Var_X2 <- n2 * p * (1 - p)

# PMF and tail probability
k_vals <- 0:n2
pmf_method2 <- data.frame(
  molecules = k_vals,
  probability = dbinom(k_vals, size = n2, prob = p)
)

P_X2_ge_need <- 1 - pbinom(q = need - 1, size = n2, prob = p)  # P(X2 >= 4)

# -----------------------------
# 3) Print core results
# -----------------------------
cat("\n=== METHOD 1: 1 liposome with 20 molecules ===\n")
cat(sprintf("E[X1] = %.2f\n", E_X1))
cat(sprintf("Var[X1] = %.2f\n", Var_X1))
cat(sprintf("P(X1 >= %d) = %.3f\n", need, P_X1_ge_need))

cat("\n=== METHOD 2: 20 liposomes with 1 molecule each ===\n")
cat(sprintf("E[X2] = %.2f\n", E_X2))
cat(sprintf("Var[X2] = %.2f\n", Var_X2))
cat(sprintf("P(X2 >= %d) = %.3f\n", need, P_X2_ge_need))

# -----------------------------
# 4) Optional: PMF tables
# -----------------------------
cat("\nPMF (Method 1):\n"); print(pmf_method1)
cat("\nPMF (Method 2) first 10 rows:\n"); print(head(pmf_method2, 10))

# -----------------------------
# 5) Optional: PMF plots (base R)
# -----------------------------
# Barplot for Method 1 (two-point distribution)
barplot(height = pmf_method1$probability, names.arg = pmf_method1$molecules,
        main = "PMF — Method 1 (1 liposome with 20 molecules)",
        xlab = "Molecules delivered", ylab = "Probability")

# Barplot for Method 2 (Binomial)
barplot(height = pmf_method2$probability, names.arg = pmf_method2$molecules,
        main = "PMF — Method 2 (20 liposomes, 1 molecule each)",
        xlab = "Molecules delivered", ylab = "Probability")

