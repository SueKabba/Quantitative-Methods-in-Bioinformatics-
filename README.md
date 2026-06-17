# Quantitative Methods in Bioinformatics: Translational Cancer Therapeutics

An end-to-end computational and biostatistical analysis supporting the engineering and validation of a novel cancer therapeutic workflow. This investigation spans molecular structural genomics, stochastic pharmacokinetic simulation, discrete micro-vector optimization, and longitudinal Phase I clinical trial validation.

---

## Repository Directory Layout & Asset Blueprint
To maintain production grade modularity and clarity for engineering reviews, the computational infrastructure is divided into four standalone, highly commented execution scripts:

```text
Quantitative-Methods-in-Bioinformatics/
├── Project-1-Cancer-Therapeutics/
│   ├── part1_trinucleotide_expansion.R <- Structural genomic analysis & Fisher's Exact Test
│   ├── part2_drug_pharmacokinetics.R   <- Physiological clearance sums & Monte Carlo models
│   ├── part3_liposomal_delivery.R      <- Discrete distributions & stochastic vector analysis
│   └── part4_drug_effectiveness.R      <- Spatial baseline counts & Student's t Clinical Trial CI
└── README.md                           <- Central Portfolio Executive Report (This Document)
```

# Computational Tech Stack
- Language Environment: R(v4.0+)
- Genomic Library: Biostings (Exact string matching and vector alignment)
- Data Core: tidyverse(dplyr, readr, ggplot2)
- Statistical Engine: stats(Base R probability distributions)

# Pipeline Summary & Key Statistical Insights
Phase 1: Trinucleotide Repeat Expansion (TRE) Modeling 
- The Goal: Investigate if abnormal nucleotide elongation (the repetition of the 'CTG' motif) inside the open reading frame of the BNFO gene shares a significant link with clinical oncogenesis.
- The Logic: Unaffected healthy populations show an average of 33 repeats. Setting a strict statistical safety threshold (alpha = 0.01), an upper boundary cutoff was established at 47 repeats. Any donor exceeding 47 repeats is classified as "structurally expanded."
- The Finding: Processing 20 donor profiles with Fisher's Exact Test yielded a p-value of 0.06978. Because this sits above the standard baseline (0.05), we cannot reject the null hypothesis of independence. The sample cohort (n = 20) lacks the necessary statistical power, highlighting a clear need for expanded clinical trials.

Phase 2: Dosing Reliability & Pharmacokinetics
- The Goal: Evaluate wether a standard once daily IV bolus dose (200 ug/mL) effectively sustains a therapeutic blood concentration margin (>= 70 ug/mL) against natural bodily clearance.
- The Logic: Systemic clearance is modeled across three independent, continuous normal degradation channels: Absorption, Excretion, and Metabolism. Summing these channels reveals an average total daily clearance of 134 ug/mL/day with a standard deviation of 6.782 ug/mL/day.
- The Finding: Analytical calculations and a 100,000-iteration Monte Carlo simulation exposed an inadequate 27.8% clinical reliability rating—meaning blood concentrations fail and drop below the target margin 72.2% of the time.
- The Solution: Sensitivity calibration rules prove that moving the initial IV bolus target upward into a 213 to 215 ug/mL range instantly scales clinical stability to a dependable >= 90% safety window.

Phase 3: Liposomal Vector Drug Delivery Optimization
- The Goal: Choose an optimal micro-encapsulation packaging setup to maximize the probability that a target cell successfully absorbs an effective treatment payload of >= 4 drug molecules (internalization probability p = 0.25).
- The Logic: I compared two distinct structural formulations across an identical payload size of 20 active molecules:
  * Method 1 (Single Carrier Vector): 1 large vesicle containing all 20 molecules. This acts as an all-or-nothing Bernoulli trial.
    
  * Method 2 (Distributed Fleet): 20 individual micro-liposomes carrying 1 molecule each, operating as an independent Binomial distribution.
- The Finding: While both approaches yield an identical expected average of 5 molecules, Method 1 creates extreme, volatile variance (Var = 75.0), leading to a low 25% delivery success rate. By spreading risk across independent vehicles, Method 2 slashes variance down to 3.75 and more than tripling clinical reliability to a 77.5% success rate. Method 2 is heavily recommended for development.

Phase 4: Longitudinal Efficacy & Clinical Trial Validation
- The Goal: Untreated tumor cell density follows a spatial Poisson distribution with an expected mean of 52.6 cells/mm². Post-treatment tissue biopsies revealed a remaining cell density of 42 cells.
- The Finding: Untreated tumor cell density follows a spatial Poisson distribution with an expected mean of 52.6 cells/mm². Post-treatment tissue biopsies revealed a remaining cell density of 42 cells.
- The Verification: Processing the full trial cohort (n = 50) using a continuous Student's t-distribution confidence interval framework with 49 degrees of freedom created an explicit decision matrix. A confidence interval landing strictly above zero mathematically proves that the treatment effect is robust, systematic, and fully reproducible.
