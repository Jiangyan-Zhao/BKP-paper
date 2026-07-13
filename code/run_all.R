## -------------------------------------------------------------------------
## Reproduce all analyses for the BKP paper
##
## Run this script from the repository root:
##
##   source("code/run_all.R")
## -------------------------------------------------------------------------

## Create output directories if needed
dir.create("code/figure", recursive = TRUE, showWarnings = FALSE)
dir.create("code/result", recursive = TRUE, showWarnings = FALSE)


## -------------------------------------------------------------------------
## Section 4: Illustrative examples
## -------------------------------------------------------------------------

cat("\n===== Section 4: Example 1 =====\n")
source("code/s4_ex1_bkp_1d_logistic.R")

cat("\n===== Section 4: Example 2 =====\n")
source("code/s4_ex2_bkp_1d_nonlinear.R")

cat("\n===== Section 4: Example 3 =====\n")
source("code/s4_ex3_bkp_2d_goldstein_price.R")

cat("\n===== Section 4: Example 4 =====\n")
source("code/s4_ex4_bkp_two_spirals_classification.R")

cat("\n===== Section 4: Example 5 =====\n")
source("code/s4_ex5_dkp_1d_multinomial.R")

cat("\n===== Section 4: Example 6 =====\n")
source("code/s4_ex6_dkp_2d_multinomial.R")

cat("\n===== Section 4: Example 7 =====\n")
source("code/s4_ex7_dkp_iris_classification.R")

cat("\n===== Section 4: Example 8 =====\n")
source("code/s4_ex8_twinbkp_1d_nonlinear.R")

cat("\n===== Section 4: Example 9 =====\n")
source("code/s4_ex9_twindkp_1d_multinomial.R")


## -------------------------------------------------------------------------
## Section 5: Real-data applications
## -------------------------------------------------------------------------

cat("\n===== Section 5: Loa loa prevalence mapping =====\n")
source("code/s5_app1_loaloa_prevalence_mapping.R")

cat("\n===== Section 5: Mourning Warbler distribution modeling =====\n")
source("code/s5_app2_mourning_warbler_sdm.R")


## -------------------------------------------------------------------------
## Appendix: Coverage simulation
## -------------------------------------------------------------------------

cat("\n===== Appendix: Predictive interval coverage =====\n")
source("code/a_coverage.R")


cat("\n=============================================================\n")
cat("All analyses completed successfully.\n")
cat("Figures are available in code/figure/.\n")
cat("Numerical results are available in code/result/.\n")
cat("=============================================================\n")
