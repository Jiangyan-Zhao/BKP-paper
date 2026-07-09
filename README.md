# BKP-paper

This repository contains the reproduction scripts, data-processing code, generated figures, and manuscript source files for the paper:

**BKP: An R Package for Beta Kernel Process Modeling**

The repository is organized to reproduce the illustrative examples and real-data applications in the manuscript.

## Repository structure

``` text
BKP-paper/
├── code/
│   ├── run_all.R
│   ├── s4_ex1_bkp_1d_logistic.R
│   ├── s4_ex2_bkp_1d_nonlinear.R
│   ├── s4_ex3_bkp_2d_goldstein_price.R
│   ├── s4_ex4_bkp_two_spirals_classification.R
│   ├── s4_ex5_dkp_1d_multinomial.R
│   ├── s4_ex6_dkp_2d_multinomial.R
│   ├── s4_ex7_dkp_iris_classification.R
│   ├── s4_ex8_twinbkp_1d_nonlinear.R
│   ├── s4_ex9_twindkp_1d_multinomial.R
│   ├── s5_app1_loaloa_prevalence_mapping.R
│   ├── s5_app2_mourning_warbler_sdm.R
│   ├── data/
│   ├── figure/
│   └── result/
├── paper/
│   ├── TR_BKP.tex
│   ├── TR_BKP.pdf
│   ├── refs.bib
│   ├── jss.cls
│   └── jss.bst
├── renv.lock
└── README.md
```

## Requirements

The analyses were run in R using the package environment recorded in `renv.lock`.

To restore the package environment, open the repository as an R project and run:

``` r
renv::restore()
```

The main package used in the paper is `BKP`. The reproduction scripts also use several supporting packages, including `tgp`, `gplite`, `kernlab`, `pROC`, `mlbench`, `ggplot2`, `gridExtra`, `RiskMap`, `sf`, `terra`, `maps`, `rnaturalearth`, and related dependencies.

## Reproducing the examples

To reproduce all examples and real-data applications, run:

``` r
source("code/run_all.R")
```

The script `code/run_all.R` executes the Section 4 simulation examples and the Section 5 real-data applications.

Section 4 includes:

``` text
s4_ex1_bkp_1d_logistic.R
s4_ex2_bkp_1d_nonlinear.R
s4_ex3_bkp_2d_goldstein_price.R
s4_ex4_bkp_two_spirals_classification.R
s4_ex5_dkp_1d_multinomial.R
s4_ex6_dkp_2d_multinomial.R
s4_ex7_dkp_iris_classification.R
s4_ex8_twinbkp_1d_nonlinear.R
s4_ex9_twindkp_1d_multinomial.R
```

Section 5 includes:

``` text
s5_app1_loaloa_prevalence_mapping.R
s5_app2_mourning_warbler_sdm.R
```

Generated figures are saved to:

``` text
code/figure/
```

Precomputed or intermediate numerical results are saved to:

``` text
code/result/
```

## Timing experiment

The timing comparison in Example 3 can be computationally expensive. By default, the script reads the precomputed average timing results from:

``` text
code/result/elapsed_time_avg.csv
```

To rerun the full timing experiment, edit the following line in `code/s4_ex3_bkp_2d_goldstein_price.R`:

``` r
run_elapsed_time <- TRUE
```

The full timing experiment may take substantial time, especially for the optimized logistic Gaussian process benchmark.

## Real-data applications

The Loa loa prevalence mapping example uses the `loaloa` dataset from the `RiskMap` package.

The Mourning Warbler species distribution modeling example uses the data and climate raster files stored under:

``` text
code/data/
```

This example compares BKP, TwinBKP, and a logistic Gaussian process model. It may take longer to run because it includes model fitting, ROC evaluation, and raster-level spatial prediction.

## Manuscript

The LaTeX source files are in:

``` text
paper/
```

The main manuscript file is:

``` text
paper/TR_BKP.tex
```

The manuscript is configured to read generated figures from:

``` text
code/figure/
```

via the LaTeX graphic path setting in `TR_BKP.tex`.

## Notes

All scripts are written to be run from the repository root. For example:

``` r
source("code/s4_ex1_bkp_1d_logistic.R")
```

not from inside the `code/` directory.

For a clean reproduction, avoid relying on saved R workspaces such as `.RData`. These files are ignored by `.gitignore`.
