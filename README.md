# BKP: An R Package for Beta Kernel Process Modeling

This repository contains the reproducibility materials, manuscript source
files, and presentation slides for the BKP software paper.

It includes the analysis scripts, data-processing code, generated figures,
numerical results, manuscript files, and Beamer presentation materials.

The repository is organized to reproduce:

1. the illustrative examples in Section 4;
2. the real-data applications in Section 5;
3. the predictive-coverage simulation reported in the appendix.

A compiled version of the manuscript is available at
[`paper/TR_BKP.pdf`](paper/TR_BKP.pdf), and the presentation slides are
available at [`slides/BKP_Slides.pdf`](slides/BKP_Slides.pdf).

## Project links

- **Interactive website:** [BKP project website](https://jiangyan-zhao.github.io/BKP-website/)
- **R package:** [BKP on CRAN](https://cran.r-project.org/web/packages/BKP/BKP.html)
- **Package development:** [Jiangyan-Zhao/BKP](https://github.com/Jiangyan-Zhao/BKP)
- **Paper:** [arXiv:2508.10447](https://arxiv.org/abs/2508.10447)
- **Presentation slides:** [`slides/BKP_Slides.pdf`](slides/BKP_Slides.pdf)

## Repository structure

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
│   ├── a_coverage.R
│   ├── data/
│   ├── figure/
│   └── result/
├── paper/
│   ├── TR_BKP.tex
│   ├── TR_BKP.pdf
│   ├── refs.bib
│   ├── jss.cls
│   └── jss.bst
├── slides/
│   ├── BKP_Slides.Rmd
│   ├── BKP_Slides.pdf
│   ├── preamble.tex
│   ├── ecnu.sty
│   ├── ecnu_logo.png
│   ├── ecnu_title.png
│   ├── Thanks.png
│   └── qr_website.png
├── BKP-paper.Rproj
├── renv.lock
└── README.md
```

## Requirements

The analyses were run in R using the package environment recorded in `renv.lock`. 
The current lockfile records `BKP` version 0.3.1.

To restore the package environment, open the repository as an R project and run:

``` r
renv::restore()
```

The main package used in the paper is `BKP`. 
The reproduction scripts also use several supporting packages, 
including `tgp`, `gplite`, `kernlab`, `pROC`, `mlbench`, 
`ggplot2`, `gridExtra`, `RiskMap`, `sf`, `terra`, `maps`, `rnaturalearth`, 
and related dependencies.

## Reproducing all analyses

To reproduce all illustrative examples, real-data applications, and the appendix simulation, run:

``` r
source("code/run_all.R")
```

The script `code/run_all.R` executes, in order:

1.  the nine illustrative examples in Section 4;
2.  the two real-data applications in Section 5;
3.  the predictive-coverage simulation in the appendix.

Generated figures are saved to:

``` text
code/figure/
```

Precomputed, intermediate, and summary numerical results are saved to:

``` text
code/result/
```

## Section 4: Illustrative examples

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

Each script can also be run separately. For example:

``` r
source("code/s4_ex1_bkp_1d_logistic.R")
```

## Section 5: Real-data applications

Section 5 includes:

``` text
s5_app1_loaloa_prevalence_mapping.R
s5_app2_mourning_warbler_sdm.R
```

The Loa loa prevalence mapping application can be reproduced using:

``` r
source("code/s5_app1_loaloa_prevalence_mapping.R")
```

This application uses the `loaloa` dataset from the `RiskMap` package and 
compares BKP with a logistic Gaussian process model.

The Mourning Warbler species distribution application can be reproduced using:

``` r
source("code/s5_app2_mourning_warbler_sdm.R")
```

The observation data and WorldClim raster files required for this application are stored under:

``` text
code/data/
```

The analysis compares BKP, TwinBKP, and a logistic Gaussian process model using eight bioclimatic covariates. 
It includes model fitting, test-set evaluation, ROC analysis, and raster-level geographic projection.

## Appendix: Predictive-coverage simulation

The appendix simulation can be reproduced using:

``` r
source("code/a_coverage.R")
```

The script compares:

- standard BKP without effective-sample-size calibration;
- BKP with Shepard effective-sample-size calibration;
- a logistic Gaussian process model.

The simulation considers sample sizes `n = 30` and `n = 100`. 
Pointwise interval coverage is evaluated at 2,000 fixed grid locations 
and at the original training locations over 100 independent simulation replications.

The appendix script generates:

``` text
code/figure/ex2_combined.pdf
code/result/coverage_summary.csv
```

## Timing experiment

The timing comparison in Example 3 can be computationally expensive. 
By default, the script reads the precomputed average timing results from:

``` text
code/result/elapsed_time_avg.csv
```

To rerun the full timing experiment, edit the following line in `code/s4_ex3_bkp_2d_goldstein_price.R`:

``` r
run_elapsed_time <- TRUE
```

The full timing experiment may take substantial time, 
especially for the optimized logistic Gaussian process benchmark.

## Computationally intensive components

A complete run of `code/run_all.R` may take substantial time.

The most computationally intensive components are:

- the full timing experiment in Example 3, when enabled;
- the Mourning Warbler raster-level prediction;
- the repeated appendix coverage simulation.

The precomputed timing results are included 
so that the default Example 3 workflow can be reproduced without rerunning the full timing experiment.

## Manuscript

The LaTeX source files are stored in:

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

through the graphic path setting in `TR_BKP.tex`.

## Presentation slides

The accompanying Beamer presentation is maintained in:

``` text
slides/
```

The main source file is:

``` text
slides/BKP_Slides.Rmd
```

A compiled version of the presentation is provided at:

``` text
slides/BKP_Slides.pdf
```

The presentation uses the custom Beamer configuration in `slides/preamble.tex` and `slides/ecnu.sty`. 
Slide-specific assets, including the ECNU logos, the closing image, 
and the BKP website QR code, are also stored in `slides/`.

The slides reuse the generated figures from:

``` text
code/figure/
```

and the bibliography from:

``` text
paper/refs.bib
```

To render the slides from the repository root, run:

``` r
rmarkdown::render("slides/BKP_Slides.Rmd")
```

Before rendering, restore the recorded R package environment:

``` r
renv::restore()
```

A working LaTeX distribution, such as TinyTeX or TeX Live, is also required.

## Working directory

Run the analysis scripts and render the presentation from the repository root, 
namely the directory containing `renv.lock`, `code/`, `paper/`, and `slides/`.

For example, run:

``` r
source("code/run_all.R")
```

or render the slides using:

``` r
rmarkdown::render("slides/BKP_Slides.Rmd")
```

Do not change the working directory to `code/` before running the scripts, 
because paths such as `code/data/`, `code/figure/`, and `code/result/` 
are defined relative to the repository root.

## Reproducibility notes

For a clean reproduction, avoid relying on saved R workspaces such as `.RData`. 
These files are ignored by `.gitignore`.

After restoring the package environment, the synchronization status can be checked using:

``` r
renv::status()
```

A synchronized project should report that no package versions differ 
between the project library and `renv.lock`.
