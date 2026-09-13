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

```text
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
├── renv/
│   └── activate.R
├── .Rprofile
├── BKP-paper.Rproj
├── renv.lock
└── README.md
```

All paths used by the reproduction scripts are defined relative to the
repository root.

## Computational environment

The analyses were run using the package environment recorded in
`renv.lock`. The current lockfile records R 4.6.1 and BKP 0.3.1.

The main package used in the paper is `BKP`. The reproduction scripts also
use several supporting packages, including `tgp`, `gplite`, `kernlab`,
`pROC`, `mlbench`, `ggplot2`, `gridExtra`, `RiskMap`, `sf`, `terra`,
`maps`, `rnaturalearth`, and related dependencies.

### Restoring the environment

When the repository is opened as an R project, `.Rprofile` activates the
local `renv` environment automatically. The recorded package environment can
then be restored using:

```r
renv::restore()
```

For a clean terminal session using `--vanilla`, `.Rprofile` is not loaded.
The `renv` activation script should therefore be sourced explicitly:

```bash
Rscript --vanilla -e "source('renv/activate.R'); renv::restore()"
```

After restoration, the synchronization status of the project library can be
checked using:

```r
renv::status()
```

## Reproducing all analyses

All reproduction commands should be run from the repository root.

To reproduce all illustrative examples, real-data applications, and the
appendix simulation from an interactive R session, run:

```r
source("code/run_all.R")
```

For a clean terminal session using `--vanilla`, explicitly activate the
project library before executing the complete analysis:

```bash
Rscript --vanilla -e "source('renv/activate.R'); source('code/run_all.R')"
```

The script `code/run_all.R` executes, in order:

1. the nine illustrative examples in Section 4;
2. the two real-data applications in Section 5;
3. the predictive-coverage simulation in the appendix.

Generated figures are saved to:

```text
code/figure/
```

Precomputed, intermediate, and summary numerical results are saved to:

```text
code/result/
```

## Section 4: Illustrative examples

Section 4 includes the following scripts:

```text
code/s4_ex1_bkp_1d_logistic.R
code/s4_ex2_bkp_1d_nonlinear.R
code/s4_ex3_bkp_2d_goldstein_price.R
code/s4_ex4_bkp_two_spirals_classification.R
code/s4_ex5_dkp_1d_multinomial.R
code/s4_ex6_dkp_2d_multinomial.R
code/s4_ex7_dkp_iris_classification.R
code/s4_ex8_twinbkp_1d_nonlinear.R
code/s4_ex9_twindkp_1d_multinomial.R
```

Each script can also be run separately from the repository root. For
example:

```r
source("code/s4_ex1_bkp_1d_logistic.R")
```

The examples cover one- and two-dimensional BKP and DKP models,
classification, posterior visualization, and the scalable TwinBKP and
TwinDKP approximations.

## Section 5: Real-data applications

Section 5 contains two real-data applications:

```text
code/s5_app1_loaloa_prevalence_mapping.R
code/s5_app2_mourning_warbler_sdm.R
```

### Loa loa prevalence mapping

The Loa loa prevalence-mapping application can be reproduced using:

```r
source("code/s5_app1_loaloa_prevalence_mapping.R")
```

This application uses the `loaloa` dataset from the `RiskMap` package and
compares BKP with a logistic Gaussian process model.

The analysis evaluates out-of-sample predictive performance and constructs
probability-scale prevalence and uncertainty maps over the study region.

### Mourning Warbler distribution modeling

The Mourning Warbler species-distribution application can be reproduced
using:

```r
source("code/s5_app2_mourning_warbler_sdm.R")
```

The observation data and WorldClim raster files required for this application
are stored under:

```text
code/data/
```

The analysis compares BKP, TwinBKP, and a logistic Gaussian process model
using eight bioclimatic covariates. Predictive performance is evaluated on
the predefined testing set using AUC and Brier score, followed by
raster-level geographic projection.

The numerical comparison reported in Table 3 is saved to:

```text
code/result/mourning_warbler_results.csv
```

The application also generates the corresponding geographic prediction and
uncertainty figures under:

```text
code/figure/
```

## Appendix: Predictive-coverage simulation

The appendix simulation can be reproduced using:

```r
source("code/a_coverage.R")
```

The simulation compares:

- standard BKP without effective-sample-size calibration;
- BKP with Shepard effective-sample-size calibration;
- a logistic Gaussian process model.

The simulation considers sample sizes `n = 30` and `n = 100`. Pointwise
interval coverage is evaluated at 2,000 fixed grid locations and at the
original training locations over 100 independent simulation replications.

The appendix script generates:

```text
code/figure/ex2_combined.pdf
code/result/coverage_summary.csv
```

## Timing experiment

The computational timing comparison in Example 3 can be expensive and is
therefore not rerun by default.

The control variable is defined in:

```text
code/s4_ex3_bkp_2d_goldstein_price.R
```

and defaults to:

```r
run_elapsed_time <- FALSE
```

With the default setting, the script reads the precomputed average timing
results from:

```text
code/result/elapsed_time_avg.csv
```

This allows the timing figure to be reproduced without rerunning the full
benchmark.

To rerun the complete timing experiment, change the setting in
`code/s4_ex3_bkp_2d_goldstein_price.R` to:

```r
run_elapsed_time <- TRUE
```

The full benchmark evaluates five methods, five sample sizes, and twenty
independent repetitions.

When enabled, the newly computed average timing results are saved to:

```text
code/result/elapsed_time_avg_new.csv
```

The `_new` suffix prevents the precomputed benchmark file

```text
code/result/elapsed_time_avg.csv
```

from being overwritten.

The corresponding timing figure is generated at:

```text
code/figure/elapsed_time.pdf
```

The full timing experiment may take substantial time, particularly for the
optimized logistic Gaussian process benchmark.

## Computationally intensive components

A complete execution of `code/run_all.R` may take substantial time.

The most computationally intensive components are:

- the full timing experiment in Example 3 when `run_elapsed_time <- TRUE`;
- the Mourning Warbler raster-level prediction;
- the repeated appendix coverage simulation.

The precomputed timing results are included so that the default Example 3
workflow and timing figure can be reproduced without rerunning the full
timing benchmark.

## Manuscript

The LaTeX source files are stored in:

```text
paper/
```

The main manuscript file is:

```text
paper/TR_BKP.tex
```

The manuscript is configured to read generated figures from:

```text
code/figure/
```

through the graphic-path setting in `TR_BKP.tex`.

The numerical results and manuscript figures should therefore be regenerated
from the analysis scripts before rebuilding the manuscript.

## Presentation slides

The accompanying Beamer presentation is maintained in:

```text
slides/
```

The main source file is:

```text
slides/BKP_Slides.Rmd
```

A compiled version of the presentation is provided at:

```text
slides/BKP_Slides.pdf
```

The presentation uses the custom Beamer configuration in
`slides/preamble.tex` and `slides/ecnu.sty`.

Slide-specific assets, including the ECNU logos, closing image, and BKP
website QR code, are also stored under `slides/`.

The slides reuse generated figures from:

```text
code/figure/
```

and the bibliography from:

```text
paper/refs.bib
```

To render the slides from the repository root, run:

```r
rmarkdown::render("slides/BKP_Slides.Rmd")
```

A working LaTeX distribution, such as TinyTeX or TeX Live, is also required.

## Working directory

All reproduction commands should be run from the repository root, namely the
directory containing:

```text
renv.lock
code/
paper/
slides/
```

For the complete analysis, run:

```r
source("code/run_all.R")
```

Do not change the working directory to `code/` before running the scripts,
because paths such as:

```text
code/data/
code/figure/
code/result/
```

are defined relative to the repository root.

For a clean terminal session, use:

```bash
Rscript --vanilla -e "source('renv/activate.R'); source('code/run_all.R')"
```

## Reproducibility notes

For a clean reproduction, avoid relying on saved R workspaces such as
`.RData`.

After restoring the package environment, synchronization can be checked
using:

```r
renv::status()
```

The replication archive and this repository use the same directory structure
and relative file paths.

All numerical results and manuscript figures are generated from the scripts
under `code/`, so the submitted code, saved numerical results, distributed
figures, and manuscript graphics can be regenerated from a single consistent
workflow.
