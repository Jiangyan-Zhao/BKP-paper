## -------------------------------------------------------------------------
## Section 5, Application 2: Mourning Warbler distribution modeling
##
## This script reproduces the second real-data application in Section 5 of
## the manuscript. It analyzes presence-absence observations for the Mourning
## Warbler using eight bioclimatic predictors.
##
## Geographic coordinates are used only for visualization. Model fitting is
## based exclusively on environmental covariates, so the fitted response
## surfaces describe the association between species occurrence and climate
## rather than direct spatial smoothing over longitude and latitude.
##
## The observation data contain a predefined spatially separated training and
## testing split. Three models are compared:
##
##   1. A full Beta Kernel Process (BKP) classifier.
##
##   2. A scalable TwinBKP classifier based on a twinning-selected global
##      subset and prediction-specific local updates.
##
##   3. A logistic Gaussian process (LGP) classifier fitted using gplite.
##
## Predictive performance is evaluated on the withheld testing set using:
##
##   - area under the ROC curve (AUC);
##   - Brier score;
##   - elapsed fitting-and-prediction time.
##
## The fitted environmental relationships are subsequently projected onto a
## North American bioclimatic raster to produce geographic maps of predicted
## presence probability and predictive uncertainty.
##
## All paths are relative to the repository root. Generated figures are saved
## to code/figure/.
## -------------------------------------------------------------------------


## -------------------------------------------------------------------------
## Required packages
## -------------------------------------------------------------------------

## BKP provides the full BKP and scalable TwinBKP fitting and prediction
## methods.
##
## gplite provides the logistic Gaussian process comparison model.
##
## terra is used to read, crop, manipulate, and project the WorldClim raster
## layers.
##
## pROC provides ROC curves and AUC calculations.
##
## ggplot2, gridExtra, maps, dplyr, and RColorBrewer support map construction
## and figure assembly.
library(BKP)
library(gplite)
library(ggplot2)
library(dplyr)
library(maps)
library(pROC)
library(gridExtra)
library(terra)
library(RColorBrewer)


## Set console formatting options so that printed output follows the style
## used in the manuscript's R code examples.
options(
  prompt = "R> ",
  continue = "+  ",
  width = 70,
  useFancyQuotes = FALSE
)

## Create the output directory if it does not already exist. This allows the
## script to be run independently from a clean clone of the repository.
dir.create(
  "code/figure",
  recursive = TRUE,
  showWarnings = FALSE
)

## Set the random seed for reproducible stochastic operations during model
## fitting and hyperparameter optimization.
set.seed(123)


## -------------------------------------------------------------------------
## Observation data and predefined spatial split
## -------------------------------------------------------------------------

## Read the prepared Mourning Warbler dataset.
##
## The file contains:
##
##   - route identifiers;
##   - longitude and latitude;
##   - predefined data-split labels;
##   - eight WorldClim bioclimatic variables;
##   - the binary presence-absence response y.
data <- read.csv(
  "code/data/Mourning_Warbler.csv",
  stringsAsFactors = FALSE
)

## Print basic information about the dataset and its predefined split.
cat(
  sprintf(
    "Data dimensions: %d rows x %d columns\n",
    nrow(data),
    ncol(data)
  )
)

cat("Split distribution:\n")
print(table(data$split))


## Identify observations assigned to the predefined training and testing sets.
##
## The split is spatially structured rather than randomly generated. This
## produces a more demanding evaluation because the test observations occupy
## geographically separated regions.
train_idx <- which(
  data$split == "training"
)

test_idx <- which(
  data$split == "testing"
)

cat(
  sprintf(
    "Training set: %d | Test set: %d\n",
    length(train_idx),
    length(test_idx)
  )
)


## Identify the eight WorldClim covariates from their common column-name
## prefix. Longitude and latitude are deliberately excluded from the model
## input matrix and retained only for geographic visualization.
bio_cols <- grep(
  "^wc2",
  names(data),
  value = TRUE
)

cat(
  sprintf(
    "Covariates: %s\n",
    paste(
      bio_cols,
      collapse = ", "
    )
  )
)


## Construct the environmental input matrix.
X_all <- as.matrix(
  data[, bio_cols, drop = FALSE]
)

## Construct training inputs and binary responses.
X_train <- X_all[
  train_idx,
  ,
  drop = FALSE
]

y_train <- data$y[
  train_idx
]

## Each observation is a single Bernoulli presence-absence outcome, so the
## binomial trial size is one at every training location.
m_train <- rep(
  1L,
  length(y_train)
)

## Construct testing inputs and responses.
X_test <- X_all[
  test_idx,
  ,
  drop = FALSE
]

y_test <- data$y[
  test_idx
]

## Define the normalization bounds using only the training covariates. Each row
## of Xbounds contains the minimum and maximum of one bioclimatic variable.
Xbounds <- t(
  apply(
    X_train,
    2,
    range
  )
)


## -------------------------------------------------------------------------
## Bioclimatic raster preparation
## -------------------------------------------------------------------------

cat("\n===== Loading climate rasters =====\n")

## Construct the path to the directory containing the WorldClim raster layers.
tif_dir <- "code/data/climate/wc2.1_10m"

## Match each model covariate to its corresponding raster file. The raster
## layers must appear in the same order as the columns of X_train.
tif_files <- file.path(
  tif_dir,
  paste0(
    bio_cols,
    ".tif"
  )
)

## Read the eight bioclimatic raster layers as a multi-layer SpatRaster.
na_clim <- rast(
  tif_files
)

## Crop the global WorldClim layers to the North American study region.
study_extent <- ext(
  -170.15,
  -46.70,
  22.23,
  71.90
)

na_clim <- crop(
  na_clim,
  study_extent
)

## Assign the model covariate names to the raster layers to ensure that raster
## values are supplied to the fitted models in the correct variable order.
names(na_clim) <- bio_cols

cat(
  sprintf(
    "Cropped raster: %d rows x %d columns x %d layers\n",
    nrow(na_clim),
    ncol(na_clim),
    nlyr(na_clim)
  )
)


## -------------------------------------------------------------------------
## Plotting helpers and color definitions
## -------------------------------------------------------------------------

## Define a helper that copies the geometry of the first climate layer and
## inserts a vector of model predictions. This converts model output back into
## a SpatRaster suitable for geographic plotting.
make_rast <- function(vals, template) {
  r <- rast(
    template[[1]]
  )
  
  values(r) <- vals
  
  return(r)
}


## Construct the palettes used for presence probabilities, predictive
## variances, and observed training presences and absences.
purples <- colorRampPalette(
  brewer.pal(
    9,
    "Purples"
  )[-(1:2)]
)

greens <- colorRampPalette(
  brewer.pal(
    9,
    "Greens"
  )[-(1:2)]
)

prescol <- purples(3)[3]
abscol <- purples(3)[2]


## Separate the geographic coordinates into training presences, training
## absences, and withheld testing locations.
pts_train_pres <- data[
  data$y == 1 & data$split == "training",
  c("lon", "lat"),
  drop = FALSE
]

pts_train_abs <- data[
  data$y == 0 & data$split == "training",
  c("lon", "lat"),
  drop = FALSE
]

pts_test <- data[
  data$split == "testing",
  c("lon", "lat"),
  drop = FALSE
]


## -------------------------------------------------------------------------
## Observed distribution and evaluation split
## -------------------------------------------------------------------------

## Construct a land mask from one of the climate layers. All non-missing cells
## are assigned a constant value so that the map displays the geographic
## support of the environmental raster without displaying the climate variable
## itself.
bg_rast <- na_clim[["wc2.1_10m_bio_5"]]

bg_vals <- values(
  bg_rast,
  na.rm = FALSE
)

bg_vals[
  !is.na(bg_vals)
] <- 1

values(bg_rast) <- bg_vals

bg_df <- as.data.frame(
  bg_rast,
  xy = TRUE
)

colnames(bg_df)[3] <- "land"


## Construct the distribution map used in the manuscript.
##
## Dark purple points indicate training presences.
## Light purple points indicate training absences.
## White points indicate withheld testing observations.
map_distribution <- ggplot() +
  geom_tile(
    data = bg_df,
    aes(
      x = x,
      y = y,
      fill = land
    )
  ) +
  scale_fill_gradientn(
    colours = grey(0.8),
    na.value = "transparent",
    guide = "none"
  ) +
  geom_point(
    data = pts_test,
    aes(
      x = lon,
      y = lat
    ),
    shape = 21,
    size = 0.6,
    stroke = 0.15,
    fill = "white",
    colour = grey(0.4)
  ) +
  geom_point(
    data = pts_train_abs,
    aes(
      x = lon,
      y = lat
    ),
    shape = 21,
    size = 0.6,
    stroke = 0.15,
    fill = abscol,
    colour = abscol
  ) +
  geom_point(
    data = pts_train_pres,
    aes(
      x = lon,
      y = lat
    ),
    shape = 21,
    size = 0.6,
    stroke = 0.15,
    fill = prescol,
    colour = prescol
  ) +
  coord_quickmap(
    xlim = c(-170, -50),
    ylim = c(15, 75),
    expand = FALSE
  ) +
  labs(
    title = "Distribution data - Mourning Warbler",
    x = NULL,
    y = NULL
  ) +
  theme_void(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

## Save the distribution and evaluation-split figure.
ggsave(
  "code/figure/mourning_warbler_map_distribution.pdf",
  plot = map_distribution,
  width = 6.5,
  height = 3.8,
  dpi = 300
)


## -------------------------------------------------------------------------
## Model fitting and test-set prediction
## -------------------------------------------------------------------------

## Initialize containers for ROC objects, predicted probabilities, and elapsed
## times so that all three methods can be evaluated using a common workflow.
all_rocs <- list()
all_preds <- list()
times <- numeric(0)


## -------------------------------------------------------------------------
## Full BKP classifier
## -------------------------------------------------------------------------

cat("\n===== BKP =====\n")

## Fit the full BKP classifier using:
##
##   - a fixed symmetric beta prior with r0 = 0.1;
##   - the log-loss criterion for kernel hyperparameter tuning;
##   - an isotropic Gaussian kernel.
##
## The reported elapsed time includes both model fitting and prediction on the
## withheld test set.
t1 <- system.time({
  bkp_fit <- fit_BKP(
    X_train,
    y_train,
    m_train,
    Xbounds = Xbounds,
    prior = "fixed",
    r0 = 0.1,
    loss = "log_loss",
    kernel = "gaussian",
    isotropic = TRUE
  )
  
  bkp_pred_test <- predict(
    bkp_fit,
    Xnew = X_test
  )
})

times <- c(
  times,
  BKP = unname(
    t1["elapsed"]
  )
)

## Construct the BKP ROC curve and store its predicted probabilities.
bkp_roc <- roc(
  y_test,
  as.vector(
    bkp_pred_test$mean
  )
)

all_rocs$BKP <- bkp_roc
all_preds$BKP <- as.vector(
  bkp_pred_test$mean
)


## -------------------------------------------------------------------------
## Scalable TwinBKP classifier
## -------------------------------------------------------------------------

cat("\n===== TwinBKP =====\n")

## Fit the scalable TwinBKP classifier using:
##
##   - the same fixed prior and log-loss criterion as BKP;
##   - an isotropic Gaussian global kernel;
##   - a compactly supported Wendland local kernel;
##   - the package defaults for global-subset and local-neighbour sizes.
t2 <- system.time({
  twin_fit <- fit_TwinBKP(
    X_train,
    y_train,
    m_train,
    Xbounds = Xbounds,
    prior = "fixed",
    r0 = 0.1,
    loss = "log_loss",
    global_kernel = "gaussian",
    local_kernel = "wendland",
    isotropic = TRUE
  )
  
  twin_pred_test <- predict(
    twin_fit,
    Xnew = X_test
  )
})

times <- c(
  times,
  TwinBKP = unname(
    t2["elapsed"]
  )
)

## Construct the TwinBKP ROC curve and store its predicted probabilities.
twin_roc <- roc(
  y_test,
  as.vector(
    twin_pred_test$mean
  )
)

all_rocs$TwinBKP <- twin_roc
all_preds$TwinBKP <- as.vector(
  twin_pred_test$mean
)


## -------------------------------------------------------------------------
## Logistic Gaussian process classifier
## -------------------------------------------------------------------------

cat("\n===== LGP =====\n")

## Fit the logistic Gaussian process benchmark using:
##
##   - a squared-exponential covariance function;
##   - a Bernoulli likelihood;
##   - gplite's default optimization procedure.
##
## The elapsed time includes both fitting and prediction on the test set.
t3 <- system.time({
  gp <- gp_init(
    cf = cf_sexp(),
    lik = lik_bernoulli()
  )
  
  gp <- gp_optim(
    gp,
    X_train,
    y_train,
    verbose = FALSE
  )
  
  gp_pred_test <- gp_pred(
    gp,
    X_test,
    transform = TRUE
  )
})

times <- c(
  times,
  LGP = unname(
    t3["elapsed"]
  )
)

## Construct the LGP ROC curve and store its predicted probabilities.
gp_roc <- roc(
  y_test,
  as.vector(
    gp_pred_test$mean
  )
)

all_rocs$LGP <- gp_roc
all_preds$LGP <- as.vector(
  gp_pred_test$mean
)


## -------------------------------------------------------------------------
## Predictive performance comparison
## -------------------------------------------------------------------------

## Report AUC, Brier score, and elapsed time for each fitted model.
##
## AUC summarizes discrimination between presences and absences.
##
## The Brier score evaluates the squared error of the predicted probabilities
## and therefore reflects both probabilistic calibration and accuracy.
cat("\n=============================================================\n")
cat("  Comparison: Mourning Warbler (8 bio vars, r0 = 0.1)\n")
cat("=============================================================\n")

for (nm in names(all_rocs)) {
  auc_value <- as.numeric(
    auc(
      all_rocs[[nm]]
    )
  )
  
  brier_score <- mean(
    (all_preds[[nm]] - y_test)^2
  )
  
  cat(
    sprintf(
      "  %-12s | AUC: %.3f | Brier: %.6f | Time: %.2f s\n",
      nm,
      auc_value,
      brier_score,
      times[nm]
    )
  )
}


## Save elapsed times to CSV.
dir.create(
  "code/result",
  recursive = TRUE,
  showWarnings = FALSE
)

times_df <- data.frame(
  method = names(times),
  elapsed_seconds = unname(times),
  stringsAsFactors = FALSE
)

write.csv(
  times_df,
  "code/result/mourning_warbler_times.csv",
  row.names = FALSE
)


## -------------------------------------------------------------------------
## ROC comparison figure
## -------------------------------------------------------------------------

## Construct a common ROC figure for BKP, TwinBKP, and LGP. The legend reports
## the test-set AUC for each method.
n_methods <- length(
  all_rocs
)

cols <- c(
  "#1c61b6",
  "#fdae61",
  "#d7191c"
)

nm_pad <- format(
  names(all_rocs),
  justify = "none",
  width = max(
    nchar(
      names(all_rocs)
    )
  )
)

legend_lines <- sprintf(
  "%s  (AUC=%.3f)",
  nm_pad,
  vapply(
    all_rocs,
    function(x) {
      as.numeric(
        auc(x)
      )
    },
    numeric(1)
  )
)

pdf(
  "code/figure/mourning_warbler_roc_comparison.pdf",
  width = 6,
  height = 6
)

plot(
  all_rocs[[1]],
  col = cols[1],
  lwd = 2,
  main = "ROC Curve: Mourning Warbler"
)

for (i in seq_len(n_methods)[-1]) {
  lines(
    all_rocs[[i]],
    col = cols[i],
    lwd = 2
  )
}

legend(
  "bottomright",
  legend = legend_lines,
  col = cols[seq_len(n_methods)],
  lwd = 2,
  bty = "o",
  cex = 0.85,
  bg = "white"
)

dev.off()


## -------------------------------------------------------------------------
## Projection onto the climate raster
## -------------------------------------------------------------------------

cat("\n===== Grid prediction =====\n")

## Extract the eight bioclimatic values from every raster cell.
grid_vals <- values(
  na_clim,
  na.rm = FALSE
)

## Reinforce the required covariate names and ordering.
colnames(grid_vals) <- bio_cols

grid_n <- nrow(
  grid_vals
)

## A raster cell is eligible for prediction only when all eight bioclimatic
## variables are observed. Missing cells, primarily ocean cells, remain NA in
## the final maps.
grid_valid <- complete.cases(
  grid_vals
)

cat(
  sprintf(
    "Grid points: %d (valid land cells %d / missing cells %d)\n",
    grid_n,
    sum(grid_valid),
    sum(!grid_valid)
  )
)


## Initialize a list for storing raster-level predictive means and variances.
raster_preds <- list()

for (method in names(all_preds)) {
  cat(
    sprintf(
      "  %s grid prediction...\n",
      method
    )
  )
  
  ## Allocate full-length vectors so that predictions can later be inserted
  ## directly into the original raster geometry.
  gm <- rep(
    NA_real_,
    grid_n
  )
  
  gv <- rep(
    NA_real_,
    grid_n
  )
  
  if (any(grid_valid)) {
    ## Construct the complete-case covariate matrix in the same variable order
    ## as the training data.
    X_grid <- grid_vals[
      grid_valid,
      bio_cols,
      drop = FALSE
    ]
    
    if (method == "LGP") {
      ## Evaluate LGP raster predictions in chunks to control memory use.
      ##
      ## valid_idx maps the rows of X_grid back to their corresponding
      ## positions in the full raster vectors gm and gv.
      chunk_size <- 5000
      n_valid <- nrow(X_grid)
      valid_idx <- which(grid_valid)
      
      for (k in seq(
        from = 1,
        to = n_valid,
        by = chunk_size
      )) {
        idx_end <- min(
          k + chunk_size - 1,
          n_valid
        )
        
        idx <- k:idx_end
        
        chunk <- X_grid[
          idx,
          ,
          drop = FALSE
        ]
        
        pred_chunk <- gp_pred(
          gp,
          chunk,
          transform = TRUE,
          var = TRUE
        )
        
        gm[
          valid_idx[idx]
        ] <- as.vector(
          pred_chunk$mean
        )
        
        gv[
          valid_idx[idx]
        ] <- as.vector(
          pred_chunk$var
        )
      }
    } else {
      ## BKP and TwinBKP use the common predict() interface and process the
      ## complete matrix of valid environmental covariates directly.
      fit_obj <- if (method == "BKP") {
        bkp_fit
      } else {
        twin_fit
      }
      
      pred_grid <- predict(
        fit_obj,
        Xnew = X_grid
      )
      
      gm[
        grid_valid
      ] <- as.vector(
        pred_grid$mean
      )
      
      gv[
        grid_valid
      ] <- as.vector(
        pred_grid$variance
      )
    }
  }
  
  raster_preds[[method]] <- list(
    mean = gm,
    variance = gv
  )
}


## -------------------------------------------------------------------------
## Geographic probability and uncertainty maps
## -------------------------------------------------------------------------

cat("\n===== Generating maps =====\n")

## Define panel labels following the order used in the combined manuscript
## figure:
##
##   (a) BKP predicted distribution;
##   (b) BKP prediction uncertainty;
##   (c) TwinBKP predicted distribution;
##   (d) TwinBKP prediction uncertainty;
##   (e) LGP predicted distribution;
##   (f) LGP prediction uncertainty.
sub_labels <- c(
  "(a)",
  "(b)",
  "(c)",
  "(d)",
  "(e)",
  "(f)"
)

sub_idx <- 1

## Initialize lists for storing model-specific probability and variance plots.
all_probs <- list()
all_vars <- list()


for (method in names(raster_preds)) {
  cat(
    sprintf(
      "  Generating %s maps...\n",
      method
    )
  )
  
  rp <- raster_preds[[method]]
  
  
  ## -----------------------------------------------------------------------
  ## Predicted probability-of-presence map
  ## -----------------------------------------------------------------------
  
  ## Convert the vector of predicted probabilities back to a raster.
  r_prob <- make_rast(
    rp$mean,
    na_clim
  )
  
  ## Convert the raster to a data frame for ggplot2.
  df_prob <- as.data.frame(
    r_prob,
    xy = TRUE
  )
  
  colnames(df_prob)[3] <- "prob"
  
  ## Construct the predicted probability-of-presence map. All models use the
  ## common probability range [0, 1].
  p_prob <- ggplot() +
    geom_tile(
      data = df_prob,
      aes(
        x = x,
        y = y,
        fill = prob
      )
    ) +
    scale_fill_gradientn(
      colours = purples(1000),
      name = "Probability of presence",
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        0.2
      ),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = grid::unit(
          0.2,
          "cm"
        ),
        barheight = grid::unit(
          4,
          "cm"
        ),
        frame.colour = NA,
        ticks.colour = NA,
        ticks.linewidth = 0,
        title.position = "right"
      )
    ) +
    coord_quickmap(
      xlim = c(-170, -50),
      ylim = c(15, 75),
      expand = FALSE
    ) +
    labs(
      title = paste(
        sub_labels[sub_idx],
        method,
        "- Predicted distribution"
      ),
      x = NULL,
      y = NULL
    ) +
    theme_void(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.key.width = grid::unit(
        0.2,
        "cm"
      ),
      legend.key.height = grid::unit(
        4,
        "cm"
      ),
      legend.ticks.length = grid::unit(
        0.2,
        "cm"
      ),
      legend.title = element_text(
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        size = 10
      )
    )
  
  sub_idx <- sub_idx + 1
  
  ## Individual probability maps can be saved by uncommenting the following
  ## code.
  # ggsave(
  #   sprintf(
  #     "code/figure/mourning_warbler_map_%s_prob.pdf",
  #     tolower(
  #       gsub(
  #         " ",
  #         "_",
  #         method
  #       )
  #     )
  #   ),
  #   p_prob,
  #   width = 6.5,
  #   height = 3.8,
  #   dpi = 300
  # )
  
  
  ## -----------------------------------------------------------------------
  ## Predictive variance map
  ## -----------------------------------------------------------------------
  
  ## Convert the vector of predictive variances back to a raster.
  r_var <- make_rast(
    rp$variance,
    na_clim
  )
  
  ## Convert the variance raster to a data frame for ggplot2.
  df_var <- as.data.frame(
    r_var,
    xy = TRUE
  )
  
  colnames(df_var)[3] <- "var"
  
  ## Determine a method-specific upper limit for the variance map.
  vmax <- max(
    rp$variance,
    na.rm = TRUE
  )
  
  vlim <- ceiling(
    vmax * 10
  ) / 10
  
  ## Construct the predictive variance map.
  p_var <- ggplot() +
    geom_tile(
      data = df_var,
      aes(
        x = x,
        y = y,
        fill = var
      )
    ) +
    scale_fill_gradientn(
      colours = greens(1000),
      name = paste(
        "Variance of",
        method,
        "posterior"
      ),
      limits = c(0, vlim),
      breaks = seq(
        0,
        vlim,
        0.1
      ),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = grid::unit(
          0.2,
          "cm"
        ),
        barheight = grid::unit(
          4,
          "cm"
        ),
        frame.colour = NA,
        ticks.colour = NA,
        ticks.linewidth = 0,
        title.position = "right"
      )
    ) +
    coord_quickmap(
      xlim = c(-170, -50),
      ylim = c(15, 75),
      expand = FALSE
    ) +
    labs(
      title = paste(
        sub_labels[sub_idx],
        method,
        "- Prediction uncertainty"
      ),
      x = NULL,
      y = NULL
    ) +
    theme_void(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.key.width = grid::unit(
        0.2,
        "cm"
      ),
      legend.key.height = grid::unit(
        4,
        "cm"
      ),
      legend.ticks.length = grid::unit(
        0.2,
        "cm"
      ),
      legend.title = element_text(
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        size = 10
      )
    )
  
  sub_idx <- sub_idx + 1
  
  ## Individual variance maps can be saved by uncommenting the following code.
  # ggsave(
  #   sprintf(
  #     "code/figure/mourning_warbler_map_%s_var.pdf",
  #     tolower(
  #       gsub(
  #         " ",
  #         "_",
  #         method
  #       )
  #     )
  #   ),
  #   p_var,
  #   width = 6.5,
  #   height = 3.8,
  #   dpi = 300
  # )
  
  ## Store the model-specific plots for final figure assembly.
  all_probs[[method]] <- p_prob
  all_vars[[method]] <- p_var
}


## -------------------------------------------------------------------------
## Combined manuscript figure
## -------------------------------------------------------------------------

## Arrange the BKP, TwinBKP, and LGP probability and variance maps into the
## 3-by-2 panel layout used in the manuscript.
cat("\n===== Generating combined map (3x2) =====\n")

p_combined <- grid.arrange(
  all_probs$BKP,
  all_vars$BKP,
  all_probs$TwinBKP,
  all_vars$TwinBKP,
  all_probs$LGP,
  all_vars$LGP,
  ncol = 2,
  nrow = 3
)

## Save the combined spatial prediction figure.
ggsave(
  "code/figure/mourning_warbler_map_combined.pdf",
  plot = p_combined,
  width = 13.0,
  height = 11.4,
  dpi = 300
)

cat("\n===== Analysis complete =====\n")