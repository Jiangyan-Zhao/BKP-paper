## -------------------------------------------------------------------------
## Section 4, Example 4: Two Spirals binary classification
##
## This script reproduces Example 4 in Section 4 of the manuscript.
## It considers a binary classification problem based on the Two Spirals
## benchmark dataset. The example compares BKP with a logistic Gaussian
## process (LGP) model fitted by the gplite package.
##
## The script contains four main parts:
##
##   1. Generate a noisy Two Spirals dataset and split it into training and
##      testing sets.
##
##   2. Fit a BKP classifier, compute test-set predictions, plot the ROC curve,
##      and visualize the BKP predictive mean and variance surfaces.
##
##   3. Fit an LGP classifier using expectation propagation.
##
##   4. Compute LGP test-set predictions, plot the ROC curve, and visualize
##      the LGP predictive mean and variance surfaces on a regular grid.
##
## Output figures are saved to code/figure/.
## -------------------------------------------------------------------------


## Load required packages.
## BKP provides the Beta Kernel Process classifier and plotting methods.
## gplite is used to fit the logistic Gaussian process comparison model.
## mlbench provides the Two Spirals benchmark dataset generator.
## pROC is used for ROC and AUC calculations.
## gridExtra is used to arrange the LGP mean and variance plots side by side.
library(BKP)
library(tgp)
library(gplite)
library(mlbench)
library(pROC)
library(gridExtra)

## Set console formatting options so that printed output follows the style
## used in the manuscript's R code examples.
options(
  prompt = "R> ",
  continue = "+  ",
  width = 70,
  useFancyQuotes = FALSE
)

## Create the output directory for figures if it does not already exist.
## This makes the script runnable from a clean clone of the repository.
dir.create("code/figure", recursive = TRUE, showWarnings = FALSE)


## -------------------------------------------------------------------------
## Part I: Data generation and train-test split
## -------------------------------------------------------------------------

## Set the random seed to make the simulated Two Spirals dataset and
## train-test split reproducible.
set.seed(123)

## Specify the total number of observations.
n <- 120

## Generate a noisy Two Spirals dataset with two complete rotations.
## The additive Gaussian noise has standard deviation 0.05.
data <- mlbench.spirals(n, cycles = 2, sd = 0.05)

## Randomly select 70% of the observations as the training set.
train_indices <- sample(1:n, 0.7 * n)

## Extract training inputs and convert the two class labels to 0/1 coding.
## BKP expects binary responses to be encoded as 0 and 1.
X_train <- data$x[train_indices, ]
y_train <- as.numeric(data$classes[train_indices]) - 1

## Extract testing inputs and convert the test labels to 0/1 coding.
X_test <- data$x[-train_indices, ]
y_test <- as.numeric(data$classes[-train_indices]) - 1

## For binary classification, each observation corresponds to one Bernoulli
## trial, so the binomial trial size is one for every training point.
m <- rep(1, length(train_indices))

## Define the two-dimensional input domain used for normalization and plotting.
Xbounds <- rbind(c(-1.7, 1.7), c(-1.7, 1.7))


## -------------------------------------------------------------------------
## Part II: BKP classifier
## -------------------------------------------------------------------------

## Fit the BKP model for binary classification.
## A fixed beta prior with r0 = 0.1 and p0 = 0.5 is used. The small value of
## r0 gives the prior weak influence while keeping the posterior update
## numerically stable. The log-loss criterion is used for kernel
## hyperparameter tuning because this is a classification task.
BKP_model_Class <- fit_BKP(
  X_train, y_train, m,
  Xbounds = Xbounds,
  prior = "fixed",
  r0 = 0.1,
  p0 = 0.5,
  loss = "log_loss"
)

## Compute BKP posterior predictions on the test inputs.
## The returned posterior mean probabilities are used for ROC analysis.
prediction <- predict(BKP_model_Class, X_test)

## Plot the ROC curve for BKP and report the corresponding AUC in the title.
pdf("code/figure/ex4roc.pdf", width = 6, height = 6)
roc_curve_BKP <- roc(y_test, prediction$mean)
plot(
  roc_curve_BKP,
  main = paste(
    "ROC curve for BKP (AUC = ",
    round(auc(roc_curve_BKP), 3),
    ")",
    sep = ""
  ),
  col = "#1c61b6"
)
dev.off()

## Visualize the BKP predictive mean and variance surfaces.
## For two-dimensional binary classification, the plot() method displays
## posterior summaries of the class probability over a regular grid, with
## the training observations overlaid.
pdf("code/figure/ex4.pdf", width = 13, height = 6)
plot(BKP_model_Class)
dev.off()


## -------------------------------------------------------------------------
## Part III: Logistic Gaussian process classifier
## -------------------------------------------------------------------------

## Fit an LGP model for comparison.
## The covariance function is the squared exponential covariance cf_sexp(),
## the likelihood is Bernoulli, and posterior approximation is performed
## using expectation propagation through approx_ep().
gp <- gp_init(cf = cf_sexp(), lik = lik_bernoulli())
gp <- gp_optim(gp, X_train, y_train, verbose = FALSE)

## Compute LGP posterior predictions on the test inputs.
prediction_gp <- gp_pred(gp, as.matrix(X_test), transform = TRUE)

## Plot the ROC curve for LGP and report the corresponding AUC in the title.
pdf("code/figure/ex4rocLGP.pdf", width = 6, height = 6)
roc_curve_LGP <- roc(y_test, prediction_gp$mean)
plot(
  roc_curve_LGP,
  main = paste(
    "ROC curve for LGP (AUC = ",
    round(auc(roc_curve_LGP), 3),
    ")",
    sep = ""
  ),
  col = "#1c61b6"
)
dev.off()


## -------------------------------------------------------------------------
## Part IV: LGP predictive surface visualization
## -------------------------------------------------------------------------

## Construct a regular two-dimensional grid over the plotting domain.
## This grid is used to visualize the LGP predictive mean and variance
## surfaces in the same style as the BKP plots.
x1_seq <- seq(-1.7, 1.7, length.out = 80)
x2_seq <- seq(-1.7, 1.7, length.out = 80)
grid <- expand.grid(x1 = x1_seq, x2 = x2_seq)

## Predict LGP posterior mean and variance on the visualization grid.
prediction <- gp_pred(gp, as.matrix(grid), transform = TRUE, var = TRUE)

## Store the grid coordinates and LGP posterior summaries in a data frame
## compatible with the BKP two-dimensional plotting helper.
df <- data.frame(
  x1 = grid$x1,
  x2 = grid$x2,
  Mean = prediction$mean,
  Variance = prediction$var
)

## Use BKP's internal two-dimensional plotting helper to produce contour-style
## plots for the LGP predictive mean and variance. The function is non-exported
## but documented in the package; it is used here only to match the manuscript
## figure style.
p1 <- BKP:::my_2D_plot_fun("Mean", "Posterior Mean", df, X = X_train, y = y_train)
p2 <- BKP:::my_2D_plot_fun("Variance", "Posterior Variance", df, X = X_train, y = y_train)

## Arrange the LGP predictive mean and variance plots side by side.
pdf("code/figure/ex4LGP.pdf", width = 13, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()