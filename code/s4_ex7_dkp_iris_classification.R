## -------------------------------------------------------------------------
## Section 4, Example 7: Iris multi-class classification
##
## This script reproduces Example 7 in Section 4 of the manuscript.
## It illustrates the Dirichlet Kernel Process (DKP) model for a three-class
## classification problem based on the Iris dataset.
##
## For visualization, only the first two covariates are used:
## sepal length and sepal width. The example compares DKP with a Gaussian
## process classifier fitted by kernlab::gausspr().
##
## The script contains four main parts:
##
##   1. Prepare the Iris data, construct one-hot class indicators, and split
##      the data into training and testing sets.
##
##   2. Fit a DKP classifier, visualize the fitted decision regions, and
##      compute one-vs-rest ROC curves.
##
##   3. Fit a Gaussian process classifier using gausspr() for comparison.
##
##   4. Visualize the Gaussian process classifier's predicted classes and
##      maximum predicted probabilities on a regular grid.
##
## Output figures are saved to code/figure/.
## -------------------------------------------------------------------------


## Load required packages.
## BKP provides the Dirichlet Kernel Process fitting, prediction, and plotting
## methods. pROC is used for one-vs-rest ROC and multiclass AUC calculations.
## kernlab provides gausspr(), which is used as the Gaussian process classifier
## comparison model. gridExtra is used to arrange multiple ggplot objects.
library(BKP)
library(pROC)
library(kernlab)
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
## Part I: Data preparation
## -------------------------------------------------------------------------

## Set the random seed to make the train-test split reproducible.
set.seed(123)

## Load the built-in Iris dataset.
data(iris)

## Use the first two covariates for two-dimensional visualization:
## sepal length and sepal width.
X <- as.matrix(iris[, 1:2])

## Define the input bounds used for normalization and plotting.
## The bounds slightly extend the observed range of sepal length and sepal
## width, producing a stable rectangular plotting domain.
Xbounds <- rbind(c(4.2, 8), c(1.9, 4.5))

## Extract the three species labels.
labels <- iris$Species

## Convert the factor labels into a one-hot response matrix.
## Each row of Y contains a single 1 indicating the observed species and
## zeros for the other two species. This is the response format expected by
## fit_DKP() for multi-class classification.
Y <- model.matrix(~ labels - 1)

## Randomly select 70% of the observations as the training set.
train_indices <- sample(1:nrow(iris), 0.7 * nrow(iris))

## Construct training inputs, one-hot responses, and labels.
X_train <- X[train_indices, ]
Y_train <- Y[train_indices, ]
labels_train <- labels[train_indices]

## Construct testing inputs, one-hot responses, and labels.
X_test <- X[-train_indices, ]
Y_test <- Y[-train_indices, ]
labels_test <- labels[-train_indices]


## -------------------------------------------------------------------------
## Part II: DKP classifier
## -------------------------------------------------------------------------

## Fit the DKP model for three-class classification.
## The log-loss criterion is used because this is a classification task.
## A fixed Dirichlet prior with equal class probabilities is specified through
## p0 = rep(1/3, 3). The prior precision r0 controls the strength of this
## baseline class-probability vector.
DKP_model_Class <- fit_DKP(
  X_train, Y_train,
  Xbounds = Xbounds,
  loss = "log_loss",
  prior = "fixed",
  r0 = 0.01,
  p0 = rep(1/3, 3)
)

## Visualize the DKP classification result.
## The plot() method displays the predicted class regions and the maximum
## predicted class probability over the two-dimensional input space.
pdf("code/figure/ex7.pdf", width = 13, height = 6)
plot(DKP_model_Class)
dev.off()

## Compute DKP posterior class-probability predictions on the test inputs.
dkp_pred_probs <- predict(DKP_model_Class, X_test)$mean

## Assign class names to the columns of the predicted probability matrix.
## This is needed for the one-vs-rest ROC calculations below.
class_levels <- levels(labels_test)
colnames(dkp_pred_probs) <- class_levels

## Compute the multiclass ROC object and macro-average AUC using pROC.
multiclass_roc_dkp <- multiclass.roc(labels_test, dkp_pred_probs)

## Compute one-vs-rest ROC curves for each Iris species.
all_rocs <- list()
for (class_name in class_levels) {
  
  ## Convert the current class into a binary one-vs-rest label.
  labels_binary <- ifelse(labels_test == class_name, 1, 0)
  
  ## Extract the DKP predicted probability for the current class.
  probabilities <- dkp_pred_probs[, class_name]
  
  ## Compute the one-vs-rest ROC curve.
  roc_curve <- roc(labels_binary, probabilities)
  
  ## Store the ROC curve using the class name.
  all_rocs[[class_name]] <- roc_curve
}

## Plot the one-vs-rest ROC curves for DKP.
pdf("code/figure/ex7roc.pdf", width = 6, height = 6)
plot(
  all_rocs[[1]],
  col = "blue",
  lwd = 2,
  lty = 1,
  main = paste(
    "One-vs-Rest ROC curve for DKP (AUC =",
    round(auc(multiclass_roc_dkp), 3),
    ")",
    sep = ""
  )
)
lines(all_rocs[[2]], col = "red", lwd = 2, lty = 2)
lines(all_rocs[[3]], col = "black", lwd = 2, lty = 4)
legend(
  "bottomright",
  legend = class_levels,
  col = c("blue", "red", "black"),
  lwd = 2,
  lty = c(1, 2, 4),
  cex = 1.2
)
dev.off()


## -------------------------------------------------------------------------
## Part III: Gaussian process classifier via kernlab
## -------------------------------------------------------------------------

## Prepare a data frame for fitting gausspr().
## The Gaussian process classifier is fitted using the same two covariates
## and the original factor-valued species response.
iris_data <- data.frame(
  Sepal.Length = iris$Sepal.Length,
  Sepal.Width = iris$Sepal.Width,
  Species = iris$Species
)

## Construct training and testing data frames using the same split as the DKP
## model.
iris_train <- iris_data[train_indices, ]
iris_test <- iris_data[-train_indices, ]

## Fit the Gaussian process classifier.
## The radial basis function kernel is used, and kpar = "automatic" lets
## kernlab select the kernel parameter automatically.
gausspr_model <- gausspr(
  Species ~ .,
  data = iris_train,
  kernel = "rbfdot",
  kpar = "automatic"
)

## Compute predicted class probabilities on the test set.
lgp_pred_probs <- predict(
  gausspr_model,
  newdata = iris_test,
  type = "probabilities"
)

## Compute the multiclass ROC object and macro-average AUC for the Gaussian
## process classifier.
multiclass_roc_lgp <- multiclass.roc(iris_test$Species, lgp_pred_probs)

## Compute one-vs-rest ROC curves for each class.
all_rocs <- list()
for (class_name in class_levels) {
  
  ## Convert the current class into a binary one-vs-rest label.
  labels_binary <- ifelse(labels_test == class_name, 1, 0)
  
  ## Extract the Gaussian process predicted probability for the current class.
  probabilities <- lgp_pred_probs[, class_name]
  
  ## Compute the one-vs-rest ROC curve.
  roc_curve <- roc(labels_binary, probabilities)
  
  ## Store the ROC curve using the class name.
  all_rocs[[class_name]] <- roc_curve
}

## Plot the one-vs-rest ROC curves for the Gaussian process classifier.
pdf("code/figure/ex7rocLGP.pdf", width = 6, height = 6)
plot(
  all_rocs[[1]],
  col = "blue",
  lwd = 2,
  lty = 1,
  main = paste(
    "One-vs-Rest ROC curve for LGP (AUC =",
    round(auc(multiclass_roc_lgp), 3),
    ")",
    sep = ""
  )
)
lines(all_rocs[[2]], col = "red", lwd = 2, lty = 2)
lines(all_rocs[[3]], col = "black", lwd = 2, lty = 4)
legend(
  "bottomright",
  legend = class_levels,
  col = c("blue", "red", "black"),
  lty = c(1, 2, 4),
  lwd = 2,
  cex = 1.2
)
dev.off()


## -------------------------------------------------------------------------
## Part IV: Gaussian process predictive surface visualization
## -------------------------------------------------------------------------

## Construct a regular two-dimensional grid over the plotting domain.
## This grid is used to visualize predicted class regions and maximum
## predicted class probabilities for the Gaussian process classifier.
grid <- expand.grid(
  Sepal.Length = seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 80),
  Sepal.Width = seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 80)
)

## Predict class probabilities on the grid.
grid_predictions <- predict(gausspr_model, newdata = grid, type = "prob")

## Convert the predicted probability matrix into hard predicted classes by
## selecting the class with the largest predicted probability at each grid
## location.
class <- max.col(grid_predictions)

## Store the grid coordinates, predicted class labels, and maximum predicted
## probabilities in a data frame compatible with the BKP plotting helper.
df <- data.frame(
  x1 = grid$Sepal.Length,
  x2 = grid$Sepal.Width,
  class = class,
  max_prob = apply(grid_predictions, 1, max)
)

## Use BKP's internal class-plotting helper to match the contour-style figure
## layout used in the manuscript. The function is non-exported but documented
## in the package and is used here only for reproducing the manuscript figure.
p1 <- BKP:::my_2D_plot_fun_class(
  "class",
  "Predicted Classes",
  df,
  X_train,
  Y_train
)

## Plot the maximum predicted class probability. Here classification = FALSE
## indicates that max_prob should be visualized as a continuous probability
## surface rather than as a discrete class label.
p2 <- BKP:::my_2D_plot_fun_class(
  "max_prob",
  "Maximum Predicted Probability",
  df,
  X_train,
  Y_train,
  classification = FALSE
)

## Arrange the Gaussian process predicted class regions and maximum predicted
## probability surface side by side.
pdf("code/figure/ex7LGP.pdf", width = 13, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()