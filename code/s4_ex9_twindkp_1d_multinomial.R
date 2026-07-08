## -------------------------------------------------------------------------
## Section 4, Example 9: One-dimensional TwinDKP multinomial example
##
## This script reproduces Example 9 in Section 4 of the manuscript.
## It revisits the one-dimensional three-class multinomial response setting
## from Example 5, but increases the sample size and replaces the full DKP
## model with the scalable TwinDKP approximation.
##
## The goal is to illustrate that TwinDKP uses the same multinomial response
## format and downstream workflow as DKP while replacing full kernel
## aggregation with a twinning-based global-local approximation.
##
## Output figures are saved to code/figure/.
## -------------------------------------------------------------------------


## Load required packages.
## BKP provides the TwinDKP fitting and plotting methods.
## The tgp package is used for Latin hypercube sampling via lhs().
library(BKP)
library(tgp)

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
## Data generation
## -------------------------------------------------------------------------

## Set the random seed to make the simulated design points, multinomial trial
## sizes, and response counts reproducible.
set.seed(123)

## Define the true three-class probability function.
##
## This is the same class-probability construction used in Example 5.
## The first component is based on the logistic probability curve:
##
##   p1(x) = 1 / {1 + exp(-3x)}.
##
## The second component is based on the nonlinear oscillatory probability
## curve:
##
##   p2(x) = {1 + exp(-x^2) cos[10(1 - exp(-x))/(1 + exp(-x))]} / 2.
##
## The three class probabilities are then defined as
##
##   pi_1(x) = p1(x) / 2,
##   pi_2(x) = p2(x) / 2,
##   pi_3(x) = 1 - {p1(x) + p2(x)} / 2.
##
## This construction ensures that the three probabilities are nonnegative
## and sum to one at every input location.
true_pi_fun <- function(X) {
  p1 <- 1 / (1 + exp(-3 * X))
  p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
  return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
}

## Specify the number of training input locations.
## Compared with Example 5, the sample size is increased to demonstrate the
## scalable TwinDKP interface.
n <- 500

## Define the one-dimensional input domain [-2, 2].
Xbounds <- matrix(c(-2, 2), nrow = 1)

## Generate n input locations using a Latin hypercube design over [-2, 2].
X <- lhs(n = n, rect = Xbounds)

## Evaluate the true three-class probability vector at the training locations.
## The resulting object is an n by 3 matrix, where each row sums to one.
true_pi <- true_pi_fun(X)

## Generate heterogeneous multinomial trial sizes.
## Each input location has a total count sampled from 1, ..., 150.
m <- sample(150, n, replace = TRUE)

## Generate multinomial count vectors.
## For each input location x_i, the response vector Y_i is drawn from
## Multinomial(m_i, true_pi_i), where true_pi_i is the corresponding row
## of the true class-probability matrix.
Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))


## -------------------------------------------------------------------------
## TwinDKP model fitting
## -------------------------------------------------------------------------

## Fit the TwinDKP model to the simulated multinomial response data.
## The function uses a twinning-selected global subset and prediction-specific
## local neighbours. The response format is the same as in fit_DKP(), namely
## an n by q matrix of multinomial counts.
TwinDKP_model_1D <- fit_TwinDKP(X, Y, Xbounds = Xbounds)


## -------------------------------------------------------------------------
## Posterior summaries and true class-probability curves
## -------------------------------------------------------------------------

## Construct a dense prediction grid over the input domain.
## This grid is used to plot the true class-probability functions for
## comparison with the fitted TwinDKP posterior summaries.
Xnew <- matrix(seq(-2, 2, length = 100), ncol = 1)

## Evaluate the true class-probability functions on the prediction grid.
true_pi <- true_pi_fun(Xnew)

## Save the fitted TwinDKP posterior summaries and the true class-probability
## curves to a PDF file.
pdf("code/figure/ex9.pdf", width = 8, height = 8)

## The plot() method for a fitted one-dimensional TwinDKP object displays the
## posterior mean class-probability curves, pointwise credible intervals,
## observed class proportions, and the twinning-selected global subset.
plot(TwinDKP_model_1D)

## Add a separate panel showing the true class-probability functions.
plot(
  Xnew, true_pi[, 1],
  type = "l",
  col = "black",
  xlab = "x",
  ylab = "Probability",
  ylim = c(0, 1),
  main = "True Probability",
  lwd = 2
)

## Add the second and third true class-probability curves.
lines(Xnew, true_pi[, 2], col = "red", lwd = 2)
lines(Xnew, true_pi[, 3], col = "blue", lwd = 2)

## Add a legend identifying the three classes.
legend(
  "topright",
  bty = "n",
  legend = c("Class 1", "Class 2", "Class 3"),
  col = c("black", "red", "blue"),
  lty = 1,
  lwd = 2
)

## Close the PDF graphics device.
dev.off()