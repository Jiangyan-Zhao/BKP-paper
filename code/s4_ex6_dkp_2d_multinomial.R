## -------------------------------------------------------------------------
## Section 4, Example 6: Two-dimensional DKP multinomial example
##
## This script reproduces Example 6 in Section 4 of the manuscript.
## It illustrates the Dirichlet Kernel Process (DKP) model for a
## two-dimensional three-class multinomial response problem.
##
## The first class probability is constructed from a rescaled
## Goldstein--Price-type surface, the second class probability is based on a
## smooth sinusoidal surface, and the third class probability is defined as
## the remaining probability mass. At each input location, a multinomial count
## vector is generated and the DKP model is fitted to estimate the
## class-probability surfaces.
##
## Output figures are saved to code/figure/.
## -------------------------------------------------------------------------


## Load required packages.
## BKP provides the Dirichlet Kernel Process fitting and plotting methods.
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

## Define the two-dimensional true three-class probability function.
##
## The first component is based on a rescaled Goldstein--Price-type latent
## function. The input X is assumed to lie in [0, 1]^2. The variables x1 and
## x2 are first mapped from [0, 1] to [-2, 2], and then used to evaluate the
## latent surface. The transformed value pnorm(f) gives a probability-scale
## function p1(x).
##
## The second component is a smooth sinusoidal function,
##
##   p2(x) = sin(pi x1) sin(pi x2).
##
## The three class probabilities are then defined as
##
##   pi_1(x) = p1(x) / 2,
##   pi_2(x) = p2(x) / 2,
##   pi_3(x) = 1 - {p1(x) + p2(x)} / 2.
##
## This construction ensures that the class probabilities are nonnegative
## and sum to one at every input location.
true_pi_fun <- function(X) {
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  m <- 8.6928
  s <- 2.4269
  x1 <- 4*X[,1]- 2
  x2 <- 4*X[,2]- 2
  a <- 1 + (x1 + x2 + 1)^2 *
    (19- 14*x1 + 3*x1^2- 14*x2 + 6*x1*x2 + 3*x2^2)
  b <- 30 + (2*x1- 3*x2)^2 *
    (18- 32*x1 + 12*x1^2 + 48*x2- 36*x1*x2 + 27*x2^2)
  f <- (log(a*b)- m)/s
  p1 <- pnorm(f)
  p2 <- sin(pi * X[,1]) * sin(pi * X[,2])
  return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
}

## Specify the number of training input locations.
n <- 100

## Define the two-dimensional input domain [0, 1]^2.
## Each row of Xbounds gives the lower and upper bounds of one input dimension.
Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)

## Generate n input locations using a Latin hypercube design over [0, 1]^2.
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
## DKP model fitting
## -------------------------------------------------------------------------

## Fit the standard DKP model to the simulated two-dimensional multinomial
## response data. Since theta is not supplied, the kernel length-scale
## parameter is selected internally by leave-one-out cross-validation using
## the default loss.
DKP_model_2D <- fit_DKP(X, Y, Xbounds = Xbounds)


## -------------------------------------------------------------------------
## DKP posterior summary plots
## -------------------------------------------------------------------------

## Save the DKP posterior summary plots to separate PDF files.
## For a two-dimensional three-class DKP model, plot() produces class-specific
## posterior summary figures. With onefile = FALSE and the "%d" pattern in
## the file name, R writes the class-specific figures to separate files.
pdf(file = "code/figure/ex6_class%d.pdf", width = 9, height = 8, onefile = FALSE)
plot(DKP_model_2D)
dev.off()


## -------------------------------------------------------------------------
## True class-probability surfaces
## -------------------------------------------------------------------------

## Construct a regular 200 by 200 prediction grid over [0, 1]^2.
## This grid is used to visualize the true class-probability surfaces for
## comparison with the fitted DKP posterior summaries.
Xnew1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 200)
Xnew2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 200)
Xnew <- expand.grid(Xnew1 = Xnew1, Xnew2 = Xnew2)

## Evaluate the true class-probability functions on the prediction grid.
true_pi <- true_pi_fun(Xnew)

## Store the grid coordinates and true class probabilities in a data frame
## compatible with the BKP two-dimensional plotting helper.
df <- data.frame(
  x1 = Xnew$Xnew1,
  x2 = Xnew$Xnew2,
  True1 = true_pi[,1],
  True2 = true_pi[,2],
  True3 = true_pi[,3]
)

## Plot the true probability surface for Class 1.
## The internal plotting helper is used here only to match the contour-style
## figure format used elsewhere in the manuscript.
pdf("code/figure/ex6_class1_true.pdf", width = 4.5, height = 4)
print(BKP:::my_2D_plot_fun("True1", title = "True Probability", data = df))
dev.off()

## Plot the true probability surface for Class 2.
pdf("code/figure/ex6_class2_true.pdf", width = 4.5, height = 4)
print(BKP:::my_2D_plot_fun("True2", title = "True Probability", data = df))
dev.off()

## Plot the true probability surface for Class 3.
pdf("code/figure/ex6_class3_true.pdf", width = 4.5, height = 4)
print(BKP:::my_2D_plot_fun("True3", title = "True Probability", data = df))
dev.off()