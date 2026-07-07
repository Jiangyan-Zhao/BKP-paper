rm(list = ls())
# setwd("XXXX") # Set your working directory here
library(BKP)
library(tgp)
library(gplite)
library(kernlab)
library(mlbench)
library(pROC)
library(gridExtra)
library(RiskMap)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
options(prompt = "R> ", continue = "+  ", width = 70,
        useFancyQuotes = FALSE)
# ============================================================== #
# ========================= BKP Examples ======================= #
# ============================================================== #
#-------------------------- 1D Example ---------------------------
#-------------------------- Example 1 ----------------------------
set.seed(123)
# Define true success probability function
true_pi_fun <- function(x) {
  1/(1+exp(-3*x))
}
# Data points
n <- 7
Xbounds <- matrix(c(-2,2), nrow = 1)
X <- lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(100, n, replace = TRUE)
y <- rbinom(n, size = m, prob = true_pi)

# Fit BKP model
BKP_model_1D_1 <- fit_BKP(X, y, m, Xbounds = Xbounds)

# Print BKP model
print(BKP_model_1D_1)

# New data points
Xnew = matrix(seq(-2, 2, length = 100), ncol = 1)
true_pi <- true_pi_fun(Xnew)

# Plot results
pdf("ex1.pdf", width = 6, height = 6)
plot(BKP_model_1D_1)
lines(Xnew,true_pi, col = "black", lwd = 2)
legend(x = -2.24,y = 0.89,
       legend = "True Probability",
       col = "black",
       lwd = 2,
       bty = "n",
       inset = 0.02)
dev.off()

# Simulate
sim <- simulate(BKP_model_1D_1, nsim = 3, Xnew = Xnew)
pdf("ex1_sim.pdf", width = 6, height = 6)
matplot(Xnew, sim$samples, type = "l", lty = 1, lwd = 1.5,
        col = rainbow(ncol(sim$samples)),
        xlab = "x", ylab = "Probability",
        main = "Simulated Probability Curves")
legend("topleft", legend = paste("sample", 1:ncol(sim$samples)),
       col = rainbow(ncol(sim$samples)), lty = 1, lwd = 2, bty = "n")
dev.off()

### Classification
set.seed(123)
# Data points
n <- 20
Xbounds <- matrix(c(-2,2), nrow = 1)
X <- lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- rep(1,n)
y <- as.numeric(true_pi > 0.5)

# Fit BKP model with r0 = 0.01
BKP_model_1D_1_class_1 <- fit_BKP(
  X, y, m, Xbounds = Xbounds,
  prior = "fixed", r0 = 0.01, loss = "log_loss")

# New data points
Xnew = matrix(seq(-2, 2, length = 100), ncol = 1)
true_pi <- true_pi_fun(Xnew)

# Plot results
pdf("ex1class001.pdf", width = 6, height = 6)
plot(BKP_model_1D_1_class_1)
lines(Xnew,true_pi, col = "black", lwd = 2)
dev.off()

# Fit BKP model with r0 = 2
BKP_model_1D_1_class_2 <- fit_BKP(
  X, y, m, Xbounds = Xbounds,
  prior = "fixed", r0 = 2, loss = "log_loss")

# Plot results
pdf("ex1class2.pdf", width = 6, height = 6)
plot(BKP_model_1D_1_class_2)
lines(Xnew,true_pi, col = "black", lwd = 2)
dev.off()

#-------------------------- Example 2 ----------------------------
set.seed(123)
# Define true success probability function
true_pi_fun <- function(x) {
  (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
}
# Data points
n <- 30
Xbounds <- matrix(c(-2,2), nrow = 1)
X <- lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(100, n, replace = TRUE)
y <- rbinom(n, size = m, prob = true_pi)

# Fit BKP model
BKP_model_1D_2 <- fit_BKP(X, y, m, Xbounds = Xbounds)

# New data points
Xnew = matrix(seq(-2, 2, length = 100), ncol=1)
true_pi <- true_pi_fun(Xnew)

# Plot results
pdf("ex2.pdf", width = 6, height = 6)
plot(BKP_model_1D_2)
lines(Xnew,true_pi, col = "black", lwd = 2)
legend(x = -2.24,y = 0.89,
       legend = "True Probability",
       col = "black",
       lwd = 2,
       bty = "n",
       inset = 0.02)
dev.off()

# Simulate
sim <- simulate(BKP_model_1D_2, nsim = 3, Xnew = Xnew)
pdf("ex2_sim.pdf", width = 6, height = 6)
matplot(Xnew, sim$samples, type = "l", lty = 1, lwd = 1.5,
        col = rainbow(ncol(sim$samples)),
        xlab = "x", ylab = "Probability",
        main = "Simulated Probability Curves")
legend("topleft", legend = paste("sample", 1:ncol(sim$samples)),
       col = rainbow(ncol(sim$samples)), lty = 1, lwd = 2, bty = "n")
dev.off()

#-------------------------- 2D Example ---------------------------
#-------------------------- Example 3 ----------------------------
set.seed(123)
# Define 2D latent function and probability transformation
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
  f <- log(a*b)
  f <- (f- m)/s
  return(pnorm(f))  # Transform to probability
}
# Data points
n <- 100
Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
X <- lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(100, n, replace = TRUE)
y <- rbinom(n, size = m, prob = true_pi)

# Fit BKP model
BKP_model_2D <- fit_BKP(X, y, m, Xbounds=Xbounds)

# Print BKP model
print(BKP_model_2D)

# Plot results
pdf("ex3.pdf", width = 9, height = 8)
plot(BKP_model_2D)
dev.off()

# Plot True distribution
pdf("ex3_true.pdf", width = 4.5, height = 4)
Xnew1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 200)
Xnew2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 200)
Xnew <- expand.grid(Xnew1 = Xnew1, Xnew2 = Xnew2)
true_pi <- true_pi_fun(Xnew)
df <- data.frame(x1 = Xnew$Xnew1, x2 = Xnew$Xnew2,
                 True = true_pi)
# Use BKP's internal 2D plotting helper to reproduce the contour-style plot.
# The function is non-exported; see ?my_2D_plot_fun for its internal help page.
print(BKP:::my_2D_plot_fun("True", title = "True Probability", data = df))
dev.off()

# Predict via BKP model
Xnew <- lhs(n = 10, rect = Xbounds)
predict(BKP_model_2D, Xnew)


# ============================================================
# Execution mode:
# run_elapsed_time = FALSE means using the pre-computed average
# elapsed time matrix, without re-running the time measurement
# experiments.
#
# Setting run_elapsed_time = TRUE will run the full benchmarking
# procedure for 5 methods × 5 n values × 20 repetitions.
# This can take many hours (potentially tens of hours) depending
# on your machine's performance.
#
# It is recommended to keep this FALSE unless you need to
# fully reproduce the timing experiments.
# ============================================================

run_elapsed_time = FALSE                # Flag to control whether to run timing experiments
n_vals <- c(200, 500, 1000, 2000, 5000) # Different sample sizes to benchmark
time_file <- "code/elapsed_time_avg.csv"     # File to save/read average elapsed time

if(run_elapsed_time){
  # Initialize matrices to store elapsed times for each method
  elapsed_time_given <- matrix(NA, nrow = 20, ncol = 5)
  elapsed_time_single_start <- matrix(NA, nrow = 20, ncol = 5)
  elapsed_time_multi_start <- matrix(NA, nrow = 20, ncol = 5)
  elapsed_time_given_gp <- matrix(NA, nrow = 20, ncol = 5)
  elapsed_time_optim_gp <- matrix(NA, nrow = 20, ncol = 5)

  # Loop over repetitions and sample sizes to benchmark
  for (i in 1:20) {
    for (j in 1:5) {
      set.seed((i - 1) * 5 + j)
      cat("i =", i, "j =", j, ":")
      n <- n_vals[j]
      Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
      X <- lhs(n = n, rect = Xbounds)
      true_pi <- true_pi_fun(X)
      m <- sample(100, n, replace = TRUE)
      y <- rbinom(n, size = m, prob = true_pi)

      # Measure elapsed time for BKP with given theta
      elapsed_time_given[i, j] <- system.time({
        fit <- fit_BKP(X, y, m, Xbounds = Xbounds, theta = 1)
      })["elapsed"]
      cat(" BKP_given =", elapsed_time_given[i, j])

      # Measure elapsed time for BKP with single start
      elapsed_time_single_start[i, j] <- system.time({
        fit <- fit_BKP(X, y, m, Xbounds = Xbounds, n_multi_start = 1)
      })["elapsed"]
      cat(", BKP_single_start =", elapsed_time_single_start[i, j])

      # Measure elapsed time for BKP with multiple starts
      elapsed_time_multi_start[i, j] <- system.time({
        fit <- fit_BKP(X, y, m, Xbounds = Xbounds)
      })["elapsed"]
      cat(", BKP_multi_start =", elapsed_time_multi_start[i, j])

      # Measure elapsed time for GP with given hyperparameters
      elapsed_time_given_gp[i, j] <- system.time({
        gp <- gp_init(cf = cf_sexp(), lik = lik_binomial())
        gp <- gp_fit(gp, X, y, trials = m)
      })["elapsed"]
      cat(", GP_given =", elapsed_time_given_gp[i, j])

      # Measure elapsed time for GP with optimized hyperparameters
      elapsed_time_optim_gp[i, j] <- system.time({
        gp <- gp_init(cf = cf_sexp(), lik = lik_binomial())
        gp <- gp_optim(gp, X, y, trials = m, verbose = FALSE)
      })["elapsed"]
      cat(", GP_optim =", elapsed_time_optim_gp[i, j], "\n")
    }
  }

  # ---- Compute average time for each n ----
  avg_time <- matrix(NA, nrow = 5, ncol = 5)
  colnames(avg_time) <- paste0("n=", n_vals)
  rownames(avg_time) <- c("given", "single_start", "multi_start", "given_gp", "optim_gp")
  avg_time[1, ] <- colMeans(elapsed_time_given)
  avg_time[2, ] <- colMeans(elapsed_time_single_start)
  avg_time[3, ] <- colMeans(elapsed_time_multi_start)
  avg_time[4, ] <- colMeans(elapsed_time_given_gp)
  avg_time[5, ] <- colMeans(elapsed_time_optim_gp)

  # Save the computed average time matrix to a CSV file
  time_file_new <- sub("\\.csv$", "_new.csv", time_file)
  write.csv(avg_time, file = time_file_new)
  cat("Timing results saved to:", time_file_new, "\n")
}else if(file.exists(time_file)){
  # If run_elapsed_time = FALSE, read the pre-computed CSV file
  avg_time <- as.matrix(read.csv(time_file, row.names = 1))
  rownames(avg_time) <- c("given", "single_start", "multi_start", "given_gp", "optim_gp")
  colnames(avg_time) <- c("n=200", "n=500", "n=1000", "n=2000", "n=5000")
} else {
  # Error if CSV file does not exist
  stop("Pre-computed timing file not found. Please run with run_elapsed_time = TRUE first.")
}

pdf("elapsed_time.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4,4,3,1))  # Adjust margins for better appearance
cols <- c("blue", "skyblue", "darkblue", "red", "orange")
pchs <- c(16, 17, 16, 15, 15)

# Left plot: given and given_gp
plot(n_vals, avg_time["given", ], log = "xy", type = "b", col = cols[1], pch = pchs[1], lwd = 2,
     xlab = "n", ylab = "Time (seconds, log scale)",
     main = "Fixed Hyperparameter",
     ylim = range(avg_time[c("given", "given_gp"), ]))
lines(n_vals, avg_time["given_gp", ], type = "b", col = cols[4], pch = pchs[4], lwd = 2)

# Complexity reference lines
ref_bkp <- avg_time["given", 1] * (n_vals / n_vals[1])^2
ref_gp  <- avg_time["given_gp", 1] * (n_vals / n_vals[1])^3
lines(n_vals, ref_bkp, col = "blue", lty = 2)
lines(n_vals, ref_gp, col = "red", lty = 4)

legend("topleft", bty = "n",
       legend = c("BKP", "LGP", "O(n^2)", "O(n^3)"),
       col = c(cols[c(1,4)], "blue", "red"),
       pch = c(pchs[c(1,4)], NA, NA),
       lty = c(1,1,2,4),
       lwd = c(2,2,1,1),
       bg = "white")

# Right plot: single_start, multi_start, optim_gp
plot(n_vals, avg_time["single_start", ], log = "xy", type = "b", col = cols[2], pch = pchs[2], lwd = 2,
     xlab = "n", ylab = "Time (seconds, log scale)",
     main = "Optimization-based Methods",
     ylim = range(avg_time[c("single_start", "multi_start", "optim_gp"), ]))
lines(n_vals, avg_time["multi_start", ], type = "b", col = cols[3], pch = pchs[3], lwd = 2)
lines(n_vals, avg_time["optim_gp", ], type = "b", col = cols[5], pch = pchs[5], lwd = 2)

# Complexity reference lines
ref_bkp2 <- avg_time["single_start", 1] * (n_vals / n_vals[1])^2
ref_gp2  <- avg_time["optim_gp", 1] * (n_vals / n_vals[1])^3
lines(n_vals, ref_bkp2, col = "skyblue", lty = 2)
lines(n_vals, ref_gp2, col = "orange", lty = 4)

legend("topleft", bty = "n",
       legend = c("BKP (single start)", "BKP (multi-start)", "LGP", "O(n^2)", "O(n^3)"),
       col = c(cols[c(2,3,5)], "skyblue", "orange"),
       pch = c(pchs[c(2,3,5)], NA, NA),
       lty = c(1,1,1,2,4),
       lwd = c(2,2,2,1,1),
       bg = "white")

# Restore default plotting parameters (optional)
par(mfrow = c(1,1))
dev.off()

#-------------------------- Classification ------------------------
#-------------------------- Example 4 ----------------------------
set.seed(123)
# Data points
n <- 120
data <- mlbench.spirals(n, cycles = 2, sd = 0.05)
train_indices <- sample(1:n, 0.7 * n)

X_train <- data$x[train_indices, ]
y_train <- as.numeric(data$classes[train_indices]) - 1 # Convert to 0/1 for BKP
X_test <- data$x[-train_indices, ]
y_test <- as.numeric(data$classes[-train_indices]) - 1 # Convert to 0/1 for BKP
m <- rep(1, length(train_indices))
Xbounds <- rbind(c(-1.7, 1.7), c(-1.7, 1.7))

# Fit BKP model
BKP_model_Class <- fit_BKP(
  X_train, y_train, m, Xbounds = Xbounds,
  prior = "fixed", r0 = 0.1, p0 = 0.5, loss = "log_loss")

# Predict via BKP
prediction <- predict(BKP_model_Class, X_test)

# Plot ROC
pdf("ex4roc.pdf", width = 6, height = 6)
roc_curve_BKP <- roc(y_test, prediction$mean)
plot(roc_curve_BKP, main = paste("ROC curve for BKP (AUC = ", round(auc(roc_curve_BKP), 3), ")", sep = ""), col = "#1c61b6")
dev.off()
pdf("ex4.pdf", width = 13, height = 6)
plot(BKP_model_Class)
dev.off()

# Fit LGP model with exponential approx
gp <- gp_init(cf = cf_sexp(), lik = lik_bernoulli())
gp <- gp_optim(gp, X_train, y_train, method = method_full(),
               approx = approx_ep(), verbose = FALSE)

# Predict via LGP
prediction_gp <- gp_pred(gp, as.matrix(X_test), transform = TRUE)

# Plot ROC
pdf("ex4rocLGP.pdf", width = 6, height = 6)
roc_curve_LGP <- roc(y_test, prediction_gp$mean)
plot(roc_curve_LGP, main = paste("ROC curve for LGP (AUC = ", round(auc(roc_curve_LGP), 3), ")", sep = ""), col = "#1c61b6")
dev.off()

# Generate the grid for predictions
x1_seq <- seq(-1.7, 1.7, length.out = 80)
x2_seq <- seq(-1.7, 1.7, length.out = 80)
grid <- expand.grid(x1 = x1_seq, x2 = x2_seq)

# Predict the grid for plot
prediction <- gp_pred(gp, as.matrix(grid), transform = TRUE, var = TRUE)
df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                 Mean = prediction$mean,
                 Variance = prediction$var)
p1 <- BKP:::my_2D_plot_fun("Mean", "Predictive Mean", df, X = X_train, y = y_train)
p2 <- BKP:::my_2D_plot_fun("Variance", "Predictive Variance", df, X = X_train, y = y_train)

# Plot results of LGP
pdf("ex4LGP.pdf", width = 13, height = 6)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

#-----------------------------------------------------------------
# ============================================================== #
# ========================= DKP Examples ======================= #
# ============================================================== #
#-------------------------- 1D Example ---------------------------
#-------------------------- Example 5 ----------------------------
set.seed(123)
# Define true class probability function (3-class)
true_pi_fun <- function(X) {
  p1 <- 1/(1+exp(-3*X))
  p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
  return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
}
# Data points
n <- 30
Xbounds <- matrix(c(-2, 2), nrow = 1)
X <- lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(150, n, replace = TRUE)
Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

# Fit DKP model
DKP_model_1D <- fit_DKP(X, Y, Xbounds = Xbounds)

# New data points
Xnew <- matrix(seq(-2,2, length = 100), ncol=1)
true_pi <- true_pi_fun(Xnew)

# Plot results
pdf("ex5.pdf", width = 8,height = 8)
plot(DKP_model_1D)
plot(Xnew, true_pi[, 1], type = "l", col = "black",
     xlab = "x", ylab = "Probability", ylim = c(0, 1),
     main = "True Probability", lwd = 2)
lines(Xnew, true_pi[, 2], col = "red", lwd = 2)
lines(Xnew, true_pi[, 3], col = "blue", lwd = 2)
legend("topright",
       bty = "n",
       legend = c("Class 1", "Class 2", "Class 3"),
       col = c("black", "red", "blue"), lty = 1, lwd = 2)
dev.off()

#-------------------------- 2D Example ---------------------------
#-------------------------- Example 6 ----------------------------
set.seed(123)
# Define 2D latent function and probability transformation (3-class)
true_pi_fun <- function(X){
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
  p1 <- pnorm(f)  # Transform to probability
  p2 <- sin(pi * X[,1]) * sin(pi * X[,2])
  return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
}
# Data points
n <- 100
Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
X <- lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(150, n, replace = TRUE)
Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

# Fit DKP model
DKP_model_2D <- fit_DKP(X, Y, Xbounds=Xbounds)

# Plot results
pdf(file = "ex6_class%d.pdf", width = 9, height = 8, onefile = FALSE)
plot(DKP_model_2D)
dev.off()

# New data points
Xnew1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 200)
Xnew2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 200)
Xnew <- expand.grid(Xnew1 = Xnew1, Xnew2 = Xnew2)
true_pi <- true_pi_fun(Xnew)
df <- data.frame(x1 = Xnew$Xnew1, x2 = Xnew$Xnew2,
                 True1 = true_pi[,1],
                 True2 = true_pi[,2],
                 True3 = true_pi[,3])
# Plot True distribution
pdf("ex6_class1_true.pdf", width = 4.5, height = 4)
print(BKP:::my_2D_plot_fun("True1", title = "True Probability", data = df))
dev.off()

pdf("ex6_class2_true.pdf", width = 4.5, height = 4)
print(BKP:::my_2D_plot_fun("True2", title = "True Probability", data = df))
dev.off()

pdf("ex6_class3_true.pdf", width = 4.5, height = 4)
print(BKP:::my_2D_plot_fun("True3", title = "True Probability", data = df))
dev.off()


#-------------------------- Classification ------------------------
#-------------------------- Example 7 ----------------------------
set.seed(123)
data(iris)
X <- as.matrix(iris[, 1:2])
Xbounds <- rbind(c(4.2, 8), c(1.9, 4.5))
labels <- iris$Species
Y <- model.matrix(~ labels - 1) # expand factors to a set of dummy variables

train_indices <- sample(1:nrow(iris), 0.7 * nrow(iris))
X_train <- X[train_indices, ]
Y_train <- Y[train_indices, ]
labels_train <- labels[train_indices]

X_test <- X[-train_indices, ]
Y_test <- Y[-train_indices, ]
labels_test <- labels[-train_indices]

# Fit DKP model
DKP_model_Class <- fit_DKP(
  X_train, Y_train, Xbounds = Xbounds, loss = "log_loss",
  prior = "fixed", r0 = 0.01, p0 = rep(1/3, 3))
pdf("ex7.pdf", width = 13, height = 6)
plot(DKP_model_Class)
dev.off()

# Test
dkp_pred_probs <- predict(DKP_model_Class, X_test)$mean
class_levels <- levels(labels_test)
colnames(dkp_pred_probs) <- class_levels
multiclass_roc_dkp <- multiclass.roc(labels_test, dkp_pred_probs)

all_rocs <- list()
for (class_name in class_levels) {
  # convert to binary label
  labels_binary <- ifelse(labels_test == class_name, 1, 0)
  probabilities <- dkp_pred_probs[, class_name]
  roc_curve <- roc(labels_binary, probabilities)
  all_rocs[[class_name]] <- roc_curve
}

# Plot ROC curve
pdf("ex7roc.pdf", width = 6, height = 6)
plot(all_rocs[[1]], col = "blue", lwd = 2, lty = 1,
     main = paste("One-vs-Rest ROC curve for BKP (AUC =", round(auc(multiclass_roc_dkp), 3), ")", sep = ""))
lines(all_rocs[[2]], col = "red", lwd = 2, lty = 2)
lines(all_rocs[[3]], col = "black", lwd = 2, lty = 4)
legend("bottomright",
       legend = class_levels,
       col = c("blue", "red", "black"),
       lwd = 2,
       lty = c(1,2,4),
       cex = 1.2)
dev.off()

### LGP model
iris_data <- data.frame(
  Sepal.Length = iris$Sepal.Length,
  Sepal.Width = iris$Sepal.Width,
  Species = iris$Species
)

iris_train <- iris_data[train_indices, ]
iris_test <- iris_data[-train_indices, ]

# Fit LGP model
gausspr_model <- gausspr(Species ~ ., data = iris_train,
                         kernel = "rbfdot", kpar = "automatic")

# Test LGP
lgp_pred_probs <- predict(gausspr_model, newdata = iris_test,
                          type = "probabilities")
multiclass_roc_lgp <- multiclass.roc(iris_test$Species, lgp_pred_probs)

# Plot ROC
all_rocs <- list()
for (class_name in class_levels) {
  # binary label
  labels_binary <- ifelse(labels_test == class_name, 1, 0)
  probabilities <- lgp_pred_probs[, class_name]
  roc_curve <- roc(labels_binary, probabilities)
  all_rocs[[class_name]] <- roc_curve
}

pdf("ex7rocLGP.pdf", width = 6, height = 6)
plot(all_rocs[[1]], col = "blue", lwd = 2, lty = 1,
     main = paste("One-vs-Rest ROC curve for LGP (AUC =", round(auc(multiclass_roc_lgp), 3), ")", sep = ""))
lines(all_rocs[[2]], col = "red", lwd = 2, lty = 2)
lines(all_rocs[[3]], col = "black", lwd = 2, lty = 4)
legend("bottomright",
       legend = class_levels,
       col = c("blue", "red", "black"),
       lty = c(1,2,4),
       lwd = 2,
       cex = 1.2)
dev.off()

# New data points for plot
grid <- expand.grid(
  Sepal.Length = seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 80),
  Sepal.Width = seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 80)
)
grid_predictions <- predict(gausspr_model, newdata = grid, type = "prob")
class <- max.col(grid_predictions)
df<- data.frame(x1 = grid$Sepal.Length,x2=grid$Sepal.Width,
                class = class,
                max_prob = apply(grid_predictions, 1, max))
p1 <- BKP:::my_2D_plot_fun_class("class", "Predicted Classes", df, X_train, Y_train)
p2 <- BKP:::my_2D_plot_fun_class("max_prob", "Maximum Predicted Probability", df, X_train, Y_train, classification = FALSE)
pdf("ex7LGP.pdf", width = 13, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()



#-------------------------- Example 8 ----------------------------
set.seed(123)
# True probability function from Example 2
true_pi_fun <- function(x) {
  x <- as.numeric(x)
  0.5 * (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x))))
}

# Simulate binomial data
n <- 500
Xbounds <- matrix(c(-2, 2), nrow = 1)

X <- lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)

m <- sample(100, n, replace = TRUE)
y <- rbinom(n, size = m, prob = true_pi)

## Fit TwinBKP model
TwinBKP_model_1D_2 <- fit_TwinBKP(X, y, m, Xbounds = Xbounds)

## Print fitted model
print(TwinBKP_model_1D_2)

## Inspect realized TwinBKP control settings
TwinBKP_model_1D_2$control[c("g", "l", "twins")]

## Prediction grid
Xnew <- matrix(seq(-2, 2, length = 100), ncol = 1)

## Posterior prediction
TwinBKP_pred_1D_2 <- predict(TwinBKP_model_1D_2, Xnew = Xnew)
print(TwinBKP_pred_1D_2)


## Plot fitted TwinBKP posterior summaries
pdf("ex8.pdf", width = 9, height = 6)
plot(TwinBKP_model_1D_2)
true_pi <- true_pi_fun(Xnew)
lines(Xnew,true_pi, col = "black", lwd = 2)
legend(x = -2.22,y = 0.85,
       legend = "True Probability",
       col = "black",
       lwd = 2,
       bty = "n",
       inset = 0.02)
dev.off()


#-----------------------------------------------------------------
# ============================================================== #
# ========================= Real Example ======================= #
# ============================================================== #

# Load the data
data("Loaloa")
# Extract input variables (X), response variable (y), and trial counts (m)
X <- as.matrix(Loaloa[, 1:2])
rownames(X) <- NULL
y <- Loaloa$npos
m <- Loaloa$ntot

# Randomly split into training (70%) and testing (30%) sets
set.seed(123)
train_idx <- sample(1:nrow(Loaloa), 0.7 * nrow(Loaloa))
X_train <- X[train_idx, ]
y_train <- y[train_idx]
m_train <- m[train_idx]
X_test <- X[-train_idx, ]
y_test <- y[-train_idx]
m_test <- m[-train_idx]


p <- y / m # Infection rate
Xbounds <- matrix(c(7.8, 15.3, 3.1, 7.0), ncol = 2, byrow = TRUE)
df <- data.frame(
  lon = X[,1], lat = X[,2], p = p, n = m,
  set = ifelse(1:nrow(Loaloa) %in% train_idx, "Train", "Test")
)

# Obtain African map data (sf object)
africa <- ne_countries(continent = "africa", scale = "medium", returnclass = "sf")

pmap <- ggplot() +
  geom_sf(data = africa, fill = "gray95", color = "gray60") +
  geom_point(data = df,
             aes(x = lon, y = lat, color = p, size = n, shape = set),
             alpha = 0.8) +
  scale_color_viridis_c(name = "Proportion y/m") +
  scale_size_continuous(name = "Trial count m") +
  scale_shape_manual(name = "Dataset", values = c("Train" = 16, "Test" = 17)) +
  coord_sf(xlim = Xbounds[1, ], ylim = Xbounds[2, ], expand = FALSE) +
  theme_minimal() +
  # theme(aspect.ratio = 1) +
  labs(title = "Loaloa infection proportion (y/m): Train vs Test",
       x = "Longitude", y = "Latitude")
ggsave("Loaloa_map.pdf", plot = pmap, width = 10, height = 5)


# Model fitting
Loaloa_bkp_model <- fit_BKP(
  X_train, y_train, m_train, Xbounds, loss = "brier",
  prior = "adaptive", r0 = mean(m_train))

# Print the model summary
summary(Loaloa_bkp_model)

# Plot the fitted model
pdf("Loaloa_bkp.pdf", width = 12, height = 12)
plot(Loaloa_bkp_model)
dev.off()

# LGP model for comparison
Loaloa_gp_model <- gp_init(cf = cf_sexp(), lik = lik_binomial())
Loaloa_gp_model <- gp_optim(
  Loaloa_gp_model, X_train, y_train, trials = m_train,
  method = method_full(), approx = approx_ep(), verbose = FALSE)

# Plot the fitted LGP model
pdf("Loaloa_gp.pdf", width = 12, height = 12)
# Generate the grid for predictions
x1_seq <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 80)
x2_seq <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 80)
grid <- expand.grid(x1 = x1_seq, x2 = x2_seq)

# Predict the grid for plot
prediction <- gp_pred(Loaloa_gp_model, as.matrix(grid),
                      transform = TRUE, var = TRUE, quantiles = c(0.025,0.975))
df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                 Mean = prediction$mean,
                 Upper = prediction$quantiles[,2],
                 Lower = prediction$quantiles[,1],
                 Variance = prediction$var)
dims <- c(1,2)
# Width = prediction$upper - prediction$lower)
p1 <- BKP:::my_2D_plot_fun("Mean", "Predictive Mean", df, dims= dims)
p3 <- BKP:::my_2D_plot_fun("Variance", "Predictive Variance", df, dims= dims)
p2 <- BKP:::my_2D_plot_fun("Upper", "95% CI Upper", df, dims= dims)
p4 <- BKP:::my_2D_plot_fun("Lower", "95% CI Lower", df, dims= dims)
# Arrange into 2×2 layout
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

# Predict the mean infection rate for the test data
predict_on_test_bkp <- predict(Loaloa_bkp_model, Xnew = X_test)
predict_on_test_gp <- gp_pred(Loaloa_gp_model, xnew = X_test,
                              transform = TRUE, var = TRUE)

# Empirical success rate for the test data
pi_tilde_test <- y_test / m_test

# Mean Squared Error (Brier Score)
mse_bkp <- mean((predict_on_test_bkp$mean - pi_tilde_test)^2)
mse_gp <- mean((predict_on_test_gp$mean - pi_tilde_test)^2)

cat("Mean Squared Error (BKP) on the test data:", mse_bkp, "\n")
cat("Mean Squared Error (LGP) on the test data:", mse_gp, "\n")













# =========================
# 0) Packages
# =========================
library(BKP)
library(ggplot2)
library(sf)
library(viridis)
library(patchwork)
library(rnaturalearth)

# （可选但推荐）避免 s2 几何带来的 st_intersects/st_union 小概率报错
sf::sf_use_s2(FALSE)

# =========================
# 1) Load & split data
# =========================
data("Loaloa")

X <- as.matrix(Loaloa[, 1:2])
rownames(X) <- NULL
y <- Loaloa$npos
m <- Loaloa$ntot

set.seed(123)
train_idx <- sample(1:nrow(Loaloa), 0.7 * nrow(Loaloa))

X_train <- X[train_idx, ]
y_train <- y[train_idx]
m_train <- m[train_idx]

X_test <- X[-train_idx, ]
y_test <- y[-train_idx]
m_test <- m[-train_idx]

p <- y / m
Xbounds <- matrix(c(7.8, 15.3, 3.1, 7.0), ncol = 2, byrow = TRUE)

# 观测点（不要再用 df 这个名字，避免后面被覆盖）
df_pts <- data.frame(
  lon = X[,1], lat = X[,2], p = p, n = m,
  set = ifelse(1:nrow(Loaloa) %in% train_idx, "Train", "Test")
)

# =========================
# 2) Africa sf + crop + land polygon (mask)
# =========================
africa <- ne_countries(continent = "africa", scale = "medium", returnclass = "sf")

africa_crop <- st_crop(
  africa,
  xmin = Xbounds[1,1], xmax = Xbounds[1,2],
  ymin = Xbounds[2,1], ymax = Xbounds[2,2]
)

africa_crop <- st_make_valid(africa_crop)
africa_crop <- st_transform(africa_crop, 4326)

# 合并成一个“陆地整体”多边形（去掉内部国界）
land_poly <- st_union(africa_crop)

# =========================
# 3) Fit models
# =========================
Loaloa_bkp_model <- fit_BKP(
  X_train, y_train, m_train, Xbounds#,
  #loss = "brier", prior = "adaptive", r0 = mean(m_train)
)

library(gplite)  # 你用的 LGP / gp_* 函数来自这里（如果不是，请按你实际包名调整）
Loaloa_gp_model <- gp_init(cf = cf_sexp(), lik = lik_binomial())
Loaloa_gp_model <- gp_optim(
  Loaloa_gp_model, X_train, y_train, trials = m_train,
  method = method_full(), approx = approx_ep(), verbose = FALSE
)

# =========================
# 4) Grid + inside_land
# =========================
lon_seq <- seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 140)
lat_seq <- seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 140)
grid_df <- expand.grid(lon = lon_seq, lat = lat_seq)
X_grid  <- as.matrix(grid_df)

grid_sf <- st_as_sf(grid_df, coords = c("lon","lat"), crs = 4326, remove = FALSE)
inside_land <- lengths(st_intersects(grid_sf, land_poly)) > 0

# =========================
# 5) Predict BKP / LGP on grid
# =========================
pred_bkp <- predict(Loaloa_bkp_model, Xnew = X_grid, CI_level = 0.95)
surf_bkp <- data.frame(
  lon = grid_df$lon,
  lat = grid_df$lat,
  Mean     = pred_bkp$mean,
  Variance = pred_bkp$variance,
  Lower    = pred_bkp$lower,
  Upper    = pred_bkp$upper
)

pred_gp <- gp_pred(
  Loaloa_gp_model, X_grid,
  transform = TRUE, var = TRUE,
  quantiles = c(0.025, 0.975)
)
surf_gp <- data.frame(
  lon = grid_df$lon,
  lat = grid_df$lat,
  Mean     = pred_gp$mean,
  Variance = pred_gp$var,
  Lower    = pred_gp$quantiles[,1],
  Upper    = pred_gp$quantiles[,2]
)

# =========================
# 6) Mask: keep predictions ONLY on land
# =========================
mask_cols <- c("Mean","Variance","Lower","Upper")
surf_bkp[!inside_land, mask_cols] <- NA
surf_gp[!inside_land,  mask_cols] <- NA

# =========================
# 7) Plot helper: geom_sf outline + raster (NA transparent)
# =========================
plot_metric_land <- function(surface_df, metric, title,
                             africa_crop, pts_df, Xbounds,
                             show_points = TRUE,
                             show_contour = FALSE,
                             contour_bins = 12) {

  g <- ggplot() +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "gray98", color = NA)
    ) +
    # 预测面（海里是 NA -> 透明）
    geom_raster(
      data = surface_df,
      aes(x = lon, y = lat, fill = .data[[metric]]),
      interpolate = TRUE,
      alpha = 0.95
    ) +
    scale_fill_viridis_c(name = metric, na.value = "transparent") +

    # 可选等高线（只会在陆地的有值区域画）
    {if (show_contour)
      geom_contour(
        data = surface_df,
        aes(x = lon, y = lat, z = .data[[metric]]),
        bins = contour_bins,
        linewidth = 0.25,
        alpha = 0.7
      )
    } +

    # 海岸线/国界线：放最上层，保证清晰
    geom_sf(
      data = africa_crop,
      fill = NA,
      color = "gray10",
      linewidth = 0.45
    ) +

    coord_sf(xlim = Xbounds[1, ], ylim = Xbounds[2, ], expand = FALSE) +
    labs(title = title, x = "Longitude", y = "Latitude")

  if (show_points) {
    g <- g +
      geom_point(
        data = pts_df,
        aes(x = lon, y = lat, shape = set),
        color = "black", size = 2, alpha = 0.85
      ) +
      scale_shape_manual(values = c("Train" = 16, "Test" = 17), name = "Dataset")
  }

  g
}

# =========================
# 8) BKP plots (2×2) + save
# =========================
p_bkp_mean  <- plot_metric_land(surf_bkp, "Mean",     "BKP: Predictive Mean",
                                africa_crop, df_pts, Xbounds, show_contour = TRUE)
p_bkp_upper <- plot_metric_land(surf_bkp, "Upper",    "BKP: 95% CI Upper",
                                africa_crop, df_pts, Xbounds)
p_bkp_var   <- plot_metric_land(surf_bkp, "Variance", "BKP: Predictive Variance",
                                africa_crop, df_pts, Xbounds)
p_bkp_lower <- plot_metric_land(surf_bkp, "Lower",    "BKP: 95% CI Lower",
                                africa_crop, df_pts, Xbounds)

p_bkp_all <- (p_bkp_mean | p_bkp_upper) / (p_bkp_var | p_bkp_lower)
ggsave("Loaloa_bkp_sf_land.pdf", plot = p_bkp_all, width = 16, height = 8)

# =========================
# 9) LGP plots (2×2) + save
# =========================
p_gp_mean  <- plot_metric_land(surf_gp, "Mean",     "LGP (EP): Predictive Mean",
                               africa_crop, df_pts, Xbounds, show_contour = TRUE)
p_gp_upper <- plot_metric_land(surf_gp, "Upper",    "LGP (EP): 95% CI Upper",
                               africa_crop, df_pts, Xbounds)
p_gp_var   <- plot_metric_land(surf_gp, "Variance", "LGP (EP): Predictive Variance",
                               africa_crop, df_pts, Xbounds)
p_gp_lower <- plot_metric_land(surf_gp, "Lower",    "LGP (EP): 95% CI Lower",
                               africa_crop, df_pts, Xbounds)

p_gp_all <- (p_gp_mean | p_gp_upper) / (p_gp_var | p_gp_lower)
ggsave("Loaloa_gp_sf_land.pdf", plot = p_gp_all, width = 16, height = 8)
