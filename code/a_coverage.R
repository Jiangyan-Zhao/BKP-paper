# ======================================================================
# Coverage Rate for Example 2
# BKP (none) vs BKP (shepard) vs LGP — 2 sample sizes
# ======================================================================

library(BKP)
library(gplite)
library(tgp)

# ======================== Tunable parameters ===========================
n_reps   <- 2000
ngrid    <- 2000
CI_level <- 0.95
alpha_ci <- (1 - CI_level) / 2
n_train_values <- c(30, 100)

# ======================== Define Example 2 =============================
set.seed(123)
Xbounds <- matrix(c(-2, 2), nrow = 1)

true_pi <- function(x) {
  (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
}

generate_data <- function(n, Xbounds) {
  X <- lhs(n = n, rect = Xbounds)
  pi <- true_pi(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = pi)
  list(X = X, y = y, m = m)
}

# ======================================================================
# Metrics — Coverage rate only
# ======================================================================

# Compute coverage from prediction vs truth
compute_coverage <- function(lower, upper, truth) {
  mean((truth >= lower) & (truth <= upper), na.rm = TRUE)
}

# ======================================================================
# Prep grid
# ======================================================================
Xgrid <- matrix(seq(Xbounds[1,1], Xbounds[1,2], length.out = ngrid), ncol = 1)
true_pi_grid <- true_pi(Xgrid)

results <- list()

# Store plot data for combined figure
plot_preds   <- list()
plot_methods <- c("BKP", "BKP-ESS", "LGP")

# ======================================================================
# Loop over sample sizes
# ======================================================================
for (ntr in n_train_values) {
  config_name <- sprintf("Example 2, n=%d", ntr)
  cat("\n============================================================\n")
  cat(sprintf("Sample size n = %d\n", ntr))
  cat("============================================================\n\n")

  cat("  Fitting methods for plot ...\n")
  dat_plot <- generate_data(ntr, Xbounds)
  p0 <- mean(dat_plot$y / dat_plot$m)

  fit_none <- fit_BKP(dat_plot$X, dat_plot$y, dat_plot$m,
                      Xbounds = Xbounds, r0 = 2, prior = "noninformative",
                      p0 = p0, loss = "brier", ess = "none")
  fit_shep <- fit_BKP(dat_plot$X, dat_plot$y, dat_plot$m,
                      Xbounds = Xbounds, r0 = 2, prior = "noninformative",
                      p0 = p0, loss = "brier", ess = "shepard")

  gp_plot <- gp_init(cf = cf_sexp(), lik = lik_binomial())
  gp_plot <- suppressWarnings(gp_optim(gp_plot, dat_plot$X, dat_plot$y,
                                        trials = dat_plot$m, verbose = FALSE))

  pred_none <- predict(fit_none, Xnew = Xgrid, CI_level = CI_level)
  pred_shep <- predict(fit_shep, Xnew = Xgrid, CI_level = CI_level)
  pred_gp   <- gp_pred(gp_plot, Xgrid,
                       quantiles = c(alpha_ci, 1 - alpha_ci), transform = TRUE)

  # Save plot data
  obs_y <- dat_plot$y / dat_plot$m
  plot_preds[[config_name]] <- list(
    X       = dat_plot$X,
    obs_y   = obs_y,
    none    = list(mean = pred_none$mean,  lower = pred_none$lower,  upper = pred_none$upper),
    shepard = list(mean = pred_shep$mean,  lower = pred_shep$lower,  upper = pred_shep$upper),
    lgp     = list(mean = pred_gp$mean    , lower = pred_gp$quantiles[,1], upper = pred_gp$quantiles[,2])
  )

  cat(sprintf("  Running %d repetitions ...\n", n_reps))

  # Storage — coverage only
  bn_grid_cov <- bn_tr_cov <- numeric(n_reps)
  bs_grid_cov <- bs_tr_cov <- numeric(n_reps)
  lg_grid_cov <- lg_tr_cov <- numeric(n_reps)

  for (rep in seq_len(n_reps)) {
    dat <- generate_data(ntr, Xbounds)
    true_pi_tr <- true_pi(dat$X)

    # ---- BKP ----
    p0_rep <- mean(dat$y / dat$m)
    m <- fit_BKP(dat$X, dat$y, dat$m, Xbounds = Xbounds, r0 = 2,
                 prior = "noninformative", p0 = p0_rep, loss = "brier",
                 ess = "none")
    pg <- predict(m, Xnew = Xgrid, CI_level = CI_level)
    pt <- predict(m, Xnew = NULL,  CI_level = CI_level)
    bn_grid_cov[rep] <- compute_coverage(pg$lower, pg$upper, true_pi_grid)
    bn_tr_cov[rep]   <- compute_coverage(pt$lower, pt$upper, true_pi_tr)

    # ---- BKP-ESS ----
    m <- fit_BKP(dat$X, dat$y, dat$m, Xbounds = Xbounds, r0 = 2,
                 prior = "noninformative", p0 = p0_rep, loss = "brier",
                 ess = "shepard")
    pg <- predict(m, Xnew = Xgrid, CI_level = CI_level)
    pt <- predict(m, Xnew = NULL,  CI_level = CI_level)
    bs_grid_cov[rep] <- compute_coverage(pg$lower, pg$upper, true_pi_grid)
    bs_tr_cov[rep]   <- compute_coverage(pt$lower, pt$upper, true_pi_tr)

    # ---- LGP ----
    g <- gp_init(cf = cf_sexp(), lik = lik_binomial())
    g <- suppressWarnings(gp_optim(g, dat$X, dat$y, trials = dat$m, verbose = FALSE))
    pg <- gp_pred(g, Xgrid, quantiles = c(alpha_ci, 1 - alpha_ci), transform = TRUE)
    pt <- gp_pred(g, dat$X, quantiles = c(alpha_ci, 1 - alpha_ci), transform = TRUE)
    lg_grid_cov[rep] <- compute_coverage(pg$quantiles[,1], pg$quantiles[,2], true_pi_grid)
    lg_tr_cov[rep]   <- compute_coverage(pt$quantiles[,1], pt$quantiles[,2], true_pi_tr)
  }

  # Summarise
  sumf <- function(x) c(mu = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))

  sg_bn <- sumf(bn_grid_cov)
  sg_bs <- sumf(bs_grid_cov)
  sg_lg <- sumf(lg_grid_cov)
  st_bn <- sumf(bn_tr_cov)
  st_bs <- sumf(bs_tr_cov)
  st_lg <- sumf(lg_tr_cov)

  # Print
  cat(sprintf("\n  %-16s  %12s\n", "", "Coverage"))
  hr <- function(lbl, s) {
    cat(sprintf("  %-16s  %5.1f%%(SD=%4.1f)\n",
                lbl, 100 * s["mu"], 100 * s["sd"]))
  }
  cat("\n  --- Grid ---\n")
  hr("BKP",     sg_bn)
  hr("BKP-ESS", sg_bs)
  hr("LGP",     sg_lg)
  cat("\n  --- Training ---\n")
  hr("BKP",     st_bn)
  hr("BKP-ESS", st_bs)
  hr("LGP",     st_lg)

  # Store
  results[[config_name]] <- list(
    n = ntr,
    grid_bn  = sg_bn,
    grid_bs  = sg_bs,
    grid_lg  = sg_lg,
    train_bn = st_bn,
    train_bs = st_bs,
    train_lg = st_lg
  )
}

# ======================================================================
# Combined 3×2 figure with subfigure labels (a)–(f)
# ======================================================================
pdf("code/figure/ex2_combined.pdf", width = 8, height = 12)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

sub_labels <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
sub_idx <- 1

for (mi in seq_along(plot_methods)) {
  method <- plot_methods[mi]
  key <- switch(mi, "1"="none", "2"="shepard", "3"="lgp")
  for (ntr in n_train_values) {
    config_name <- sprintf("Example 2, n=%d", ntr)
    pd <- plot_preds[[config_name]]
    mean  <- pd[[key]]$mean
    lower <- pd[[key]]$lower
    upper <- pd[[key]]$upper
    obs_y <- pd$obs_y
    ylim  <- c(max(0, min(c(lower, obs_y)) * 0.9), min(1, max(c(upper, obs_y)) * 1.1))

    plot(Xgrid, mean, type = "n",
         xlab = "x", ylab = "Probability",
         main = sprintf("%s %s, N=%d", sub_labels[sub_idx], method, ntr),
         ylim = ylim)
    polygon(c(Xgrid, rev(Xgrid)), c(lower, rev(upper)),
            col = "lightgrey", border = NA)
    lines(Xgrid, mean, col = "blue", lwd = 2)
    points(pd$X, obs_y, pch = 20, col = "red")
    lines(Xgrid, true_pi_grid, col = "black", lwd = 2)
    legend("topleft", bty = "n", cex = 1,
           legend = c("Estimated Probability", "95% Credible Interval",
                      "Observed Proportions", "True Probability"),
           col = c("blue", "lightgrey", "red", "black"),
           lwd = c(2, 8, NA, 2), pch = c(NA, NA, 20, NA),
           lty = c(1, 1, NA, 1))
    sub_idx <- sub_idx + 1
  }
}
dev.off()
cat("  Combined figure saved: ex2_combined.pdf\n\n")

# ======================================================================
# Summary table — Coverage only
# ======================================================================
cat("\n================================================================\n")
cat("Summary — Coverage Rate (nominal 95%)\n")
cat("Example 2: 1D Oscillating function\n")
cat("================================================================\n\n")

for (ev in c("Grid", "Training")) {
  tag <- if (ev == "Grid") "grid" else "train"
  cat(sprintf("--- %s ---\n", ev))
  cat(sprintf("  %-30s  %12s\n", "Config", "Coverage"))
  cat(sprintf("  %-30s  %12s\n", "------", "--------"))
  for (nm in names(results)) {
    r <- results[[nm]]
    cat(sprintf("\n  %s\n", nm))
    for (method in c("bn", "bs", "lg")) {
      s <- r[[paste0(tag, "_", method)]]
      lbl <- switch(method, bn = "BKP", bs = "BKP-ESS", lg = "LGP")
      cat(sprintf("    %-13s  %5.1f%%(SD=%4.1f)\n",
                  lbl, 100 * s["mu"], 100 * s["sd"]))
    }
  }
  cat("\n")
}
cat("================================================================\n")
cat("Done.\n")