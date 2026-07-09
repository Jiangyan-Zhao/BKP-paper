# =============================================================
# Mourning Warbler: BKP vs TwinBKP vs LGP
# Species distribution model on 8 bioclimatic covariates
# Classification task, r0 = 0.1
# =============================================================

library(BKP)
library(gplite)
library(ggplot2)
library(dplyr)
library(maps)
library(pROC)
library(gridExtra)
library(terra)
library(RColorBrewer)
set.seed(123)

# ---- Load data ---------------------------------------------------
data <- read.csv("code/data/Mourning_Warbler.csv", stringsAsFactors = FALSE)
cat(sprintf("Data dimensions: %d rows x %d cols\n", nrow(data), ncol(data)))
cat("Split distribution:\n"); print(table(data$split))

train_idx <- which(data$split == "training")
test_idx  <- which(data$split == "testing")
cat(sprintf("Training set: %d | Test set: %d\n", length(train_idx), length(test_idx)))

bio_cols <- grep("^wc2", names(data), value = TRUE)
cat(sprintf("Covariates: %s\n", paste(bio_cols, collapse = ", ")))

X_all <- as.matrix(data[, bio_cols])
X_train <- X_all[train_idx, ]
y_train <- data$y[train_idx]; m_train <- rep(1L, length(y_train))
X_test  <- X_all[test_idx, ]
y_test  <- data$y[test_idx]
Xbounds <- apply(X_train, 2, range); Xbounds <- t(Xbounds)

# ---- Load climate rasters ----------------------------------------
cat("\n===== Loading climate rasters =====\n")
tif_dir <- "code/data/climate/wc2.1_10m"
tif_files <- file.path(tif_dir, paste0(bio_cols, ".tif"))
na_clim <- rast(tif_files)
e <- ext(-170.15, -46.70, 22.23, 71.90)
na_clim <- crop(na_clim, e)
names(na_clim) <- bio_cols
cat(sprintf("Cropped raster: %d rows x %d cols x %d layers\n",
            nrow(na_clim), ncol(na_clim), nlyr(na_clim)))

# ---- Colors and base map prep ------------------------------------
make_rast <- function(vals, template) { 
  r <- rast(template[[1]])
  values(r) <- vals
  return(r) 
}

purples <- colorRampPalette(brewer.pal(9, "Purples")[-(1:2)])
greens  <- colorRampPalette(brewer.pal(9, "Greens")[-(1:2)])
prescol <- purples(3)[3]; abscol <- purples(3)[2]

pts_train_pres <- data[data$y == 1 & data$split == "training", c("lon", "lat")]
pts_train_abs  <- data[data$y == 0 & data$split == "training", c("lon", "lat")]
pts_test       <- data[data$split == "testing", c("lon", "lat")]

world_map <- map_data("world")
na_map <- world_map %>%
  filter(long >= -170, long <= -50, lat >= 15, lat <= 75)

# land mask from climate raster
bg_rast <- na_clim[["wc2.1_10m_bio_5"]]
bg_vals <- values(bg_rast, na.rm = FALSE)
bg_vals[!is.na(bg_vals)] <- 1
values(bg_rast) <- bg_vals
bg_df <- as.data.frame(bg_rast, xy = TRUE)
colnames(bg_df)[3] <- "land"

# ---- Distribution data map ---------------------------------------
map_distribution <- ggplot() +
  geom_tile(data = bg_df, aes(x = x, y = y, fill = land)) +
  scale_fill_gradientn(colours = grey(0.8), na.value = "transparent", guide = "none") +
  geom_point(data = pts_test, aes(x = lon, y = lat),
             shape = 21, size = 0.6, stroke = 0.15,
             fill = "white", colour = grey(0.4)) +
  geom_point(data = pts_train_abs, aes(x = lon, y = lat),
             shape = 21, size = 0.6, stroke = 0.15,
             fill = abscol, colour = abscol) +
  geom_point(data = pts_train_pres, aes(x = lon, y = lat),
             shape = 21, size = 0.6, stroke = 0.15,
             fill = prescol, colour = prescol) +
  coord_quickmap(xlim = c(-170, -50), ylim = c(15, 75), expand = FALSE) +
  labs(title = "Distribution data - Mourning Warbler", x = NULL, y = NULL) +
  theme_void(12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("code/figure/mourning_warbler_map_distribution.pdf", 
       plot = map_distribution,
       width = 6.5, height = 3.8, dpi = 300)

# ---- Model fitting -----------------------------------------------
all_rocs <- list(); all_preds <- list(); times <- c()

cat("\n===== BKP =====\n")
t1 <- system.time({
  bkp_fit <- fit_BKP(X_train, y_train, m_train, Xbounds = Xbounds,
                     prior = "fixed", r0 = 0.1, loss = "log_loss",
                     kernel = "gaussian", isotropic = TRUE)
  bkp_pred_test <- predict(bkp_fit, Xnew = X_test)
})
times <- c(times, BKP = unname(t1["elapsed"]))
bkp_roc <- roc(y_test, bkp_pred_test$mean)
all_rocs$BKP <- bkp_roc; all_preds$BKP <- bkp_pred_test$mean


cat("\n===== TwinBKP =====\n")
t2 <- system.time({
  twin_fit <- fit_TwinBKP(X_train, y_train, m_train, Xbounds = Xbounds,
                          prior = "fixed", r0 = 0.1, loss = "log_loss",
                          global_kernel = "gaussian", local_kernel = "wendland",
                          isotropic = TRUE)
  twin_pred_test <- predict(twin_fit, Xnew = X_test)
})
times <- c(times, TwinBKP = unname(t2["elapsed"]))
twin_roc <- roc(y_test, twin_pred_test$mean)
all_rocs$TwinBKP <- twin_roc; all_preds$TwinBKP <- twin_pred_test$mean


cat("\n===== LGP =====\n")
t3 <- system.time({ 
  gp <- gp_init(cf = cf_sexp(), lik = lik_bernoulli())
  gp <- gp_optim(gp, X_train, y_train, verbose = FALSE) 
  gp_pred_test <- gp_pred(gp, X_test, transform = TRUE)
})
times <- c(times, LGP = unname(t3["elapsed"]))
gp_roc <- roc(y_test, gp_pred_test$mean)
all_rocs$LGP <- gp_roc; all_preds$LGP <- gp_pred_test$mean


# ---- Comparison ---------------------------------------------------
cat("\n=============================================================\n")
cat("  Comparison: Mourning Warbler (8 bio vars, r0=0.1)\n")
cat("=============================================================\n")
for (nm in names(all_rocs)) {
  cat(sprintf("  %-12s | AUC: %.3f | Brier: %.6f | Time: %.2f s\n",
              nm, auc(all_rocs[[nm]]), mean((all_preds[[nm]] - y_test)^2), times[nm]))
}

# ---- ROC curve ----------------------------------------------------
n_methods <- length(all_rocs)
cols <- c("#1c61b6", "#fdae61", "#d7191c")
nm_pad <- format(names(all_rocs), justify = "none",
                 width = max(nchar(names(all_rocs))))
legend_lines <- sprintf("%s  (AUC=%.3f)", nm_pad, sapply(all_rocs, auc))
pdf("code/figure/mourning_warbler_roc_comparison.pdf", 6, 6)
plot(all_rocs[[1]], col = cols[1], lwd = 2,
     main = "ROC Curve: Mourning Warbler")
for (i in seq_len(n_methods)[-1])
  lines(all_rocs[[i]], col = cols[i], lwd = 2)
legend("bottomright", legend = legend_lines,
       col = cols[1:n_methods], lwd = 2,
       bty = "o", cex = 0.85, bg = "white")
dev.off()

# ---- Grid prediction (climate raster) -----------------------------
cat("\n===== Grid prediction =====\n")

grid_vals <- values(na_clim, na.rm = FALSE)
colnames(grid_vals) <- bio_cols
grid_n <- nrow(grid_vals)
grid_valid <- complete.cases(grid_vals)
cat(sprintf("Grid points: %d (land %d / ocean %d)\n",
            grid_n, sum(grid_valid), sum(!grid_valid)))

raster_preds <- list()
for (method in names(all_preds)) {
  cat(sprintf("  %s grid prediction...\n", method))
  gm <- rep(NA_real_, grid_n); gv <- rep(NA_real_, grid_n)
  if (any(grid_valid)) {
    X_grid <- grid_vals[grid_valid, bio_cols, drop = FALSE]
    if (method == "LGP") {
      chunk_size <- 5000
      n_valid <- nrow(X_grid)
      valid_idx <- which(grid_valid)
      
      for (k in seq(1, n_valid, by = chunk_size)) {
        idx_end <- min(k + chunk_size - 1, n_valid)
        idx <- k:idx_end
        chunk <- X_grid[idx, , drop = FALSE]
        
        pred_chunk <- gp_pred(gp, chunk, transform = TRUE, var = TRUE)
        
        gm[valid_idx[idx]] <- as.vector(pred_chunk$mean)
        gv[valid_idx[idx]] <- as.vector(pred_chunk$var)
      }
    } else {
      fit_obj <- if (method == "BKP") bkp_fit else twin_fit
      pred_grid <- predict(fit_obj, Xnew = X_grid)
      gm[grid_valid] <- as.vector(pred_grid$mean)
      gv[grid_valid] <- as.vector(pred_grid$variance)
    }
  }
  raster_preds[[method]] <- list(mean = gm, variance = gv)
}

# ---- Per‑model maps (probability + variance) ----------------------
cat("\n===== Generating maps =====\n")
all_probs <- list()
all_vars  <- list()
for (method in names(raster_preds)) {
  cat(sprintf("  Generating %s maps...\n", method))
  
  rp <- raster_preds[[method]]
  r_prob <- make_rast(rp$mean, na_clim)
  df_prob <- as.data.frame(r_prob, xy = TRUE); colnames(df_prob)[3] <- "prob"
  p_prob <- ggplot() +
    geom_tile(data = df_prob, aes(x = x, y = y, fill = prob)) +
    scale_fill_gradientn(
      colours = purples(1000),
      name = "Probability of presence",
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = unit(0.2, "cm"), barheight = unit(4, "cm"),
        frame.colour = NA, ticks.colour = NA, ticks.linewidth = 0,
        title.position = "right"
      )
    ) +
    coord_quickmap(xlim = c(-170, -50), ylim = c(15, 75), expand = FALSE) +
    labs(title = paste(method, "- Predicted distribution"), x = NULL, y = NULL) +
    theme_void(12) + theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.key.width  = unit(0.2, "cm"),
      legend.key.height = unit(4, "cm"),
      legend.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
    )
  # ggsave(sprintf("mourning_warbler_map_%s_prob.pdf", tolower(gsub(" ", "_", method))),
         # p_prob, width = 6.5, height = 3.8, dpi = 300)
  
  r_var <- make_rast(rp$variance, na_clim)
  df_var <- as.data.frame(r_var, xy = TRUE); colnames(df_var)[3] <- "var"
  vmax <- max(rp$variance, na.rm = TRUE)
  vlim <- ceiling(vmax * 10) / 10
  p_var <- ggplot() +
    geom_tile(data = df_var, aes(x = x, y = y, fill = var)) +
    scale_fill_gradientn(
      colours = greens(1000),
      name = paste("Variance of", method, "posterior"),
      limits = c(0, vlim), breaks = seq(0, vlim, 0.1),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = unit(0.2, "cm"), barheight = unit(4, "cm"),
        frame.colour = NA, ticks.colour = NA, ticks.linewidth = 0,
        title.position = "right"
      )
    ) +
    coord_quickmap(xlim = c(-170, -50), ylim = c(15, 75), expand = FALSE) +
    labs(title = paste(method, "- Prediction uncertainty"), x = NULL, y = NULL) +
    theme_void(12) + theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.key.width  = unit(0.2, "cm"),
      legend.key.height = unit(4, "cm"),
      legend.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
    )
  # ggsave(sprintf("mourning_warbler_map_%s_var.pdf", tolower(gsub(" ", "_", method))),
  #        p_var, width = 6.5, height = 3.8, dpi = 300)
  
  all_probs[[method]] <- p_prob
  all_vars[[method]]  <- p_var
}

# ---- Combined 3×2 map ----------------------------------------------
cat("\n===== Generating combined map (3x2) =====\n")
p_combined <- grid.arrange(
  all_probs$BKP,     all_vars$BKP,
  all_probs$TwinBKP, all_vars$TwinBKP,
  all_probs$LGP,     all_vars$LGP,
  ncol = 2, nrow = 3
)
ggsave("code/figure/mourning_warbler_map_combined.pdf", plot = p_combined,
       width = 13.0, height = 11.4, dpi = 300)

cat("\n===== Analysis complete =====\n")