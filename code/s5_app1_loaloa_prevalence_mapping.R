# ==============================================================
# Real Data Example: Loaloa (BKP vs LGP)
# ==============================================================
library(BKP)
library(gplite)
library(ggplot2)
library(gridExtra)
library(pROC)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(RiskMap)
library(viridis)

options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)

set.seed(123)

# ---- Load data -------------------------------------------------
data("loaloa", package = "RiskMap") 

# Extract input variables (X), response variable (y), and trial counts (m)
X <- as.matrix(loaloa[, c("LONGITUDE", "LATITUDE")])
rownames(X) <- NULL
y <- loaloa$NO_INF
m <- loaloa$NO_EXAM

# Randomly split into training (70%) and testing (30%) sets
train_idx <- sample(1:nrow(loaloa), 0.7 * nrow(loaloa))
X_train <- X[train_idx, ]
y_train <- y[train_idx]
m_train <- m[train_idx]
X_test  <- X[-train_idx, ]
y_test  <- y[-train_idx]
m_test  <- m[-train_idx]

p        <- y / m # Infection rate
Xbounds  <- matrix(c(7.8, 15.3, 3.1, 7.0), ncol = 2, byrow = TRUE) 
pi_test  <- y_test / m_test

# ---- Colors and base map ---------------------------------------
africa <- ne_countries(continent = "africa", scale = "medium", returnclass = "sf")

# ---- Distribution map ---------------------------------------------
df_map <- data.frame(
  lon = X[,1], lat = X[,2], p = p, n = m,
  set = ifelse(1:nrow(loaloa) %in% train_idx, "Train", "Test")
)

pmap <- ggplot() +
  geom_sf(data = africa, fill = "gray95", color = "gray60") +
  geom_point(data = df_map,
             aes(x = lon, y = lat, color = p, size = n, shape = set),
             alpha = 0.8) +
  scale_color_viridis_c(name = "Proportion y/m") +
  scale_size_continuous(name = "Trial count m") +
  scale_shape_manual(name = "Dataset", values = c("Train" = 16, "Test" = 17)) +
  coord_sf(xlim = Xbounds[1, ], ylim = Xbounds[2, ], expand = FALSE) +
  theme_minimal() +
  labs(title = "Loaloa infection proportion (y/m): Train vs Test",
       x = "Longitude", y = "Latitude")
ggsave("code/figure/Loaloa_map.pdf", plot = pmap, width = 10, height = 5, dpi = 300)

# ---- Model fitting ------------------------------------------------
all_preds <- list()

# BKP 
bkp_fit <- fit_BKP(X_train, y_train, m_train, Xbounds = Xbounds)
all_preds$BKP <- predict(bkp_fit, Xnew = X_test)$mean
summary(bkp_fit)

# LGP 
gp <- gp_init(cf = cf_sexp(), lik = lik_binomial())
gp <- gp_optim(gp, X_train, y_train, trials = m_train, verbose = FALSE)
all_preds$LGP <- gp_pred(gp, X_test, transform = TRUE)$mean

# ---- Comparison ----------------------------------------------------
cat("\n=============================================================\n")
cat("  Comparison: Loaloa (BKP vs LGP)\n")
cat("=============================================================\n")
for (nm in names(all_preds)) {
  cat(sprintf("  %-8s | Brier: %.6f \n",
              nm, mean((all_preds[[nm]] - pi_test)^2)))
}

# ---- Land‑only grid prediction ------------------------------------
x1_seq <- seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 300)
x2_seq <- seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 260)
grid   <- expand.grid(lon = x1_seq, lat = x2_seq)

grid_sf   <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)
on_land   <- lengths(st_intersects(grid_sf, africa)) > 0
grid_land <- grid[on_land, , drop = FALSE]
X_grid    <- as.matrix(grid_land)

# BKP grid prediction
pred_bkp <- predict(bkp_fit, Xnew = X_grid)
grid_land$BKP_mean <- as.vector(pred_bkp$mean)
grid_land$BKP_var  <- as.vector(pred_bkp$variance)

# LGP grid prediction (batched)
chunk_size <- 5000
gm <- rep(NA_real_, nrow(X_grid))
gv <- rep(NA_real_, nrow(X_grid))
for (k in seq(1, nrow(X_grid), by = chunk_size)) {
  idx_end <- min(k + chunk_size - 1, nrow(X_grid))
  chunk   <- X_grid[k:idx_end, , drop = FALSE]
  pc      <- gp_pred(gp, chunk, transform = TRUE, var = TRUE)
  gm[k:idx_end] <- as.vector(pc$mean)
  gv[k:idx_end] <- as.vector(pc$var)
}
grid_land$LGP_mean <- gm
grid_land$LGP_var  <- gv

raster_preds <- list(
  BKP = list(mean = grid_land$BKP_mean, variance = grid_land$BKP_var),
  LGP = list(mean = grid_land$LGP_mean, variance = grid_land$LGP_var)
)

# ---- Unified color scale limits across methods ---------------------
all_pmax <- sapply(raster_preds, function(x) max(x$mean, na.rm = TRUE))
plim <- ceiling(max(all_pmax) * 10) / 10

all_vmax <- sapply(raster_preds, function(x) max(x$variance, na.rm = TRUE))
vlim <- ceiling(max(all_vmax) * 10) / 10

# ---- Per‑model maps (probability + variance) ----------------------
sub_labels <- c("(a)", "(b)", "(c)", "(d)")
sub_idx <- 1

all_probs <- list()
all_vars  <- list()
for (method in names(raster_preds)) {
  rp <- raster_preds[[method]]
  
  df_method <- data.frame(
    lon      = grid_land$lon,
    lat      = grid_land$lat,
    mean     = rp$mean,
    variance = rp$variance
  )
  
  # Probability map
  p_prob <- ggplot() +
    geom_sf(data = africa, fill = grey(0.8), colour = "gray60") +
    geom_tile(data = df_method, aes(x = lon, y = lat, fill = mean)) +
    scale_fill_viridis_c(
      option  = "C",
      name    = "Probability of prevalence",
      limits  = c(0, plim), breaks = seq(0, plim, 0.2),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = unit(0.2, "cm"), barheight = unit(4, "cm"),
        frame.colour = NA, ticks.colour = NA, ticks.linewidth = 0,
        title.position = "right"
      )
    ) +
    coord_sf(xlim = Xbounds[1, ], ylim = Xbounds[2, ], expand = FALSE) +
    labs(title = paste(sub_labels[sub_idx], method, "- Predicted distribution"),
         x = "Longitude", y = "Latitude") +
    theme_minimal(12) + theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.key.width  = unit(0.2, "cm"),
      legend.key.height = unit(4, "cm"),
      legend.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
    )
  sub_idx <- sub_idx + 1
  # ggsave(sprintf("Loaloa_%s_prob.pdf", tolower(method)),
  #        p_prob, width = 6.5, height = 3.8, dpi = 300)
  
  # Variance map
  p_var <- ggplot() +
    geom_sf(data = africa, fill = grey(0.8), colour = "gray60") +
    geom_tile(data = df_method, aes(x = lon, y = lat, fill = variance)) +
    scale_fill_viridis_c(
      option  = "C",
      name    = paste("Variance of", method, "posterior"),
      limits  = c(0, vlim), breaks = seq(0, vlim, 0.1),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = unit(0.2, "cm"), barheight = unit(4, "cm"),
        frame.colour = NA, ticks.colour = NA, ticks.linewidth = 0,
        title.position = "right"
      )
    ) +
    coord_sf(xlim = Xbounds[1, ], ylim = Xbounds[2, ], expand = FALSE) +
    labs(title = paste(sub_labels[sub_idx], method, "- Prediction uncertainty"),
         x = "Longitude", y = "Latitude") +
    theme_minimal(12) + theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.key.width  = unit(0.2, "cm"),
      legend.key.height = unit(4, "cm"),
      legend.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10)
    )
  sub_idx <- sub_idx + 1
  # ggsave(sprintf("Loaloa_%s_var.pdf", tolower(method)),
  #        p_var, width = 6.5, height = 3.8, dpi = 300)
  
  all_probs[[method]] <- p_prob
  all_vars[[method]]  <- p_var
}

# ---- Combined 2×2 map ----------------------------------------------
cat("\n===== Generating combined map (2x2) =====\n")
p_combined <- grid.arrange(
  all_probs$BKP, all_vars$BKP,
  all_probs$LGP, all_vars$LGP,
  ncol = 2, nrow = 2
)
ggsave("code/figure/Loaloa_combined.pdf", plot = p_combined,
       width = 13, height = 9, dpi = 300)

cat("\n===== Analysis complete =====\n")
