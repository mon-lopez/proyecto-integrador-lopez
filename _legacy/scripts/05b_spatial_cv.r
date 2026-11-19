# ==============================================================================
# Script: 05b_spatial_cv.R
# Proyecto: Análisis Espacial de Emisiones de CO2 en Colombia
# Sección: 3. Machine Learning - Validación Cruzada Espacial
# Objetivo: Implementar Spatial Cross-Validation con bloques geográficos para 
#           evaluar performance de modelos predictivos corrigiendo sesgos espaciales.
# Referencias: Meyer et al. (2023), blockCV 3.0.
# ==============================================================================

# 0. CONFIGURACIÓN Y LIBRERÍAS ------------------------------------------------
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(ranger)        # Random Forests
  library(blockCV)       # Spatial Cross-Validation
  library(viridis)
  library(patchwork)
})

set.seed(2026)
theme_set(theme_minimal())

dir.create("results", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# ==============================================================================
# 1. CARGA Y PREPARACIÓN DE DATOS ---------------------------------------------
# ==============================================================================
cat("=== 1. Cargando datos ===\n")

co2_raw <- read.csv("climate_trace_COL_2023.csv", stringsAsFactors = FALSE)

co2_clean <- co2_raw %>%
  filter(!is.na(longitude), !is.na(latitude),
         longitude > -82 & longitude < -66,
         latitude > -5 & latitude < 13,
         !is.na(co2), is.finite(co2), co2 >= 0) %>%
  distinct(longitude, latitude, sector, sub_sector, .keep_all = TRUE)

# Outcome y covariables
co2_clean <- co2_clean %>%
  mutate(
    log_co2 = log1p(co2),
    log_actividad = log1p(activity_data),
    sector_idx = as.integer(as.factor(sector)),
    sub_sector_idx = as.integer(as.factor(sub_sector)),
    lat = latitude,
    lon = longitude,
    lat2 = lat^2,
    lon2 = lon^2,
    lat_lon = lat * lon
  )

co2_sf <- st_as_sf(co2_clean, coords = c("longitude", "latitude"), crs = 4326)

cat(sprintf("  • %d observaciones cargadas\n", nrow(co2_sf)))

# ==============================================================================
# 2. CREACIÓN DE BLOQUES ESPACIALES (blockCV) ---------------------------------
# ==============================================================================
cat("\n=== 2. Creando bloques espaciales para CV ===\n")

# --- 2.1 Definir tamaño de bloque ---
# Tamaño recomendado: al menos 2x el rango de autocorrelación espacial
# Para Colombia, usamos bloques de ~100 km (aprox 1 grado)
block_size <- 100000  # metros (~100 km)

# Transformar a CRS métrico para blockCV
co2_metric <- st_transform(co2_sf, 9377)  # MAGNA-SIRGAS

cat(sprintf("  • Tamaño de bloque: %d metros (~%.0f km)\n", 
            block_size, block_size / 1000))

# --- 2.2 Generar bloques espaciales ---
cat("  • Generando bloques espaciales...\n")

sp_blocks <- spatialBlock(
  speciesData = co2_metric,
  theRange = block_size,
  k = 5,              # 5 folds
  selection = "random",
  iteration = 100,
  report = TRUE,
  plot = FALSE
)

cat(sprintf("  • %d bloques creados\n", length(unique(sp_blocks$foldID))))

# --- 2.3 Visualizar bloques ---
co2_metric$fold_spatial <- as.factor(sp_blocks$foldID)

p_blocks <- ggplot() +
  geom_sf(data = co2_metric, aes(fill = fold_spatial), 
          color = NA, size = 0.5) +
  scale_fill_viridis_d(name = "Fold") +
  labs(
    title = "Spatial Cross-Validation Blocks",
    subtitle = sprintf("5 folds | Block size = %.0f km", block_size / 1000),
    caption = "Cada color representa un fold en la validación cruzada"
  ) +
  theme(legend.position = "bottom")

ggsave("results/01_spatial_blocks.png", p_blocks, width = 10, height = 10, dpi = 300)
cat("  ✓ Mapa de bloques guardado\n")

# ==============================================================================
# 3. SPATIAL CV CON RANDOM FOREST ---------------------------------------------
# ==============================================================================
cat("\n=== 3. Spatial CV con Random Forest ===\n")

# Preparar datos
X_cols <- c("sector_idx", "sub_sector_idx", "log_actividad", 
            "lat", "lon", "lat2", "lon2", "lat_lon")

# Función para entrenar y evaluar RF en cada fold
spatial_cv_rf <- function(fold_id, train_idx, test_idx, data, X_cols) {
  train_data <- data[train_idx, ]
  test_data <- data[test_idx, ]
  
  # Entrenar RF
  rf <- ranger(
    formula = log_co2 ~ .,
    data = train_data[, c("log_co2", X_cols)],
    num.trees = 200,
    mtry = floor(sqrt(length(X_cols))),
    min.node.size = 5,
    num.threads = 1
  )
  
  # Predecir
  pred <- predict(rf, data = test_data[, X_cols])$predictions
  
  # Métricas
  residuals <- test_data$log_co2 - pred
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  r2 <- 1 - sum(residuals^2) / sum((test_data$log_co2 - mean(test_data$log_co2))^2)
  
  # Coverage de IC 95% (aproximación con cuartiles de predicciones)
  # Usamos predicciones de todos los árboles para estimar incertidumbre
  pred_all_trees <- predict(rf, data = test_data[, X_cols], 
                            predict.all = TRUE)$predictions
  pred_lower <- apply(pred_all_trees, 1, quantile, 0.025)
  pred_upper <- apply(pred_all_trees, 1, quantile, 0.975)
  
  coverage <- mean(test_data$log_co2 >= pred_lower & test_data$log_co2 <= pred_upper)
  
  data.frame(
    fold = fold_id,
    rmse = rmse,
    mae = mae,
    r2 = r2,
    coverage = coverage,
    n_train = length(train_idx),
    n_test = length(test_idx)
  )
}

# Ejecutar Spatial CV
cat("  • Ejecutando 5-fold Spatial CV...\n")
spatial_cv_results <- list()

for (fold in 1:5) {
  train_idx <- which(sp_blocks$foldID != fold)
  test_idx <- which(sp_blocks$foldID == fold)
  
  cat(sprintf("    Fold %d: %d train, %d test\n", 
              fold, length(train_idx), length(test_idx)))
  
  result <- spatial_cv_rf(fold, train_idx, test_idx, co2_metric, X_cols)
  spatial_cv_results[[fold]] <- result
}

spatial_cv_df <- do.call(rbind, spatial_cv_results)

cat("\n  • Resultados Spatial CV (Random Forest):\n")
print(spatial_cv_df)

# Métricas promedio
spatial_cv_summary <- spatial_cv_df %>%
  summarise(
    rmse_mean = mean(rmse),
    rmse_sd = sd(rmse),
    mae_mean = mean(mae),
    r2_mean = mean(r2),
    r2_sd = sd(r2),
    coverage_mean = mean(coverage)
  )

cat(sprintf("\n  • Resumen Spatial CV:\n"))
cat(sprintf("    RMSE: %.4f ± %.4f\n", spatial_cv_summary$rmse_mean, spatial_cv_summary$rmse_sd))
cat(sprintf("    R²: %.4f ± %.4f\n", spatial_cv_summary$r2_mean, spatial_cv_summary$r2_sd))
cat(sprintf("    Coverage IC 95%%: %.3f\n", spatial_cv_summary$coverage_mean))

# ==============================================================================
# 4. CV ALEATORIO (COMPARACIÓN) -----------------------------------------------
# ==============================================================================
cat("\n=== 4. CV Aleatorio (comparación) ===\n")

# Asignar folds aleatorios
set.seed(2026)
co2_metric$fold_random <- as.factor(sample(1:5, nrow(co2_metric), replace = TRUE))

# Ejecutar Random CV
random_cv_results <- list()

for (fold in 1:5) {
  train_idx <- which(co2_metric$fold_random != fold)
  test_idx <- which(co2_metric$fold_random == fold)
  
  result <- spatial_cv_rf(fold, train_idx, test_idx, co2_metric, X_cols)
  random_cv_results[[fold]] <- result
}

random_cv_df <- do.call(rbind, random_cv_results)

cat("\n  • Resultados Random CV (Random Forest):\n")
print(random_cv_df)

random_cv_summary <- random_cv_df %>%
  summarise(
    rmse_mean = mean(rmse),
    rmse_sd = sd(rmse),
    mae_mean = mean(mae),
    r2_mean = mean(r2),
    r2_sd = sd(r2),
    coverage_mean = mean(coverage)
  )

cat(sprintf("\n  • Resumen Random CV:\n"))
cat(sprintf("    RMSE: %.4f ± %.4f\n", random_cv_summary$rmse_mean, random_cv_summary$rmse_sd))
cat(sprintf("    R²: %.4f ± %.4f\n", random_cv_summary$r2_mean, random_cv_summary$r2_sd))
cat(sprintf("    Coverage IC 95%%: %.3f\n", random_cv_summary$coverage_mean))

# ==============================================================================
# 5. COMPARACIÓN SPATIAL CV vs RANDOM CV --------------------------------------
# ==============================================================================
cat("\n=== 5. Comparación Spatial CV vs Random CV ===\n")

comparison_cv <- data.frame(
  Metodo = c("Spatial CV", "Random CV"),
  RMSE_mean = c(spatial_cv_summary$rmse_mean, random_cv_summary$rmse_mean),
  RMSE_sd = c(spatial_cv_summary$rmse_sd, random_cv_summary$rmse_sd),
  R2_mean = c(spatial_cv_summary$r2_mean, random_cv_summary$r2_mean),
  R2_sd = c(spatial_cv_summary$r2_sd, random_cv_summary$r2_sd),
  Coverage = c(spatial_cv_summary$coverage_mean, random_cv_summary$coverage_mean)
)

cat("\n  • Tabla comparativa:\n")
print(comparison_cv)

write.csv(comparison_cv, "results/02_cv_comparison.csv", row.names = FALSE)

# Visualización comparativa
p_comparison <- ggplot() +
  geom_point(data = comparison_cv, 
             aes(x = Metodo, y = R2_mean, color = Metodo), 
             size = 4) +
  geom_errorbar(data = comparison_cv, 
                aes(x = Metodo, ymin = R2_mean - R2_sd, ymax = R2_mean + R2_sd, 
                    color = Metodo), 
                width = 0.2) +
  labs(
    title = "Spatial CV vs Random CV",
    subtitle = "R² promedio por fold (Random Forest)",
    x = "Método de CV", y = "R²",
    caption = "Spatial CV típicamente muestra menor R² (más realista)"
  ) +
  scale_color_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("results/03_cv_comparison.png", p_comparison, width = 8, height = 6, dpi = 300)
cat("  ✓ Gráfico comparativo guardado\n")

# ==============================================================================
# 6. ANÁLISIS DE PERFORMANCE POR REGIÓN ---------------------------------------
# ==============================================================================
cat("\n=== 6. Performance por región geográfica ===\n")

# Dividir en regiones (Norte/Sur/Este/Oeste)
co2_metric <- co2_metric %>%
  mutate(
    region_ns = ifelse(st_coordinates(.)[,2] > median(st_coordinates(.)[,2]), 
                       "Norte", "Sur"),
    region_ew = ifelse(st_coordinates(.)[,1] > median(st_coordinates(.)[,1]), 
                       "Oriente", "Occidente")
  )

# Evaluar performance del modelo final en cada región
# Entrenar modelo final con todos los datos
rf_final <- ranger(
  formula = log_co2 ~ .,
  data = co2_metric[, c("log_co2", X_cols)],
  num.trees = 500,
  mtry = floor(sqrt(length(X_cols))),
  min.node.size = 5,
  importance = "impurity",
  num.threads = 1
)

pred_final <- predict(rf_final, data = co2_metric[, X_cols])$predictions
co2_metric$residual <- co2_metric$log_co2 - pred_final

performance_region <- co2_metric %>%
  group_by(region_ns, region_ew) %>%
  summarise(
    n = n(),
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    r2 = 1 - sum(residual^2) / sum((log_co2 - mean(log_co2))^2),
    .groups = "drop"
  )

cat("\n  • Performance por región:\n")
print(performance_region)

write.csv(performance_region, "results/04_performance_by_region.csv", row.names = FALSE)

# Mapa de residuos
p_residuos <- ggplot() +
  geom_sf(data = co2_metric, aes(fill = residual), 
          color = NA, size = 0.5) +
  scale_fill_gradient2(low = "#1a9850", mid = "white", high = "#d73027",
                       midpoint = 0, name = "Residuo",
                       limits = quantile(co2_metric$residual, c(0.02, 0.98))) +
  labs(
    title = "Residuos del Modelo Final (Random Forest)",
    subtitle = "Spatial CV R² = 0.XX",
    caption = "Rojo = sobreestimación, Verde = subestimación"
  ) +
  theme(legend.position = "bottom")

ggsave("results/05_residuals_map.png", p_residuos, width = 10, height = 10, dpi = 300)
cat("  ✓ Mapa de residuos guardado\n")

# ==============================================================================
# 7. REPORTE DE PERFORMANCE FINAL ---------------------------------------------
# ==============================================================================
cat("\n=== 7. Reporte de performance final ===\n")

reporte_final <- data.frame(
  Metrica = c(
    "N observaciones",
    "N variables predictoras",
    "N folds (Spatial CV)",
    "Tamaño de bloque (km)",
    "RMSE (Spatial CV)",
    "SE del RMSE",
    "R² (Spatial CV)",
    "SE del R²",
    "Coverage IC 95%",
    "RMSE (Random CV)",
    "R² (Random CV)",
    "Diferencia R² (Spatial - Random)"
  ),
  Valor = c(
    nrow(co2_metric),
    length(X_cols),
    5,
    block_size / 1000,
    round(spatial_cv_summary$rmse_mean, 4),
    round(spatial_cv_summary$rmse_sd, 4),
    round(spatial_cv_summary$r2_mean, 4),
    round(spatial_cv_summary$r2_sd, 4),
    round(spatial_cv_summary$coverage_mean, 3),
    round(random_cv_summary$rmse_mean, 4),
    round(random_cv_summary$r2_mean, 4),
    round(spatial_cv_summary$r2_mean - random_cv_summary$r2_mean, 4)
  )
)

write.csv(reporte_final, "results/06_performance_report.csv", row.names = FALSE)

cat("\n  Archivos generados:\n")
cat("  ✓ results/01_spatial_blocks.png          - Mapa de bloques espaciales\n")
cat("  ✓ results/02_cv_comparison.csv           - Comparación Spatial vs Random CV\n")
cat("  ✓ results/03_cv_comparison.png           - Gráfico comparativo\n")
cat("  ✓ results/04_performance_by_region.csv   - Performance por región\n")
cat("  ✓ results/05_residuals_map.png           - Mapa de residuos\n")
cat("  ✓ results/06_performance_report.csv      - Reporte final\n")

# Guardar modelo final
saveRDS(rf_final, "data/rf_final_model.rds")
cat("  ✓ Modelo final guardado\n")

cat("\n=== Spatial Cross-Validation completado ===\n")

cat("\n⚠ NOTA METODOLÓGICA PARA EL PAPER:\n")
cat("  Spatial CV con blockCV corrige el sesgo de optimismo en la estimación\n")
cat("  de performance causado por autocorrelación espacial. La diferencia entre\n")
cat("  R² de Spatial CV y Random CV cuantifica este sesgo.\n")