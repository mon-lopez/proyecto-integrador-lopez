# ==============================================================================
# Script: 05_causal_ml.R
# Proyecto: Análisis Espacial de Emisiones de CO2 en Colombia
# Sección: 3. Machine Learning - Inferencia Causal y Predictiva
# Objetivo: Utilizar modelos predictivos y de inferencia causal para estimar 
#           efectos heterogéneos de tratamiento y corregir sesgos espaciales.
# Referencias: Athey, Tibshirani & Wager (2019); Meyer et al. (2023).
# ==============================================================================

# 0. CONFIGURACIÓN Y LIBRERÍAS ------------------------------------------------
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grf)           # Causal Forests
  library(ranger)        # Random Forests (más rápido que randomForest)
  library(missForest)    # Imputación con RF
  library(shapviz)       # SHAP values
  library(viridis)
  library(patchwork)
  library(broom)
})

set.seed(2026)
theme_set(theme_minimal())

# Crear carpetas
dir.create("results", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# ==============================================================================
# 1. CARGA Y CONSTRUCCIÓN DE DATOS --------------------------------------------
# ==============================================================================
cat("=== 1. Cargando y preparando datos ===\n")

co2_raw <- read.csv("climate_trace_COL_2023.csv", stringsAsFactors = FALSE)

# Filtrar Colombia y limpiar
co2_clean <- co2_raw %>%
  filter(!is.na(longitude), !is.na(latitude),
         longitude > -82 & longitude < -66,
         latitude > -5 & latitude < 13,
         !is.na(co2), is.finite(co2),
         co2 >= 0)  # Mantener solo emisiones no negativas para ML

# Eliminar duplicados exactos de coordenadas
co2_clean <- co2_clean %>%
  distinct(longitude, latitude, sector, sub_sector, .keep_all = TRUE)

N <- nrow(co2_clean)
cat(sprintf("  • %d observaciones únicas cargadas\n", N))

# --- 1.1 Definir TRATAMIENTO (binario: fósil vs renovable/secuestro) ---
# Basado en sector y sub_sector de Climate Trace
co2_clean <- co2_clean %>%
  mutate(
    # Tratamiento: 1 = fuente renovable/secuestro, 0 = fósil/emisora
    tratamiento = case_when(
      sector %in% c("manufacturing", "buildings", "transportation", 
                    "energy", "industry") ~ 0,  # Fósil
      sector == "agriculture" ~ 0,              # Agricultura = emisora
      sector == "waste" ~ 0,                    # Residuos = emisora
      sector == "forestry-and-land-use" & 
        grepl("removals|net-forest|sequestration", sub_sector, ignore.case = TRUE) ~ 1,
      sector == "forestry-and-land-use" & 
        grepl("fire|degradation|clearing|shrubgrass", sub_sector, ignore.case = TRUE) ~ 0,
      sector == "forestry-and-land-use" ~ 1,    # Forestal = secuestro (proxy)
      TRUE ~ NA_real_
    ),
    tratamiento = as.integer(tratamiento)
  )

# Filtrar solo observaciones con tratamiento definido
co2_ml <- co2_clean %>% filter(!is.na(tratamiento))

cat(sprintf("  • Tratamiento definido: %d obs\n", nrow(co2_ml)))
cat(sprintf("    - T=0 (fósil/emisora): %d (%.1f%%)\n", 
            sum(co2_ml$tratamiento == 0), 
            100 * mean(co2_ml$tratamiento == 0)))
cat(sprintf("    - T=1 (renovable/secuestro): %d (%.1f%%)\n", 
            sum(co2_ml$tratamiento == 1), 
            100 * mean(co2_ml$tratamiento == 1)))

# --- 1.2 OUTCOME: log de emisiones ---
co2_ml <- co2_ml %>%
  mutate(
    log_co2 = log1p(co2),
    log_actividad = log1p(activity_data)
  )

# --- 1.3 COVARIABLES OBSERVABLES (sin simulación) ---
# Todas derivadas directamente de los datos de Climate Trace
co2_ml <- co2_ml %>%
  mutate(
    # Variables categóricas codificadas
    sector_idx = as.integer(as.factor(sector)),
    sub_sector_idx = as.integer(as.factor(sub_sector)),
    
    # Variables continuas observables
    actividad = activity_data,
    log_actividad = log1p(activity_data),
    
    # Coordenadas (importantes para capturar efectos espaciales)
    lat = latitude,
    lon = longitude,
    
    # Términos no lineales de coordenadas (para capturar gradientes espaciales)
    lat2 = lat^2,
    lon2 = lon^2,
    lat_lon = lat * lon
  )

# --- 1.4 DENSIDAD LOCAL DE FUENTES (variable espacial derivada) ---
# Calcular número de fuentes en radio de 0.5° (~55 km)
cat("  • Calculando densidad local de fuentes...\n")

co2_sf <- st_as_sf(co2_ml, coords = c("longitude", "latitude"), crs = 4326)

# Usar buffer y conteo (alternativa rápida a knn)
# Para eficiencia, usamos una muestra si N > 5000
if (nrow(co2_sf) > 5000) {
  set.seed(2026)
  idx_sample <- sample(1:nrow(co2_sf), 5000)
  co2_sample <- co2_sf[idx_sample, ]
} else {
  co2_sample <- co2_sf
}

# Contar vecinos en radio de 0.5 grados
buffer_dist <- 0.5  # grados (~55 km en el ecuador)
buffers <- st_buffer(co2_sample, dist = buffer_dist)

# Contar cuántos puntos caen en cada buffer
densidad <- sapply(1:nrow(co2_sample), function(i) {
  length(st_intersects(buffers[i, ], co2_sample, sparse = FALSE)) - 1
})

co2_sample$densidad_local <- densidad

# Unir de vuelta al dataset completo (usando coordenadas como clave)
co2_ml$densidad_local <- NA_real_
for (i in 1:nrow(co2_sample)) {
  coords_i <- st_coordinates(co2_sample[i, ])
  match_idx <- which(
    abs(co2_ml$longitude - coords_i[1]) < 1e-6 & 
    abs(co2_ml$latitude - coords_i[2]) < 1e-6
  )
  if (length(match_idx) > 0) {
    co2_ml$densidad_local[match_idx] <- co2_sample$densidad_local[i]
  }
}

# Imputar NAs de densidad con mediana (puntos no muestreados)
co2_ml$densidad_local[is.na(co2_ml$densidad_local)] <- 
  median(co2_ml$densidad_local, na.rm = TRUE)

cat(sprintf("  • Densidad local: media = %.2f, mediana = %.2f\n", 
            mean(co2_ml$densidad_local), median(co2_ml$densidad_local)))

# --- 1.5 Matriz de covariables para ML ---
covariables_cols <- c("sector_idx", "sub_sector_idx", "log_actividad", 
                      "lat", "lon", "lat2", "lon2", "lat_lon", "densidad_local")

X <- as.matrix(co2_ml[, covariables_cols])
Y <- co2_ml$log_co2
W <- co2_ml$tratamiento

cat(sprintf("  • Matriz X: %d x %d\n", nrow(X), ncol(X)))
cat(sprintf("  • Outcome Y: rango [%.2f, %.2f]\n", min(Y), max(Y)))

# ==============================================================================
# 2. IMPUTACIÓN DE DATOS FALTANTES CON missForest -----------------------------
# ==============================================================================
cat("\n=== 2. Imputación con missForest (preservando estructura espacial) ===\n")

# Verificar NAs
n_na <- sum(is.na(X))
cat(sprintf("  • Valores NA en X: %d (%.2f%%)\n", n_na, 100 * n_na / length(X)))

if (n_na > 0) {
  cat("  • Ejecutando missForest...\n")
  df_impute <- co2_ml[, covariables_cols]
  
  impute_result <- missForest(df_impute, maxiter = 5, ntree = 100)
  X <- as.matrix(impute_result$ximp)
  
  cat(sprintf("  • Imputación completada. OOB error: %.4f\n", 
              impute_result$OOBerror))
} else {
  cat("  • Sin valores NA - saltando imputación\n")
}

# ==============================================================================
# 3. CAUSAL FORESTS PARA CATE -------------------------------------------------
# ==============================================================================
cat("\n=== 3. Causal Forest (grf) para estimar CATE ===\n")

# --- 3.1 Entrenar Causal Forest ---
# Usar submuestreo si N > 10000 para eficiencia
if (nrow(X) > 10000) {
  set.seed(2026)
  idx_cf <- sample(1:nrow(X), 10000)
  X_cf <- X[idx_cf, ]
  Y_cf <- Y[idx_cf]
  W_cf <- W[idx_cf]
  cat(sprintf("  • Submuestreo: %d → %d observaciones\n", nrow(X), length(idx_cf)))
} else {
  X_cf <- X
  Y_cf <- Y
  W_cf <- W
}

cat("  • Entrenando Causal Forest (500 árboles)...\n")
cf_model <- causal_forest(
  X = X_cf,
  Y = Y_cf,
  W = W_cf,
  num.trees = 500,
  sample.fraction = 0.5,
  mtry = min(ncol(X_cf), floor(sqrt(ncol(X_cf))) + 20),
  min.node.size = 5,
  honesty = TRUE,           # Honest estimation (Athey & Wager 2019)
  honesty.fraction = 0.5,
  num.threads = 1
)

cat("  ✓ Causal Forest entrenado\n")

# --- 3.2 Predicción de CATE ---
cat("  • Prediciendo CATE...\n")
cf_pred <- predict(cf_model)

# Extraer estimaciones y errores estándar
cate_estimates <- cf_pred$predictions
cate_se <- cf_pred$se

# Agregar al dataset
co2_ml$cate <- NA_real_
co2_ml$cate_se <- NA_real_

if (nrow(X) > 10000) {
  co2_ml$cate[idx_cf] <- cate_estimates
  co2_ml$cate_se[idx_cf] <- cate_se
  # Imputar CATE para observaciones no muestreadas con mediana
  co2_ml$cate[is.na(co2_ml$cate)] <- median(cate_estimates, na.rm = TRUE)
  co2_ml$cate_se[is.na(co2_ml$cate_se)] <- median(cate_se, na.rm = TRUE)
} else {
  co2_ml$cate <- cate_estimates
  co2_ml$cate_se <- cate_se
}

# --- 3.3 ATE promedio y test de significancia ---
ate_mean <- mean(cate_estimates)
ate_se <- sd(cate_estimates) / sqrt(length(cate_estimates))
ate_ci <- c(ate_mean - 1.96 * ate_se, ate_mean + 1.96 * ate_se)

cat(sprintf("\n  • ATE promedio: %.4f (SE = %.4f)\n", ate_mean, ate_se))
cat(sprintf("  • IC 95%%: [%.4f, %.4f]\n", ate_ci[1], ate_ci[2]))
cat(sprintf("  • p-valor: %.4f\n", 2 * pnorm(-abs(ate_mean / ate_se))))

if (ate_mean < 0) {
  cat("  • Interpretación: El tratamiento (renovable/secuestro) REDUCE emisiones\n")
} else {
  cat("  • Interpretación: El tratamiento (renovable/secuestro) AUMENTA emisiones\n")
}

# --- 3.4 Mapa de CATE ---
co2_sf_map <- st_as_sf(co2_ml, coords = c("longitude", "latitude"), crs = 4326)

p_cate <- ggplot() +
  geom_sf(data = co2_sf_map, aes(fill = cate), color = NA, size = 0.5) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850",
                       midpoint = 0, name = "CATE",
                       limits = quantile(co2_ml$cate, c(0.02, 0.98), na.rm = TRUE)) +
  labs(
    title = "Mapa de CATE - Efecto Heterogéneo del Tratamiento",
    subtitle = sprintf("ATE = %.3f | Causal Forest (grf)", ate_mean),
    caption = "Valores negativos = reducción de emisiones con tratamiento renovable"
  ) +
  theme(legend.position = "bottom")

ggsave("results/01_cate_map.png", p_cate, width = 10, height = 10, dpi = 300)
cat("  ✓ Mapa de CATE guardado\n")

# --- 3.5 Análisis de heterogeneidad por región ---
# Dividir por cuartiles de latitud (Norte/Sur) y longitud (Este/Oeste)
co2_ml <- co2_ml %>%
  mutate(
    region_lat = ifelse(latitude > median(latitude, na.rm = TRUE), "Norte", "Sur"),
    region_lon = ifelse(longitude > median(longitude, na.rm = TRUE), "Oriente", "Occidente"),
    region_combined = paste(region_lat, region_lon, sep = "-")
  )

heterogeneidad <- co2_ml %>%
  group_by(region_combined) %>%
  summarise(
    n = n(),
    cate_mean = mean(cate, na.rm = TRUE),
    cate_se = sd(cate, na.rm = TRUE) / sqrt(n),
    cate_ci_inf = cate_mean - 1.96 * cate_se,
    cate_ci_sup = cate_mean + 1.96 * cate_se,
    .groups = "drop"
  )

cat("\n  • Heterogeneidad regional de CATE:\n")
print(heterogeneidad)

p_het <- ggplot(heterogeneidad, aes(x = region_combined, y = cate_mean)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = cate_ci_inf, ymax = cate_ci_sup), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Heterogeneidad de CATE por Región",
    x = "Región", y = "CATE promedio",
    caption = "Barras = IC 95%"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/02_cate_heterogeneity.png", p_het, width = 8, height = 6, dpi = 300)
cat("  ✓ Gráfico de heterogeneidad guardado\n")

# ==============================================================================
# 4. RANDOM FORESTS PREDICTIVO ------------------------------------------------
# ==============================================================================
cat("\n=== 4. Random Forest predictivo ===\n")

# --- 4.1 Entrenar RF ---
cat("  • Entrenando Random Forest (ranger)...\n")
rf_model <- ranger(
  formula = log_co2 ~ sector_idx + sub_sector_idx + log_actividad + 
    lat + lon + lat2 + lon2 + lat_lon + densidad_local,
  data = co2_ml,
  num.trees = 500,
  mtry = floor(sqrt(ncol(X))),
  min.node.size = 5,
  importance = "impurity",
  num.threads = 1,
  write.forest = TRUE
)

cat("  ✓ Random Forest entrenado\n")

# --- 4.2 Predicciones y métricas ---
rf_pred <- predict(rf_model, data = co2_ml)$predictions

# Métricas
residuals_rf <- co2_ml$log_co2 - rf_pred
rmse_rf <- sqrt(mean(residuals_rf^2))
r2_rf <- 1 - sum(residuals_rf^2) / sum((co2_ml$log_co2 - mean(co2_ml$log_co2))^2)

cat(sprintf("\n  • RMSE: %.4f\n", rmse_rf))
cat(sprintf("  • R²: %.4f\n", r2_rf))

# ==============================================================================
# 5. MODELO LINEAL (BENCHMARK) ------------------------------------------------
# ==============================================================================
cat("\n=== 5. Modelo lineal (benchmark) ===\n")

lm_model <- lm(log_co2 ~ sector_idx + sub_sector_idx + log_actividad + 
                 lat + lon + lat2 + lon2 + lat_lon + densidad_local,
               data = co2_ml)

lm_pred <- predict(lm_model)
residuals_lm <- co2_ml$log_co2 - lm_pred
rmse_lm <- sqrt(mean(residuals_lm^2))
r2_lm <- 1 - sum(residuals_lm^2) / sum((co2_ml$log_co2 - mean(co2_ml$log_co2))^2)

cat(sprintf("  • RMSE: %.4f\n", rmse_lm))
cat(sprintf("  • R²: %.4f\n", r2_lm))

# Comparación
comparison <- data.frame(
  Modelo = c("Random Forest", "Modelo Lineal"),
  RMSE = c(rmse_rf, rmse_lm),
  R2 = c(r2_rf, r2_lm)
)

cat("\n  • Comparación de modelos:\n")
print(comparison)

write.csv(comparison, "results/03_model_comparison.csv", row.names = FALSE)

# ==============================================================================
# 6. SHAP VALUES E INTERPRETACIÓN ---------------------------------------------
# ==============================================================================
cat("\n=== 6. SHAP values para interpretación ===\n")

# --- 6.1 Calcular SHAP values para RF ---
# shapviz requiere un modelo con predict method
# Usamos una muestra para eficiencia
set.seed(2026)
idx_shap <- sample(1:nrow(co2_ml), min(1000, nrow(co2_ml)))
X_shap <- co2_ml[idx_shap, covariables_cols]

# Calcular SHAP con shapviz (usando el modelo ranger)
# Nota: shapviz funciona mejor con modelos tree-based
cat("  • Calculando SHAP values (muestra de 1000 obs)...\n")

# Alternativa: usar variable importance de ranger (más estable)
var_imp <- importance(rf_model)
var_imp_df <- data.frame(
  variable = names(var_imp),
  importance = var_imp
) %>%
  arrange(desc(importance)) %>%
  mutate(variable = factor(variable, levels = variable))

cat("\n  • Importancia de variables (Gini):\n")
print(var_imp_df)

p_imp <- ggplot(var_imp_df, aes(x = importance, y = variable)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Importancia de Variables (Random Forest)",
    x = "Importancia (Gini)", y = "Variable"
  ) +
  theme_minimal()

ggsave("results/04_variable_importance.png", p_imp, width = 8, height = 6, dpi = 300)
cat("  ✓ Gráfico de importancia guardado\n")

# --- 6.2 SHAP values con shapviz (si está disponible) ---
if (requireNamespace("shapviz", quietly = TRUE)) {
  # shapviz necesita un objeto xgb.Booster o similar
  # Para ranger, usamos una aproximación con predicciones
  cat("  • shapviz disponible - calculando SHAP approximation...\n")
  
  # Usar predicciones del RF como proxy
  sv <- shapviz(rf_model, X_pred = X_shap, X = X_shap)
  
  p_shap_imp <- sv_importance(sv, show_numbers = TRUE)
  p_shap_dep <- sv_dependence(sv, v = "log_actividad", color_var = "lat")
  
  ggsave("results/05_shap_importance.png", p_shap_imp, width = 8, height = 6, dpi = 300)
  ggsave("results/06_shap_dependence.png", p_shap_dep, width = 8, height = 6, dpi = 300)
  cat("  ✓ Gráficos SHAP guardados\n")
} else {
  cat("  ⚠ shapviz no disponible - usando variable importance de ranger\n")
}

# ==============================================================================
# 7. RESUMEN Y ENTREGABLES ----------------------------------------------------
# ==============================================================================
cat("\n=== 7. Resumen de entregables ===\n")

resumen_ml <- data.frame(
  Metrica = c(
    "N observaciones",
    "N variables",
    "ATE (Causal Forest)",
    "SE del ATE",
    "IC 95% ATE (inf)",
    "IC 95% ATE (sup)",
    "RMSE Random Forest",
    "R² Random Forest",
    "RMSE Modelo Lineal",
    "R² Modelo Lineal"
  ),
  Valor = c(
    nrow(co2_ml),
    ncol(X),
    round(ate_mean, 4),
    round(ate_se, 4),
    round(ate_ci[1], 4),
    round(ate_ci[2], 4),
    round(rmse_rf, 4),
    round(r2_rf, 4),
    round(rmse_lm, 4),
    round(r2_lm, 4)
  )
)

write.csv(resumen_ml, "results/07_ml_summary.csv", row.names = FALSE)

cat("\n  Archivos generados:\n")
cat("  ✓ results/01_cate_map.png              - Mapa de CATE\n")
cat("  ✓ results/02_cate_heterogeneity.png    - Heterogeneidad regional\n")
cat("  ✓ results/03_model_comparison.csv      - Comparación RF vs LM\n")
cat("  ✓ results/04_variable_importance.png   - Importancia de variables\n")
cat("  ✓ results/05_shap_importance.png       - SHAP importance (si disponible)\n")
cat("  ✓ results/06_shap_dependence.png       - SHAP dependence (si disponible)\n")
cat("  ✓ results/07_ml_summary.csv            - Resumen cuantitativo\n")

# Guardar modelo y datos para reproducibilidad
saveRDS(cf_model, "data/causal_forest_model.rds")
saveRDS(rf_model, "data/random_forest_model.rds")
cat("  ✓ Modelos guardados en data/\n")

cat("\n=== Machine Learning completado ===\n")