# ==============================================================================
# Script: 03_spatial_descriptive.R
# Proyecto: Análisis Espacial de Emisiones de CO2 en Colombia
# Sección: 1. Estadística Espacial
# Objetivo: Caracterizar la dependencia espacial de emisiones de CO2 y modelar
#           su variabilidad a múltiples escalas.
# Referencias: Anselin (1995); Fotheringham et al. (2002).
# ==============================================================================

# 0. CONFIGURACIÓN Y LIBRERÍAS ------------------------------------------------
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(spdep)        # Moran global/local, matrices de pesos
  library(gstat)        # Variogramas y kriging
  library(spatialreg)   # Modelos SAR/SEM/SLX
  library(viridis)
  library(patchwork)    # Para combinar gráficos
})

# Fallback: si sfdep está disponible, úsalo para LISA (más rápido)
if (requireNamespace("sfdep", quietly = TRUE)) {
  library(sfdep)
  USE_SFDEP <- TRUE
} else {
  USE_SFDEP <- FALSE
  message("sfdep no disponible. Usando spdep para LISA (más lento).")
}

set.seed(2026)
theme_set(theme_minimal())

# Crear carpetas de salida
dir.create("results", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# 1. CARGA Y PREPARACIÓN DE DATOS --------------------------------------------
cat("=== 1. Cargando y preparando datos ===\n")

co2_raw <- read.csv("C:/Users/monic/Downloads/emisiones_co2_espacial/emisiones_co2_espacial/datos/climate_trace_COL_2023.csv", stringsAsFactors = FALSE)

# Filtrar solo coordenadas válidas en Colombia
co2_clean <- co2_raw %>%
  filter(!is.na(lon), !is.na(lat),
         lon > -82 & lon < -66,
         lat > -5 & lat < 13) %>%
  # Agrupar por punto exacto y sumar emisiones
  group_by(lon, lat, sector, subsector) %>%
  summarise(
    co2_total = sum(co2, na.rm = TRUE),
    activity  = sum(activity, na.rm = TRUE),
    .groups = "drop"
  )

co2_sf <- st_as_sf(co2_clean, coords = c("lon", "lat"), crs = 4326)

# Proyectar a un CRS métrico para análisis de distancia (MAGNA-SIRGAS / Colombia)
# EPSG:9377 es el oficial para Colombia
co2_sf <- st_transform(co2_sf, 9377)

N <- nrow(co2_sf)
cat(sprintf("  • %d puntos únicos de emisión cargados\n", N))
cat(sprintf("  • Rango de CO2: [%.2f, %.2f] t\n", 
            min(co2_sf$co2_total), max(co2_sf$co2_total)))

# 2. ANÁLISIS DE SENSIBILIDAD MAUP (0.1° vs 0.5°) ----------------------------
cat("\n=== 2. Análisis MAUP: sensibilidad a la escala de agregación ===\n")

# Función para agregar datos en una malla y calcular Moran's I
maup_analysis <- function(sf_data, cellsize_deg, label) {
  # Volver a WGS84 temporalmente para crear grid en grados
  sf_wgs <- st_transform(sf_data, 4326)
  bbox <- st_bbox(sf_wgs)
  
  # Crear malla regular
  grid <- st_make_grid(sf_wgs, cellsize = cellsize_deg, square = TRUE)
  grid_sf <- st_sf(id = 1:length(grid), geometry = grid)
  
  # Intersectar y agregar
  joined <- st_join(grid_sf, sf_wgs, join = st_intersects)
  agg <- joined %>%
    group_by(id) %>%
    summarise(co2_agg = sum(co2_total, na.rm = TRUE), .groups = "drop")
  
  # Filtrar celdas con datos
  agg <- agg[!is.na(agg$co2_agg) & agg$co2_agg > 0, ]
  
  # Matriz de pesos (reina = 8 vecinos)
  nb <- poly2nb(agg, queen = TRUE)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Moran's I global
  moran_res <- moran.test(agg$co2_agg, lw, zero.policy = TRUE)
  
  list(
    label = label,
    cellsize = cellsize_deg,
    n_cells = nrow(agg),
    moran_I = moran_res$estimate[1],
    p_value = moran_res$p.value[1],
    data = agg,
    listw = lw
  )
}

maup_01 <- maup_analysis(co2_sf, 0.1, "0.1° (~11 km)")
maup_05 <- maup_analysis(co2_sf, 0.5, "0.5° (~55 km)")

cat(sprintf("  • Resolución 0.1°: %d celdas | Moran's I = %.4f (p = %.4f)\n",
            maup_01$n_cells, maup_01$moran_I, maup_01$p_value))
cat(sprintf("  • Resolución 0.5°: %d celdas | Moran's I = %.4f (p = %.4f)\n",
            maup_05$n_cells, maup_05$moran_I, maup_05$p_value))

# Visualización MAUP
df_maup <- data.frame(
  resolucion = c(maup_01$label, maup_05$label),
  n_celdas   = c(maup_01$n_cells, maup_05$n_cells),
  moran_I    = c(maup_01$moran_I, maup_05$moran_I)
)

p_maup <- ggplot(df_maup, aes(x = resolucion, y = moran_I, fill = resolucion)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("I=%.3f\nn=%d", moran_I, n_celdas)), 
            vjust = -0.5, size = 4) +
  coord_cartesian(ylim = c(0, max(df_maup$moran_I) * 1.3)) +
  scale_fill_viridis_d(begin = 0.3, end = 0.8) +
  labs(
    title = "Efecto MAUP: Moran's I a diferentes resoluciones",
    subtitle = "Agregación de emisiones de CO2 en grids de 0.1° y 0.5°",
    x = "Resolución espacial", y = "Índice de Moran Global (I)",
    fill = "Resolución"
  ) +
  theme(legend.position = "none")

ggsave("results/01_maup_moran.png", p_maup, width = 8, height = 6, dpi = 300)

# 3. ÍNDICE DE MORAN GLOBAL (a nivel de punto) -------------------------------
cat("\n=== 3. Autocorrelación espacial global (nivel punto) ===\n")

# Limpiar puntos duplicados primero
co2_sf <- co2_sf %>%
  filter(!is.na(log_co2), is.finite(log_co2)) %>%
  distinct(geometry, .keep_all = TRUE)  # Eliminar coordenadas idénticas

cat("Observaciones después de limpiar:", nrow(co2_sf), "\n")

# Recalcular matriz de pesos
coords <- st_coordinates(co2_sf)
nb_knn <- knn2nb(knearneigh(coords, k = 8))
lw_knn <- nb2listw(nb_knn, style = "W", zero.policy = TRUE)

# Moran's I con permutaciones Monte Carlo (999 permutaciones)
set.seed(2026)
moran_global <- moran.mc(co2_sf$log_co2, lw_knn, nsim = 999, zero.policy = TRUE)

cat("Moran's I observado:", moran_mc$statistic, "\n")
cat("p-valor (permutaciones):", moran_mc$p.value, "\n")

# Interpretación
if (moran_mc$p.value < 0.05) {
  cat("✓ Autocorrelación espacial significativa\n")
} else {
  cat("✗ No hay evidencia de autocorrelación espacial\n")
}
cat(sprintf("  • Moran's I = %.4f\n", moran_global$estimate[1]))
cat(sprintf("  • Estadístico Z = %.4f\n", moran_global$estimate[3]))
cat(sprintf("  • p-valor = %.2e\n", moran_global$p.value[1]))

# 4. LISA: ÍNDICE DE MORAN LOCAL Y MAPA DE CLUSTERS --------------------------
cat("\n=== 4. LISA: Identificación de hot/cold spots ===\n")

# Usar spdep directamente (más estable)
coords <- st_coordinates(co2_sf)
nb_knn <- knn2nb(knearneigh(coords, k = 8))
lw_knn <- nb2listw(nb_knn, style = "W", zero.policy = TRUE)

# LISA con spdep
lisa_res <- localmoran(co2_sf$co2_total, lw_knn, zero.policy = TRUE)
co2_sf$lisa_Ii   <- lisa_res[, "Ii"]
co2_sf$lisa_pval <- lisa_res[, "Pr(z != E(Ii))"]

# Clasificación de clusters (Anselin 1995) 
# 1. Verificar y limpiar valores problemáticos
co2_sf <- co2_sf %>%
  filter(!is.na(co2_total), 
         is.finite(co2_total),
         co2_total > 0)  # Eliminar ceros y negativos si es necesario

# 2. Recalcular la matriz de pesos después de filtrar
coords <- st_coordinates(co2_sf)
nb_knn <- knn2nb(knearneigh(coords, k = 8))
lw_knn <- nb2listw(nb_knn, style = "W", zero.policy = TRUE)

# 3. Ahora sí, crear la variable estandarizada y el lag
co2_sf <- co2_sf %>%
  mutate(
    co2_z = as.numeric(scale(co2_total)),  # Convertir a vector numérico
  )

# 4. Verificar que no haya NAs después de escalar
co2_sf <- co2_sf %>%
  filter(!is.na(co2_z))

# 5. Calcular el lag espacial
co2_sf$lag_z <- lag.listw(lw_knn, co2_sf$co2_z, zero.policy = TRUE)

# 6. Finalmente, la clasificación de clusters
co2_sf <- co2_sf %>%
  mutate(
    cluster = case_when(
      is.na(lisa_pval) | lisa_pval > 0.05 ~ "No significativo",
      co2_z > 0 & lag_z > 0 ~ "High-High (Hot Spot)",
      co2_z < 0 & lag_z < 0 ~ "Low-Low (Cold Spot)",
      co2_z > 0 & lag_z < 0 ~ "High-Low (Outlier)",
      co2_z < 0 & lag_z > 0 ~ "Low-High (Outlier)",
      TRUE ~ "No significativo"
    ),
    cluster = factor(cluster, levels = c(
      "High-High (Hot Spot)", 
      "Low-Low (Cold Spot)",
      "High-Low (Outlier)", 
      "Low-High (Outlier)", 
      "No significativo"
    ))
  )

# Resumen de clusters
cat("\nDistribución de clusters LISA:\n")
print(table(co2_sf$cluster))

# Mapa de clusters LISA
p_lisa <- ggplot() +
  geom_sf(data = co2_sf, aes(color = cluster), size = 0.8, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "High-High (Hot Spot)" = "#d73027",
      "Low-Low (Cold Spot)"  = "#4575b4",
      "High-Low (Outlier)"   = "#fee090",
      "Low-High (Outlier)"   = "#abd9e9",
      "No significativo"     = "#d9d9d9"
    ),
    name = "Cluster LISA"
  ) +
  labs(
    title = "Mapa de Clusters LISA - Emisiones de CO2 en Colombia",
    subtitle = sprintf("Moran's I global = %.3f | k-vecinos = 8", 
                       moran_global$estimate[1]),
    caption = "Clasificación según Anselin (1995)"
  ) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

ggsave("results/02_lisa_clusters.png", p_lisa, width = 10, height = 10, dpi = 300)

# Mapa de significancia (p-values)
p_signif <- ggplot() +
  geom_sf(data = co2_sf, aes(fill = -log10(lisa_pval + 1e-10)), 
          color = NA, size = 0.5) +
  scale_fill_viridis(option = "plasma", name = "-log10(p)") +
  labs(
    title = "Significancia estadística LISA",
    subtitle = "Valores altos indican clusters espacialmente significativos"
  ) +
  theme(legend.position = "bottom")

ggsave("results/03_lisa_significance.png", p_signif, width = 10, height = 10, dpi = 300)

# 5. SEMIVARIOGRAMA Y MODELOS GEOESTADÍSTICOS POR REGIÓN ---------------------
cat("\n=== 5. Semivariograma y modelos geoestadísticos ===\n")

# Asignar región natural aproximada por cuantiles de coordenadas
# (En producción, usar shapefile oficial de regiones naturales)
coords_wgs <- st_coordinates(st_transform(co2_sf, 4326))
co2_sf$region_approx <- cut(coords_wgs[,2], 
                            breaks = quantile(coords_wgs[,2], probs = seq(0, 1, 0.25)),
                            labels = c("Sur", "Sur-Centro", "Norte-Centro", "Norte"),
                            include.lowest = TRUE)

# Variograma global
vgm_global <- variogram(co2_total ~ 1, data = co2_sf, cutoff = 500000, width = 25000)

# Ajustar modelo exponencial (típico para procesos ambientales)
fit_vgm_global <- fit.variogram(vgm_global, model = vgm(
  psill = var(co2_sf$co2_total) * 0.7,
  model = "Exp",
  range = 100000,
  nugget = var(co2_sf$co2_total) * 0.3
))

cat("\nModelo de variograma global ajustado:\n")
print(fit_vgm_global)

# Variograma por región
regiones <- unique(co2_sf$region_approx[!is.na(co2_sf$region_approx)])
vgm_list <- list()
fit_list <- list()

for (r in regiones) {
  sub <- co2_sf[co2_sf$region_approx == r, ]
  if (nrow(sub) < 30) next  # Mínimo para ajustar
  
  v <- variogram(co2_total ~ 1, data = sub, cutoff = 300000, width = 15000)
  f <- tryCatch(
    fit.variogram(v, model = vgm(psill = var(sub$co2_total)*0.7, 
                                  model = "Exp", range = 80000, 
                                  nugget = var(sub$co2_total)*0.3)),
    error = function(e) NULL
  )
  vgm_list[[r]] <- v
  fit_list[[r]] <- f
}

# Gráfico comparativo de variogramas por región
df_vgm <- data.frame()
for (r in names(vgm_list)) {
  v <- vgm_list[[r]]
  v$region <- r
  df_vgm <- rbind(df_vgm, as.data.frame(v))
}

p_vgm <- ggplot(df_vgm, aes(x = dist, y = gamma, color = region)) +
  geom_point(size = 1.5, alpha = 0.6) +
  labs(
    title = "Semivariogramas empíricos por región",
    subtitle = "Modelo exponencial ajustado por región natural",
    x = "Distancia (metros)", y = "Semivarianza (γ)",
    color = "Región"
  ) +
  scale_color_viridis_d() +
  theme(legend.position = "bottom")

ggsave("results/04_variogram_by_region.png", p_vgm, width = 10, height = 6, dpi = 300)

# 6. KRIGING: SUPERFICIE DE INTENSIDAD DE CARBONO A 1 KM ---------------------
cat("\n=== 6. Kriging ordinario - Superficie a 1 km ===\n")

# Crear malla de predicción a 1 km sobre el bounding box
bbox_proj <- st_bbox(co2_sf)
x_seq <- seq(bbox_proj["xmin"], bbox_proj["xmax"], by = 1000)
y_seq <- seq(bbox_proj["ymin"], bbox_proj["ymax"], by = 1000)
pred_grid <- expand.grid(x = x_seq, y = y_seq)
pred_grid_sf <- st_as_sf(pred_grid, coords = c("x", "y"), crs = 9377)

# Kriging ordinario usando el modelo global ajustado
cat("  • Ejecutando kriging ordinario... (esto puede tomar unos minutos)\n")
krig_result <- krige(
  co2_total ~ 1,
  locations = co2_sf,
  newdata = pred_grid_sf,
  model = fit_vgm_global,
  nmax = 30  # Limitar vecinos para eficiencia
)

# Extraer predicción y varianza
krig_result$pred_log <- log1p(krig_result$var1.pred)  # Transformación log para visualización

p_krig <- ggplot() +
  geom_sf(data = krig_result, aes(fill = pred_log), color = NA) +
  scale_fill_viridis(option = "magma", name = "log(CO2+1) [t]") +
  geom_sf(data = co2_sf, color = "white", size = 0.3, alpha = 0.3) +
  labs(
    title = "Superficie de intensidad de carbono (Kriging Ordinario)",
    subtitle = "Resolución: 1 km | Modelo: Exponencial",
    caption = "Predicción basada en datos de Climate Trace 2023"
  ) +
  theme(legend.position = "bottom")

ggsave("results/05_kriging_surface_1km.png", p_krig, width = 10, height = 12, dpi = 300)

# 7. MODELOS ESPACIALES: SAR / SEM / SLX -------------------------------------
cat("\n=== 7. Modelos de regresión espacial (SAR/SEM/SLX) ===\n")

# Variable respuesta: log(CO2) para normalizar
co2_sf$log_co2 <- log1p(co2_sf$co2_total)

# Variable explicativa sintética: densidad de actividad (proxy)
# En un análisis real, aquí irían variables como población, PIB, uso de suelo
co2_sf$log_activity <- log1p(co2_sf$activity)

# Lista de pesos para modelos (normalizada por fila)
lw_model <- nb2listw(nb_knn, style = "W")

# --- Modelo SAR (Spatial Autoregressive / Lag) ---
model_sar <- lagsarlm(log_co2 ~ log_activity, data = co2_sf, listw = lw_model)
cat("\n--- Modelo SAR (Lag espacial) ---\n")
print(summary(model_sar))

# --- Modelo SEM (Spatial Error) ---
model_sem <- errorsarlm(log_co2 ~ log_activity, data = co2_sf, listw = lw_model)
cat("\n--- Modelo SEM (Error espacial) ---\n")
print(summary(model_sem))

# --- Modelo SLX (Spatial Lag of X) ---
# SLX incluye los valores rezagados de las variables explicativas
model_slx <- lm.LMtests(log_co2 ~ log_activity, data = co2_sf, listw = lw_model, 
                        test = "all")
cat("\n--- Tests LM para selección de modelo ---\n")
print(model_slx)

# Comparación de modelos por AIC
aic_table <- data.frame(
  Modelo = c("SAR (Lag)", "SEM (Error)", "OLS (referencia)"),
  AIC = c(AIC(model_sar), AIC(model_sem), 
          AIC(lm(log_co2 ~ log_activity, data = co2_sf)))
)
cat("\nComparación AIC (menor es mejor):\n")
print(aic_table)

# Guardar tabla de resultados
write.csv(aic_table, "results/06_model_comparison_aic.csv", row.names = FALSE)

# 8. RESUMEN FINAL Y ENTREGABLES ---------------------------------------------
cat("\n=== 8. Resumen de entregables ===\n")
cat("✓ results/01_maup_moran.png           - Análisis MAUP\n")
cat("✓ results/02_lisa_clusters.png        - Mapa de clusters LISA\n")
cat("✓ results/03_lisa_significance.png    - Significancia LISA\n")
cat("✓ results/04_variogram_by_region.png  - Variogramas por región\n")
cat("✓ results/05_kriging_surface_1km.png  - Superficie de kriging 1km\n")
cat("✓ results/06_model_comparison_aic.csv - Comparación SAR/SEM/SLX\n")

# Guardar objeto sf con resultados LISA para uso posterior
st_write(co2_sf, "data/co2_lisa_results.gpkg", delete_layer = TRUE, quiet = TRUE)

cat("\n=== Análisis espacial descriptivo completado ===\n")