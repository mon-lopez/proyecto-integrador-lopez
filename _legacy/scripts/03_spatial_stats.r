#!/usr/bin/env Rscript
# ============================================================================
# 03_spatial_stats.R — Análisis Exploratorio de Datos Espaciales (ESDA)
# ============================================================================
# Proyecto : Análisis Espacial de CO₂ y Evaluación Cuasi-Experimental de FNCER
#            en Colombia
# Autor    : Grupo de Investigación en Estadística Aplicada
# Fecha    : 2026-11-17
# Propósito:
#   1. Cargar la muestra GRTS de 1750 celdas (5km x 5km)
#   2. Construir matriz de vecindad espacial (reina, k-vecinos más cercanos)
#   3. Calcular Índice de Moran global para densidad de emisiones
#   4. Calcular Índice de Moran local (LISA) para identificar clusters
#   5. Elaborar mapas de hot/cold spots
#   6. Ajustar semivariograma empírico
#   7. Exportar resultados y visualizaciones
#
# Uso      : Rscript scripts/03_spatial_stats.R
#            source("scripts/03_spatial_stats.R")
#
# Entradas : datos/processed/muestra_grts_1750celdas.rds
#            datos/processed/celdas_poblacion.rds
# Salidas  : output/maps/moran_local_lisa.png
#            output/maps/moran_scatterplot.png
#            output/figures/semivariograma.png
#            datos/processed/estadistica_espacial.rds
# ============================================================================
# ----------------------------------------------------------------------------
# 0. Cargar configuración del proyecto
# ----------------------------------------------------------------------------
source(file.path(Sys.getenv("CO2_PROJECT_ROOT", getwd()),
"scripts", "00_setup.R"))
library("tidyverse")
library("sf")
library("spdep")
library("sfdep")
library("gstat")
log_msg("Iniciando ANÁLISIS EXPLORATORIO DE DATOS ESPACIALES (ESDA)", "INFO")
inicio <- Sys.time()
# ----------------------------------------------------------------------------
# 0.1. Detectar carpeta de datos procesados (flexible)
# ----------------------------------------------------------------------------
carpeta_proc <- c(
file.path(RAIZ_PROYECTO, "datos", "processed"),
file.path(RAIZ_PROYECTO, "datos", "datos_procesados")
)
carpeta_proc <- carpeta_proc[dir.exists(carpeta_proc)][1]
if (is.na(carpeta_proc)) {
stop("No se encontró carpeta de datos procesados. Ejecuta 02_spatial_sampling.R primero.")
}
log_msg(paste("Carpeta de datos procesados:", carpeta_proc), "INFO")
# ============================================================================
# MÓDULO A — CARGA DE DATOS ESPACIALES
# ============================================================================
log_msg("MÓDULO A: Cargando datos espaciales", "INFO")
# --- A1. Cargar muestra GRTS -------------------------------------------------
path_muestra <- file.path(carpeta_proc, "muestra_grts_1750celdas.rds")
if (!file.exists(path_muestra)) {
stop("No se encontró la muestra GRTS. Ejecuta 02_spatial_sampling.R primero.")
}
muestra_sf <- readRDS(path_muestra)
log_msg(sprintf("A1: Muestra cargada: %d celdas × %d columnas",
nrow(muestra_sf), ncol(muestra_sf)), "OK")
# --- A2. Verificar que tiene geometría y variable de interés -----------------
stopifnot(inherits(muestra_sf, "sf"))
stopifnot("densidad_emisiones" %in% names(muestra_sf))
stopifnot(!is.null(st_geometry(muestra_sf)))
# --- A3. Cargar población completa (opcional, para contexto) -----------------
path_poblacion <- file.path(carpeta_proc, "celdas_poblacion.rds")
if (file.exists(path_poblacion)) {
poblacion_sf <- readRDS(path_poblacion)
log_msg(sprintf("A3: Población cargada: %d celdas", nrow(poblacion_sf)), "OK")
} else {
poblacion_sf <- NULL
log_msg("A3: Población completa no disponible (opcional)", "WARN")
}
# --- A4. Resumen estadístico de la variable de interés -----------------------
resumen_var <- muestra_sf %>%
st_drop_geometry() %>%
summarise(
n = n(),
media = mean(densidad_emisiones, na.rm = TRUE),
mediana = median(densidad_emisiones, na.rm = TRUE),
sd = sd(densidad_emisiones, na.rm = TRUE),
min = min(densidad_emisiones, na.rm = TRUE),
max = max(densidad_emisiones, na.rm = TRUE),
pct_na = 100 * sum(is.na(densidad_emisiones)) / n()
)
log_msg("A4: Resumen de densidad_emisiones:", "INFO")
print(resumen_var)
# ============================================================================
# MÓDULO B — CONSTRUCCIÓN DE MATRIZ DE VECINDAD ESPACIAL
# ============================================================================
log_msg("MÓDULO B: Construyendo matriz de vecindad espacial", "INFO")
# --- B1. Vecindad tipo reina (Queen contiguity) ------------------------------
log_msg("B1: Calculando vecindad tipo reina...", "INFO")
nb_reina <- poly2nb(muestra_sf, queen = TRUE)
log_msg(sprintf("B1: Vecindad reina: %d vecinos en total",
sum(card(nb_reina))), "OK")
# --- B2. Matriz de pesos espaciales (W) --------------------------------------
log_msg("B2: Construyendo matriz de pesos espaciales...", "INFO")
W_reina <- nb2listw(nb_reina, style = "W", zero.policy = TRUE)
log_msg("B2: Matriz W estandarizada por fila (style = 'W')", "OK")
# --- B3. Vecindad k-NN (backup, k=4) -----------------------------------------
log_msg("B3: Calculando vecindad k-NN (k=4)...", "INFO")
coords_centroides <- st_coordinates(st_centroid(muestra_sf))
nb_knn <- knearneigh(coords_centroides, k = 4)
nb_knn_list <- knn2nb(nb_knn)
W_knn <- nb2listw(nb_knn_list, style = "W", zero.policy = TRUE)
log_msg("B3: Matriz W k-NN (k=4) construida", "OK")
# ============================================================================
# MÓDULO C — ÍNDICE DE MORAN GLOBAL
# ============================================================================
log_msg("MÓDULO C: Calculando Índice de Moran global", "INFO")
# --- C1. Moran global con vecindad reina -------------------------------------
log_msg("C1: Moran global (reina)...", "INFO")
moran_test_reina <- moran.test(
muestra_sf$densidad_emisiones,
W_reina,
)
log_msg(sprintf("C1: Moran I = %.4f | E(I) = %.4f | p-value = %.6f",
moran_test_reina$estimate["Moran I"],
moran_test_reina$estimate["Expectation"],
moran_test_reina$p.value), "OK")
# --- C2. Moran global con k-NN -----------------------------------------------
log_msg("C2: Moran global (k-NN)...", "INFO")
moran_test_knn <- moran.test(
muestra_sf$densidad_emisiones,
W_knn,
)
log_msg(sprintf("C2: Moran I = %.4f | E(I) = %.4f | p-value = %.6f",
moran_test_knn$estimate["Moran I"],
moran_test_knn$estimate["Expectation"],
moran_test_knn$p.value), "OK")
# --- C3. Interpretación ------------------------------------------------------
interpretacion_moran <- if (moran_test_reina$p.value < 0.05) {
if (moran_test_reina$estimate["Moran I"] > 0) {
"Autocorrelación espacial POSITIVA significativa: valores similares tienden a agruparse"
} else {
"Autocorrelación espacial NEGATIVA significativa: valores disímiles tienden a agruparse"
}
} else {
"No hay evidencia de autocorrelación espacial significativa (p > 0.05)"
}
log_msg(paste("C3:", interpretacion_moran), "INFO")
# ============================================================================
# MÓDULO D — ÍNDICE DE MORAN LOCAL (LISA)
# ============================================================================
log_msg("MÓDULO D: Calculando Índice de Moran local (LISA)", "INFO")
# --- D1. Calcular estadísticos locales ---------------------------------------
log_msg("D1: Calculando LISA con vecindad reina...", "INFO")
local_moran <- localmoran(
muestra_sf$densidad_emisiones,
W_reina,
zero.policy = TRUE
)
# --- D2. Añadir resultados al objeto sf --------------------------------------
muestra_sf <- muestra_sf %>%
  mutate(
    Ii = local_moran[, "Ii"],
    Ii_pvalue = local_moran[, "Pr(z != E(Ii))"],
    Ii_zscore = local_moran[, "Z.Ii"],  # Corregido: "Z.Ii" en lugar de "Z(Ii)"
    # Clasificar clusters LISA
    cluster_lisa = case_when(
      Ii > 0 & Ii_pvalue < 0.05 & densidad_emisiones > mean(densidad_emisiones) ~ "HH (Hot Spot)",
      Ii > 0 & Ii_pvalue < 0.05 & densidad_emisiones <= mean(densidad_emisiones) ~ "LL (Cold Spot)",
      Ii < 0 & Ii_pvalue < 0.05 ~ "Outlier",
      TRUE ~ "No significativo"
    )
  )
# --- D3. Resumen de clusters -------------------------------------------------
resumen_clusters <- muestra_sf %>%
st_drop_geometry() %>%
count(cluster_lisa, name = "n_celdas") %>%
mutate(pct = round(100 * n_celdas / sum(n_celdas), 1))
log_msg("D3: Distribución de clusters LISA:", "INFO")
print(resumen_clusters)
# ============================================================================
# MÓDULO E — SEMIVARIOGRAMA EMPÍRICO
# ============================================================================
log_msg("MÓDULO E: Calculando semivariograma empírico", "INFO")
# --- E1. Convertir a objeto SpatialPolygonsDataFrame (requerido por gstat) ----
log_msg("E1: Convirtiendo a formato compatible con gstat...", "INFO")
muestra_sp <- as(muestra_sf, "Spatial")
# --- E2. Calcular semivariograma empírico ------------------------------------
log_msg("E2: Calculando semivariograma...", "INFO")
variograma_emp <- variogram(
densidad_emisiones ~ 1,
data = muestra_sp,
cutoff = 50000,  # 50 km máximo
width = 2500     # bins de 2.5 km
)
log_msg(sprintf("E2: Semivariograma calculado: %d puntos", nrow(variograma_emp)), "OK")
# --- E3. Ajustar modelo esférico ---------------------------------------------
log_msg("E3: Ajustando modelo esférico...", "INFO")
modelo_esf <- fit.variogram(
variograma_emp,
model = vgm(
psill = var(muestra_sf$densidad_emisiones, na.rm = TRUE),
model = "Sph",
range = 20000,  # 20 km
nugget = 0
)
)
log_msg(sprintf("E3: Modelo ajustado | Nugget: %.2f | Sill: %.2f | Range: %.0f m",
modelo_esf$psill[1],
modelo_esf$psill[1] + modelo_esf$psill[2],
modelo_esf$range[2]), "OK")
# ============================================================================
# MÓDULO F — VISUALIZACIÓN ESPACIAL
# ============================================================================
log_msg("MÓDULO F: Generando visualizaciones espaciales", "INFO")
# --- F1. Mapa de densidad de emisiones ---------------------------------------
log_msg("F1: Mapa de densidad de emisiones...", "INFO")
mapa_densidad <- ggplot() +
geom_sf(data = muestra_sf, aes(fill = densidad_emisiones),
color = NA, alpha = 0.8) +
scale_fill_viridis_c(
name = "Densidad emisiones\n(t CO₂/km²)",
option = "plasma",
direction = 1
) +
labs(
title = "Densidad de Emisiones de CO₂ por Celda (5km)",
subtitle = sprintf("Muestra GRTS: %d celdas | Media: %.1f t/km²",
nrow(muestra_sf), mean(muestra_sf$densidad_emisiones, na.rm = TRUE))
) +
theme_minimal(base_size = 12) +
theme(
plot.title = element_text(face = "bold", size = 14),
legend.position = "right"
)
path_mapa_densidad <- file.path(PATHS$output, "maps", "mapa_densidad_emisiones_03.png")
ggsave(path_mapa_densidad, mapa_densidad, width = 12, height = 8, dpi = 300)
log_msg(sprintf("F1: Mapa guardado: %s", path_mapa_densidad), "OK")
# --- F2. Mapa de clusters LISA -----------------------------------------------
log_msg("F2: Mapa de clusters LISA...", "INFO")
mapa_lisa <- ggplot() +
geom_sf(data = muestra_sf, aes(fill = cluster_lisa),
color = "white", size = 0.2, alpha = 0.85) +
scale_fill_manual(
name = "Cluster LISA",
values = c(
"HH (Hot Spot)" = "#d73027",
"LL (Cold Spot)" = "#3288bd",
"Outlier" = "#fee08b",
"No significativo" = "#f0f0f0"
),
na.value = "#f0f0f0"
) +
labs(
title = "Análisis Local de Autocorrelación Espacial (LISA)",
subtitle = sprintf("Hot/Cold Spots | Moran I = %.3f (p = %.4f)",
moran_test_reina$estimate["Moran I"],
moran_test_reina$p.value)
) +
theme_minimal(base_size = 12) +
theme(
plot.title = element_text(face = "bold", size = 14),
legend.position = "right"
)
path_mapa_lisa <- file.path(PATHS$output, "maps", "moran_local_lisa.png")
ggsave(path_mapa_lisa, mapa_lisa, width = 12, height = 8, dpi = 300)
log_msg(sprintf("F2: Mapa LISA guardado: %s", path_mapa_lisa), "OK")
# --- F3. Scatterplot de Moran ------------------------------------------------
log_msg("F3: Scatterplot de Moran...", "INFO")
# Calcular valores estandarizados
z <- scale(muestra_sf$densidad_emisiones)
z_lag <- lag.listw(W_reina, z)
df_scatter <- tibble(
z = as.numeric(z),
z_lag = as.numeric(z_lag),
cluster = muestra_sf$cluster_lisa
)
scatter_moran <- ggplot(df_scatter, aes(x = z, y = z_lag, color = cluster)) +
geom_point(alpha = 0.6, size = 2) +
geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
scale_color_manual(
name = "Cluster",
values = c(
"HH (Hot Spot)" = "#d73027",
"LL (Cold Spot)" = "#3288bd",
"Outlier" = "#fee08b",
"No significativo" = "#999999"
)
) +
labs(
title = "Scatterplot de Moran",
subtitle = sprintf("Moran I = %.3f | Pendiente = %.3f",
moran_test_reina$estimate["Moran I"],
coef(lm(z_lag ~ z, data = df_scatter))[2]),
x = "Densidad emisiones (estandarizada)",
y = "Lag espacial (promedio vecinos)"
) +
theme_minimal(base_size = 12) +
theme(
plot.title = element_text(face = "bold", size = 14),
legend.position = "right"
)
path_scatter <- file.path(PATHS$output, "maps", "moran_scatterplot.png")
ggsave(path_scatter, scatter_moran, width = 10, height = 8, dpi = 300)
log_msg(sprintf("F3: Scatterplot guardado: %s", path_scatter), "OK")
# --- F4. Semivariograma ------------------------------------------------------
log_msg("F4: Gráfico de semivariograma...", "INFO")
df_vario <- as.data.frame(variograma_emp)

# Extraer parámetros del modelo ajustado
nugget   <- modelo_esf$psill[1]
sill     <- modelo_esf$psill[1] + modelo_esf$psill[2]
range_m  <- modelo_esf$range[2]

# Calcular curva teórica del modelo esférico manualmente
dist_seq <- seq(0, max(df_vario$dist), length.out = 100)
gamma_teo <- ifelse(
  dist_seq <= range_m,
  nugget + (sill - nugget) * (1.5 * (dist_seq / range_m) - 0.5 * (dist_seq / range_m)^3),
  sill
)
df_modelo <- data.frame(dist = dist_seq, gamma = gamma_teo)

p_vario <- ggplot() +
  geom_point(data = df_vario, aes(x = dist, y = gamma),
             size = 3, color = "#3288bd") +
  geom_line(data = df_vario, aes(x = dist, y = gamma),
            linewidth = 1, color = "#3288bd") +
  geom_line(data = df_modelo, aes(x = dist, y = gamma),
            linewidth = 1, color = "#d73027", linetype = "dashed") +
  labs(
    title = "Semivariograma Empírico",
    subtitle = sprintf("Modelo esférico | Range: %.1f km | Sill: %.1f | Nugget: %.2f",
                       range_m / 1000, sill, nugget),
    x = "Distancia (m)",
    y = "Semivarianza γ(h)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14)
  )

path_vario <- file.path(PATHS$output, "figures", "semivariograma.png")
ggsave(path_vario, p_vario, width = 10, height = 7, dpi = 300)
log_msg(sprintf("F4: Semivariograma guardado: %s", path_vario), "OK")
# ============================================================================
# MÓDULO G — EXPORTACIÓN DE RESULTADOS
# ============================================================================
log_msg("MÓDULO G: Exportando resultados del análisis espacial", "INFO")
# --- G1. Guardar muestra con resultados LISA ---------------------------------
path_muestra_lisa <- file.path(carpeta_proc, "muestra_grts_con_lisa.rds")
guardar_rds(muestra_sf, path_muestra_lisa)
# --- G2. Guardar resultados del análisis espacial ----------------------------
resultados_espaciales <- list(
fecha_analisis = Sys.time(),
moran_global_reina = list(
I = moran_test_reina$estimate["Moran I"],
E_I = moran_test_reina$estimate["Expectation"],
p_value = moran_test_reina$p.value,
interpretacion = interpretacion_moran
),
moran_global_knn = list(
I = moran_test_knn$estimate["Moran I"],
E_I = moran_test_knn$estimate["Expectation"],
p_value = moran_test_knn$p.value
),
resumen_clusters = resumen_clusters,
variograma = list(
modelo = modelo_esf,
empirico = variograma_emp
),
resumen_variable = resumen_var,
figuras_generadas = c(path_mapa_densidad, path_mapa_lisa, path_scatter, path_vario)
)
path_resultados <- file.path(carpeta_proc, "estadistica_espacial.rds")
guardar_rds(resultados_espaciales, path_resultados)
# ============================================================================
# MÓDULO H — RESUMEN FINAL
# ============================================================================
tiempo_total <- round(as.numeric(difftime(Sys.time(), inicio, units = "mins")), 2)
message("\n", paste(rep("=", 72), collapse = ""))
message("ANÁLISIS EXPLORATORIO DE DATOS ESPACIALES COMPLETADO")
message(paste(rep("=", 72), collapse = ""))
message("  Tiempo total        : ", tiempo_total, " min")
message("  Celdas analizadas   : ", nrow(muestra_sf))
message("  Moran I (reina)     : ", round(moran_test_reina$estimate["Moran I"], 4))
message("  p-value             : ", format(moran_test_reina$p.value, digits = 4))
message("  Hot spots (HH)      : ", sum(muestra_sf$cluster_lisa == "HH (Hot Spot)"))
message("  Cold spots (LL)     : ", sum(muestra_sf$cluster_lisa == "LL (Cold Spot)"))
message("  Range semivariograma: ", round(modelo_esf$range[2] / 1000, 1), " km")
message(paste(rep("=", 72), collapse = ""))
message("\nHallazgos clave:")
message("  1. ", interpretacion_moran)
message("  2. Hot spots identificados: ", sum(muestra_sf$cluster_lisa == "HH (Hot Spot)"), " celdas")
message("  3. Cold spots identificados: ", sum(muestra_sf$cluster_lisa == "LL (Cold Spot)"), " celdas")
message("  4. Alcance espacial (range): ", round(modelo_esf$range[2] / 1000, 1), " km")
message("\nArchivos exportados:")
message("  - ", path_muestra_lisa)
message("  - ", path_resultados)
message("  - ", path_mapa_densidad)
message("  - ", path_mapa_lisa)
message("  - ", path_scatter)
message("  - ", path_vario)
message("\nSiguiente paso sugerido:")
message("   source('scripts/04_matching.R')  # Propensity Score Matching")
message("   source('scripts/05_ml_models.R') # Machine Learning (Random Forest)")

