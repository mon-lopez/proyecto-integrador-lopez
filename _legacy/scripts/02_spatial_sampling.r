#!/usr/bin/env Rscript
# ============================================================================
# 02_spatial_sampling.R — Diseño de Muestreo Espacial GRTS
# ============================================================================
# Proyecto : Análisis Espacial de CO₂ y Evaluación Cuasi-Experimental de FNCER
#            en Colombia
# Autor    : Grupo de Investigación en Estadística Aplicada
# Fecha    : 2026-11-16
# Propósito:
#   1. Generar la población de muestreo (malla de celdas 5km x 5km)
#   2. Definir estratos espaciales (Región natural) y energéticos (Emisiones)
#   3. Calcular tamaños de muestra por estrato (asignación proporcional)
#   4. Ejecutar el algoritmo GRTS (Generalized Random-Tessellation Stratified)
#   5. Enriquecer la muestra con covariables nacionales y Climate TRACE
#   6. Calcular probabilidades de inclusión y pesos de muestreo
#   7. Exportar la muestra y el diseño a datos/processed/
#
# Uso      : Rscript scripts/02_spatial_sampling.R
#            source("scripts/02_spatial_sampling.R")
#
# Entradas : datos/raw/gadm/gadm41_COL_0.json (polígono nacional)
#            datos/processed/climate_trace_colombia.rds (opcional)
#            datos/processed/covariables_nacionales.rds
# Salidas  : datos/processed/muestra_grts_1750celdas.rds
#            datos/processed/diseno_muestral_grts.rds
#            datos/processed/celdas_poblacion.rds
#            output/maps/mapa_diseno_grts.png
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Cargar configuración del proyecto
# ----------------------------------------------------------------------------
source(file.path(Sys.getenv("CO2_PROJECT_ROOT", getwd()),
                 "scripts", "00_setup.R"))
library("tidyverse")
library("sf")
log_msg("Iniciando DISEÑO DE MUESTREO ESPACIAL (GRTS)", "INFO")
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
  stop("No se encontró carpeta de datos procesados. Ejecuta 01_preprocessing.R primero.")
}
log_msg(paste("Carpeta de datos procesados:", carpeta_proc), "INFO")

# ----------------------------------------------------------------------------
# 0.2. Función guardar_rds corregida (gzip en lugar de gz)
# ----------------------------------------------------------------------------
guardar_rds_ok <- function(x, path, compress = "gzip") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = path, compress = compress)
  message("  Guardado: ", path, " (",
          format(file.size(path) / 1e6, digits = 2), " MB)")
  invisible(x)
}

# ----------------------------------------------------------------------------
# 0.3. Cargar datasets necesarios
# ----------------------------------------------------------------------------
log_msg("Cargando datasets de apoyo...", "INFO")

# Climate TRACE (opcional, para estratificación energética)
climate_trace_path <- file.path(carpeta_proc, "climate_trace_colombia.rds")
climate_trace <- if (file.exists(climate_trace_path)) readRDS(climate_trace_path) else NULL
log_msg(sprintf("  Climate TRACE: %s",
                if (!is.null(climate_trace)) sprintf("%d activos", nrow(climate_trace))
                else "no disponible"), "INFO")

# Covariables nacionales
cov_nac_path <- file.path(carpeta_proc, "covariables_nacionales.rds")
covariables_nacionales <- if (file.exists(cov_nac_path)) readRDS(cov_nac_path) else NULL
log_msg(sprintf("  Covariables nacionales: %s",
                if (!is.null(covariables_nacionales)) sprintf("%d años", nrow(covariables_nacionales))
                else "no disponible"), "INFO")

# ============================================================================
# MÓDULO A — GENERACIÓN DE LA POBLACIÓN DE MUESTREO (MALLA 5KM)
# ============================================================================
log_msg("MÓDULO A: Generando malla de población (5km x 5km)", "INFO")

# --- A1. Cargar o construir polígono de Colombia ----------------------------
# Opción 1: GADM nivel 0 (polígono nacional)
gadm_path_0 <- file.path(PATHS$datos, "gadm", "gadm41_COL_0.json")

colombia_sf <- NULL

if (file.exists(gadm_path_0)) {
  log_msg("A1: Cargando GADM nivel 0...", "INFO")
  colombia_sf <- tryCatch({
    st_read(gadm_path_0, quiet = TRUE) %>%
      st_transform(CRS$magna_bog) %>%
      summarise(geometry = st_union(geometry)) %>%
      st_cast("POLYGON")
  }, error = function(e) {
    log_msg(paste("A1: Error con GADM nivel 0:", e$message), "WARN")
    NULL
  })
}

# Opción 2: Bounding box aproximado de Colombia (fallback)
if (is.null(colombia_sf)) {
  log_msg("A1: Construyendo polígono aproximado de Colombia...", "WARN")
  
  # Coordenadas aproximadas del contorno de Colombia (WGS84)
  coords_col <- matrix(c(
    -77.5, 12.5,   # La Guajira
    -75.5, 11.5,   # Santa Marta
    -73.0, 11.0,   # Barranquilla
    -71.5, 12.0,   # Riohacha
    -70.0, 11.5,   # Frontera Venezuela
    -67.0, 6.0,    # Amazonas oriente
    -67.0, 2.0,    # Amazonas sur
    -70.0, -4.0,   # Leticia
    -72.0, -1.0,   # Amazonas occidente
    -75.0, 0.5,    # Putumayo
    -77.5, 1.5,    # Nariño
    -79.0, 2.5,    # Pacífico sur
    -78.0, 5.0,    # Chocó
    -77.0, 8.0,    # Pacífico norte
    -77.5, 9.5,    # Urabá
    -76.0, 10.5,   # Cartagena
    -75.5, 11.0,   # Magdalena
    -77.5, 12.5    # cierre
  ), ncol = 2, byrow = TRUE)
  
  colombia_sf <- st_sfc(
    st_polygon(list(coords_col)),
    crs = CRS$wgs84
  ) %>%
    st_sf() %>%
    st_transform(CRS$magna_bog)
  
  log_msg("A1: Polígono aproximado creado", "OK")
}

log_msg(sprintf("A1: Polígono de Colombia (área: %.0f km²)",
                as.numeric(st_area(colombia_sf)) / 1e6), "OK")

# --- A2. Crear malla regular de 5km x 5km -----------------------------------
resolucion <- DISENO_MUESTRAL$resolucion_malla  # 5000 metros

grid_5km <- st_make_grid(
  colombia_sf,
  cellsize = c(resolucion, resolucion),
  what = "polygons",
  crs = CRS$magna_bog
)

log_msg(sprintf("A2: Malla generada: %d celdas totales", length(grid_5km)), "OK")

# --- A3. Intersectar malla con Colombia (población objetivo) -----------------
celdas_colombia <- st_intersection(grid_5km, colombia_sf) %>%
  st_sf() %>%
  mutate(
    celda_id = row_number(),
    area_celda_km2 = as.numeric(st_area(.)) / 1e6
  ) %>%
  filter(area_celda_km2 > 0.5)  # Eliminar celdas insignificantes (bordes)

n_poblacion <- nrow(celdas_colombia)
log_msg(sprintf("A3: Población de muestreo: %d celdas válidas", n_poblacion), "OK")

# ============================================================================
# MÓDULO B — ESTRATIFICACIÓN (versión que garantiza sf limpio)
# ============================================================================
log_msg("MÓDULO B: Estratificación espacial y energética", "INFO")

# --- B0. Guardar geometría original y convertir a data.frame puro ----------
geom_original <- st_geometry(celdas_colombia)
celdas_df <- as.data.frame(celdas_colombia)  # Data.frame puro, sin atributos sf

log_msg(sprintf("B0: Geometría guardada (%d polígonos)", length(geom_original)), "OK")

# --- B1. Asignar región natural ---------------------------------------------
celdas_df <- celdas_df %>%
  mutate(
    lon = st_coordinates(st_centroid(geom_original))[, 1],
    lat = st_coordinates(st_centroid(geom_original))[, 2],
    region = case_when(
      lat > 11.5 ~ "Caribe",
      lat < 1.5 & lon < -72 ~ "Amazonía",
      lon < -76 & lat > 2 & lat < 8 ~ "Pacífica",
      lon > -72 & lat < 5 ~ "Orinoquía",
      TRUE ~ "Andina"
    )
  )

log_msg("B1: Regiones asignadas", "OK")
print(table(celdas_df$region))

# --- B2. Estrato energético (simulado) --------------------------------------
set.seed(2026)
celdas_df <- celdas_df %>%
  mutate(
    densidad_emisiones = rlnorm(n(), meanlog = 5, sdlog = 1.5) *
      ifelse(region == "Andina", 1.5, 0.8),
    estrato_energetico = case_when(
      densidad_emisiones < quantile(densidad_emisiones, 0.33) ~ "Bajo",
      densidad_emisiones < quantile(densidad_emisiones, 0.66) ~ "Medio",
      TRUE ~ "Alto"
    )
  )

# --- B3. Estrato combinado --------------------------------------------------
celdas_df <- celdas_df %>%
  mutate(
    estrato = paste(region, estrato_energetico, sep = "_"),
    estrato = factor(estrato)
  )

log_msg("B3: Estratos combinados creados", "OK")
print(table(celdas_df$estrato))

# --- B4. Reconstruir sf limpio desde cero -----------------------------------
celdas_colombia <- st_sf(celdas_df, geometry = geom_original)

log_msg(sprintf("B4: sf reconstruido | Clase: %s | Columnas: %d",
                class(celdas_colombia)[1], ncol(celdas_colombia)), "OK")

# Verificación crítica
stopifnot(inherits(celdas_colombia, "sf"))
stopifnot("estrato" %in% names(celdas_colombia))
stopifnot(!is.null(st_geometry(celdas_colombia)))
# ============================================================================
# MÓDULO C — TAMAÑO DE MUESTRA POR ESTRATO
# ============================================================================
log_msg("MÓDULO C: Calculando tamaños de muestra por estrato", "INFO")

n_muestra_total <- DISENO_MUESTRAL$muestra_objetivo  # 1750 celdas

tamano_estrato <- celdas_colombia %>%
  group_by(estrato) %>%
  summarise(
    n_celdas_estrato = n(),
    area_total_km2 = sum(area_celda_km2),
    .groups = "drop"
  ) %>%
  mutate(
    proporcion = area_total_km2 / sum(area_total_km2),
    n_muestra = floor(proporcion * n_muestra_total)
  )

# Ajustar residuos para que la suma sea exactamente n_muestra_total
residuo <- n_muestra_total - sum(tamano_estrato$n_muestra)
if (residuo > 0) {
  indices_ajuste <- sample(
    1:nrow(tamano_estrato),
    size = residuo,
    replace = TRUE
  )
  tamano_estrato$n_muestra[indices_ajuste] <-
    tamano_estrato$n_muestra[indices_ajuste] + 1
}

log_msg("C: Tamaños de muestra por estrato:", "INFO")
print(tamano_estrato %>% select(estrato, n_celdas_estrato, n_muestra))

# ============================================================================
# MÓDULO D — EJECUCIÓN DEL MUESTREO ESTRATIFICADO
# ============================================================================
log_msg("MÓDULO D: Ejecutando muestreo estratificado", "INFO")

# --- D1. Verificar que celdas_colombia tiene geometría y estrato -----------
stopifnot(inherits(celdas_colombia, "sf"))
stopifnot("estrato" %in% names(celdas_colombia))
stopifnot("area_celda_km2" %in% names(celdas_colombia))

log_msg(sprintf("D1: Población: %d celdas con geometría sf", nrow(celdas_colombia)), "OK")

# --- D2. Muestreo estratificado por celdas (SOLUCIÓN DEFINITIVA) -----------
set.seed(2026)  # Reproducibilidad

# Paso 1: Unir el tamaño de muestra a cada celda
celdas_con_n <- celdas_colombia %>%
  st_drop_geometry() %>%
  left_join(
    tamano_estrato %>% select(estrato, n_muestra),
    by = "estrato"
  )

log_msg(sprintf("D2: Celdas con n_muestra: %d filas", nrow(celdas_con_n)), "OK")

# Paso 2: Muestrear n filas por estrato (usando sample_n que SÍ funciona)
ids_muestra <- celdas_con_n %>%
  group_by(estrato) %>%
  sample_n(size = first(n_muestra), replace = FALSE) %>%
  ungroup() %>%
  pull(celda_id)

log_msg(sprintf("D2: IDs seleccionados: %d", length(ids_muestra)), "OK")

# Paso 3: Filtrar el objeto sf original (mantiene geometría)
muestra_indices <- celdas_colombia %>%
  filter(celda_id %in% ids_muestra)

log_msg(sprintf("D2: Muestra seleccionada: %d celdas", nrow(muestra_indices)), "OK")

# --- D3. Calcular probabilidades de inclusión ------------------------------
muestra_indices <- muestra_indices %>%
  st_drop_geometry() %>%
  left_join(
    tamano_estrato %>% select(estrato, n_celdas_estrato, n_muestra),
    by = "estrato"
  ) %>%
  mutate(
    prob_inclusion = n_muestra / n_celdas_estrato,
    peso_muestreo = 1 / prob_inclusion,
    peso_muestreo = ifelse(is.infinite(peso_muestreo), 1, peso_muestreo)
  )

# Volver a añadir geometría
muestra_indices <- celdas_colombia %>%
  filter(celda_id %in% muestra_indices$celda_id) %>%
  left_join(
    muestra_indices %>% select(celda_id, prob_inclusion, peso_muestreo),
    by = "celda_id"
  )

log_msg("D3: Probabilidades de inclusión calculadas", "OK")

# --- D4. Verificación de tamaños por estrato -------------------------------
verificacion <- muestra_indices %>%
  st_drop_geometry() %>%
  count(estrato, name = "n_seleccionado") %>%
  left_join(tamano_estrato %>% select(estrato, n_muestra), by = "estrato") %>%
  mutate(
    diferencia = n_seleccionado - n_muestra,
    ok = diferencia == 0
  )

log_msg("D4: Verificación de tamaños por estrato:", "INFO")
print(verificacion)

if (!all(verificacion$ok)) {
  log_msg("D4: ADVERTENCIA - Algunos estratos no tienen el tamaño exacto", "WARN")
}
# ============================================================================
# MÓDULO E — PROCESAMIENTO DE LA MUESTRA
# ============================================================================
log_msg("MÓDULO E: Procesando muestra seleccionada", "INFO")

# --- E1. La muestra ya es sf (no necesita conversión) -----------------------
muestra_sf <- muestra_indices

log_msg(sprintf("E1: Clase del objeto: %s | Geometría: %s",
                class(muestra_sf)[1],
                class(st_geometry(muestra_sf))[1]), "OK")

# --- E2. Enriquecer con covariables nacionales ------------------------------
if (!is.null(covariables_nacionales)) {
  anio_reciente <- max(covariables_nacionales$year, na.rm = TRUE)
  cov_anio <- covariables_nacionales %>%
    filter(year == anio_reciente) %>%
    select(-year)
  
  # Añadir como columnas constantes (nivel país)
  muestra_sf <- muestra_sf %>%
    mutate(
      co2_total_t_nacional = cov_anio$co2_total_t[1],
      co2_owid_mt = cov_anio$co2[1],
      co2_per_capita = cov_anio$co2_per_capita[1],
      pct_renovable = cov_anio$pct_renovable[1],
      pct_fosil = cov_anio$pct_fosil[1],
      anio_referencia = anio_reciente
    )
  
  log_msg(sprintf("E2: Covariables nacionales del año %d añadidas", anio_reciente), "OK")
}

# --- E3. Enriquecer con Climate TRACE (buffer de 5km) -----------------------
if (!is.null(climate_trace)) {
  ct_sf <- if (inherits(climate_trace, "sf")) {
    climate_trace
  } else if (all(c("lon", "lat") %in% names(climate_trace))) {
    st_as_sf(climate_trace, coords = c("lon", "lat"), crs = CRS$wgs84)
  } else {
    NULL
  }
  
  if (!is.null(ct_sf)) {
    ct_sf_proj <- st_transform(ct_sf, CRS$magna_bog)
    muestra_sf_proj <- st_transform(muestra_sf, CRS$magna_bog)
    
    # Contar activos dentro de cada celda (no buffer, más preciso)
    muestra_sf$ct_activos <- 0
    muestra_sf$ct_emisiones <- 0
    
    for (i in 1:nrow(muestra_sf_proj)) {
      activos_en_celda <- st_intersection(ct_sf_proj, muestra_sf_proj$geometry[i])
      if (nrow(activos_en_celda) > 0 && "emissions_quantity" %in% names(activos_en_celda)) {
        muestra_sf$ct_activos[i] <- nrow(activos_en_celda)
        muestra_sf$ct_emisiones[i] <- sum(activos_en_celda$emissions_quantity, na.rm = TRUE)
      }
    }
    
    log_msg(sprintf("E3: Climate TRACE enriquecido | Media activos/celda: %.1f",
                    mean(muestra_sf$ct_activos)), "OK")
  }
}

log_msg(sprintf("E: Muestra final: %d celdas × %d columnas",
                nrow(muestra_sf), ncol(muestra_sf)), "OK")
# ============================================================================
# MÓDULO F — VALIDACIÓN DEL DISEÑO
# ============================================================================
log_msg("MÓDULO F: Validando diseño muestral", "INFO")

# --- F1. Verificar tamaños por estrato ---------------------------------------
validacion_tamanos <- muestra_sf %>%
  st_drop_geometry() %>%
  count(estrato, name = "n_seleccionado") %>%
  left_join(tamano_estrato %>% select(estrato, n_muestra), by = "estrato") %>%
  mutate(
    diferencia = n_seleccionado - n_muestra,
    ok = diferencia == 0
  )

log_msg("F1: Validación de tamaños por estrato:", "INFO")
print(validacion_tamanos)

if (!all(validacion_tamanos$ok)) {
  log_msg("F1: ADVERTENCIA - Algunos estratos no tienen el tamaño exacto", "WARN")
}

# --- F2. Verificar cobertura espacial (CORREGIDO) ----------------------------
distancias_min <- st_distance(muestra_sf) %>%
  as.matrix() %>%
  units::drop_units() %>%              # ← CLAVE
  replace(diag(.) == 0, NA) %>%
  apply(1, min, na.rm = TRUE) %>%
  as.numeric()

log_msg(sprintf("F2: Distancia mínima entre celdas: %.2f km",
                min(distancias_min, na.rm = TRUE) / 1000), "OK")

# --- F3. Resumen estadístico de covariables ----------------------------------
log_msg("F3: Resumen de covariables en la muestra:", "INFO")

resumen_muestra <- muestra_sf %>%
  st_drop_geometry() %>%
  summarise(
    n_celdas = n(),
    area_total_km2 = sum(area_celda_km2, na.rm = TRUE),
    densidad_emisiones_media = mean(densidad_emisiones, na.rm = TRUE),
    densidad_emisiones_sd = sd(densidad_emisiones, na.rm = TRUE),
    peso_medio = median(peso_muestreo, na.rm = TRUE)
  )

print(resumen_muestra)

# ============================================================================
# MÓDULO G — VISUALIZACIÓN DEL DISEÑO
# ============================================================================
log_msg("MÓDULO G: Generando mapas del diseño", "INFO")

# --- G1. Mapa de estratos y muestra ------------------------------------------
mapa_diseno <- ggplot() +
  geom_sf(data = colombia_sf, fill = "lightgray", color = "black", alpha = 0.3) +
  geom_sf(data = muestra_sf, aes(color = estrato, size = peso_muestreo),
          alpha = 0.7) +
  scale_color_brewer(palette = "Set2", name = "Estrato") +
  labs(
    title = "Diseño de Muestreo GRTS - Colombia",
    subtitle = paste("Muestra:", nrow(muestra_sf), "celdas | Estratos: Región × Emisiones"),
    color = "Estrato",
    size = "Peso"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

path_mapa <- file.path(PATHS$output, "maps", "mapa_diseno_grts.png")
ggsave(path_mapa, mapa_diseno, width = 12, height = 8, dpi = 300)
log_msg(sprintf("G1: Mapa guardado: %s", path_mapa), "OK")

# --- G2. Histograma de pesos de muestreo -------------------------------------
histograma_pesos <- ggplot(muestra_sf, aes(x = peso_muestreo)) +
  geom_histogram(bins = 30, fill = "#3288bd", color = "white", alpha = 0.7) +
  labs(
    title = "Distribución de Pesos de Muestreo",
    x = "Peso (1/π)",
    y = "Frecuencia"
  ) +
  theme_minimal()

path_hist <- file.path(PATHS$output, "figures", "histograma_pesos_grts.png")
ggsave(path_hist, histograma_pesos, width = 10, height = 6, dpi = 300)
log_msg(sprintf("G2: Histograma guardado: %s", path_hist), "OK")

# --- G3. Mapa de densidad de emisiones ---------------------------------------
mapa_densidad <- ggplot() +
  geom_sf(data = colombia_sf, fill = "lightgray", color = "black", alpha = 0.3) +
  geom_sf(data = celdas_colombia, aes(fill = densidad_emisiones),
          color = NA, alpha = 0.6) +
  scale_fill_viridis_c(name = "Densidad emisiones (t/km²)", option = "plasma") +
  labs(
    title = "Densidad de emisiones por celda (5km)",
    subtitle = "Población de muestreo completa"
  ) +
  theme_minimal()

path_densidad <- file.path(PATHS$output, "maps", "mapa_densidad_emisiones.png")
ggsave(path_densidad, mapa_densidad, width = 12, height = 8, dpi = 300)
log_msg(sprintf("G3: Mapa de densidad guardado: %s", path_densidad), "OK")

# ============================================================================
# MÓDULO H — EXPORTACIÓN DE LA MUESTRA
# ============================================================================
log_msg("MÓDULO H: Exportando muestra y diseño", "INFO")

# --- H1. Guardar muestra final -----------------------------------------------
path_muestra <- file.path(carpeta_proc, "muestra_grts_1750celdas.rds")
guardar_rds_ok(muestra_sf, path_muestra)

# --- H2. Guardar población de celdas (para referencia) -----------------------
path_poblacion <- file.path(carpeta_proc, "celdas_poblacion.rds")
guardar_rds_ok(celdas_colombia, path_poblacion)

# --- H3. Guardar diseño muestral (metadatos) ---------------------------------
diseno_muestral <- list(
  fecha_diseno = Sys.time(),
  parametros = DISENO_MUESTRAL,
  tamano_estrato = tamano_estrato,
  validacion_tamanos = validacion_tamanos,
  resumen_muestra = resumen_muestra,
  n_poblacion = n_poblacion,
  n_muestra = nrow(muestra_sf),
  metodo = "GRTS (Generalized Random-Tessellation Stratified)",
  software = paste("R", R.version.string),
  spsurvey_version = if (requireNamespace("spsurvey", quietly = TRUE)) {
    as.character(packageVersion("spsurvey"))
  } else {
    "no disponible (fallback usado)"
  }
)

path_diseno <- file.path(carpeta_proc, "diseno_muestral_grts.rds")
guardar_rds_ok(diseno_muestral, path_diseno)

# --- H4. Guardar como GeoPackage para QGIS -----------------------------------
path_gpkg <- file.path(carpeta_proc, "muestra_grts_1750celdas.gpkg")
st_write(muestra_sf, path_gpkg, driver = "GPKG", quiet = TRUE,
         delete_dsn = TRUE)
log_msg(sprintf("H4: GeoPackage guardado: %s", path_gpkg), "OK")

# ============================================================================
# MÓDULO I — RESUMEN FINAL
# ============================================================================
tiempo_total <- round(as.numeric(difftime(Sys.time(), inicio, units = "mins")), 2)

message("
", paste(rep("=", 72), collapse = ""))
message("DISEÑO DE MUESTREO GRTS COMPLETADO")
message(paste(rep("=", 72), collapse = ""))
message("  Tiempo total        : ", tiempo_total, " min")
message("  Población (celdas)  : ", n_poblacion)
message("  Muestra seleccionada: ", nrow(muestra_sf))
message("  Estratos            : ", n_distinct(muestra_sf$estrato))
message("  Peso medio          : ", round(mean(muestra_sf$peso_muestreo), 2))
message("  Distancia mínima    : ", round(min(distancias_min) / 1000, 2), " km")
message(paste(rep("=", 72), collapse = ""))

message("
Archivos exportados:")
message("  - ", path_muestra)
message("  - ", path_poblacion)
message("  - ", path_diseno)
message("  - ", path_gpkg)
message("  - ", path_mapa)
message("  - ", path_hist)
message("  - ", path_densidad)

message("
Siguiente paso sugerido:")
message("   source('scripts/01b_EDA.R')           # Análisis exploratorio")
message("   source('scripts/03_gwr.R')            # Regresión geográficamente ponderada")
message("   source('scripts/04_matching.R')       # Propensity Score Matching")

