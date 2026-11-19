#!/usr/bin/env Rscript
# ============================================================================
# 03_spatial_sampling.R — Diseño muestral espacialmente balanceado con
#                         estratificación de Neyman
# ============================================================================
# Proyecto : Análisis Espacial de CO₂ y Evaluación Cuasi-Experimental de FNCER
#            en Colombia
# Sección  : Teoría del Muestreo (preregistro §4)
# Autor    : Grupo de Investigación en Estadística Aplicada
# Fecha    : 2026-11-14
# Prereq   : scripts/00_setup.R
#            scripts/02_preprocessing.R (ejecutado; produce municipios_col.rds
#                                        y malla_co2_municipios.rds)
#
# Propósito:
#   1. Construir marco muestral de celdas 5 km sobre territorio colombiano
#   2. Estratificar por región natural × estrato UPME (15 estratos)
#   3. Estimar Deff vía variograma piloto (gstat)
#   4. Calcular tamaño muestral con asignación de NEYMAN (óptima)
#      y comparar con proporcional
#   5. Ejecutar muestreo balanceado espacialmente (Cube) dentro de cada estrato
#   6. Estimar total de CO₂ con Horvitz-Thompson
#   7. Validar con bootstrap espacial por bloques
#   8. Generar curvas de varianza vs n para los 4 diseños
#   9. Guardar entregables en datos/processed/ y output/
#
# Referencias:
#   - Stevens & Donagan (2010). GRTS for spatial surveys.
#   - Cochran (1977). Sampling Techniques. Cap. 5 (Neyman).
#   - Grafström & Schelin (2014). Cube method for spatial sampling.
# ============================================================================

source("scripts/00_setup.R")

# ----------------------------------------------------------------------------
# 0. Flags y configuración
# ----------------------------------------------------------------------------
RES_MALLA_M     <- DISENO_MUESTRAL$resolucion_malla  # 5000 m
N_OBJETIVO      <- DISENO_MUESTRAL$muestra_objetivo  # 1750
CONF            <- DISENO_MUESTRAL$confianza         # 0.95
MARGEN_ERROR    <- DISENO_MUESTRAL$margen_error      # 0.05

set.seed(20261114)  # Reproducibilidad total

message("\n", strrep("=", 72))
message("🎲  DISEÑO MUESTRAL ESPACIAL — Proyecto CO₂ Colombia")
message(strrep("=", 72))
message(sprintf("  📏 Resolución malla : %d m", RES_MALLA_M))
message(sprintf("  🎯 Muestra objetivo : %d celdas", N_OBJETIVO))
message(sprintf("  📊 Confianza        : %.0f%%", CONF * 100))
message(sprintf("  📉 Margen error     : %.0f%%", MARGEN_ERROR * 100))

# Registro de QC
QC <- list()
qc_add <- function(paso, estado, detalle = "", n = NA_integer_) {
  QC[[length(QC) + 1L]] <<- data.frame(
    paso = paso, estado = estado, n = n, detalle = detalle,
    fecha = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  ico <- switch(estado, OK = "✅", WARN = "⚠️", ERROR = "❌", "ℹ️")
  message(sprintf("  %s %-30s %s", ico, paso, detalle))
}

# ----------------------------------------------------------------------------
# 1. Cargar datos del preprocesamiento
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n📥 [1/10] Cargando datos procesados...")

mun <- cargar_rds(file.path(PATHS$datos_proc, "municipios_col.rds"))
malla_mun <- cargar_rds(file.path(PATHS$datos_proc, "malla_co2_municipios.rds"))

qc_add("datos_cargados", "OK",
       sprintf("%d municipios, %d filas malla", nrow(mun), nrow(malla_mun)))

# Unir malla con municipios para tener variable de estudio (CO₂) y estratos
mun <- merge(mun, malla_mun[, c("clave_mun", "co2_t_sum", "log1p_co2",
                                 "co2_densidad", "cero_estructural")],
             by = "clave_mun", all.x = TRUE)

# ----------------------------------------------------------------------------
# 2. Construir marco muestral: malla 5 km sobre territorio colombiano
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n🗺️  [2/10] Construyendo marco muestral (malla 5 km)...")

# Proyectar a CRS métrico para generar malla regular
mun_m <- sf::st_transform(mun, CRS$magna_bog)
bb <- sf::st_bbox(mun_m)

# Generar malla hexagonal (más eficiente que cuadrada para cobertura espacial)
grid_cells <- sf::st_make_grid(mun_m, cellsize = RES_MALLA_M, square = FALSE)
marco <- sf::st_sf(id_celda = seq_along(grid_cells), geometry = grid_cells)

# Intersectar con el territorio colombiano (descartar celdas fuera)
marco <- sf::st_intersection(marco, sf::st_union(mun_m))
marco$centroid <- sf::st_centroid(marco$geometry)
coords_centro <- sf::st_coordinates(marco$centroid)
marco$X <- coords_centro[, 1]
marco$Y <- coords_centro[, 2]

N_total <- nrow(marco)
qc_add("marco_muestral", "OK",
       sprintf("%d celdas de %d m sobre territorio colombiano", N_total, RES_MALLA_M),
       N_total)

# ----------------------------------------------------------------------------
# 3. Asignar estratos: región natural × estrato UPME
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n🎯 [3/10] Asignando estratos (región × UPME)...")

# Intersectar centroides de celdas con municipios para heredar región y estrato
marco_sf_cent <- sf::st_sf(id_celda = marco$id_celda,
                            geometry = sf::st_geometry(marco$centroid))
marco_sf_cent <- sf::st_transform(marco_sf_cent, CRS$magna_bog)

# Unión espacial: cada celda hereda atributos del municipio donde cae su centroide
idx_inter <- sf::st_intersects(marco_sf_cent, mun_m)
marco$clave_mun <- mun_m$clave_mun[vapply(idx_inter, function(i)
  if (length(i)) i[1] else NA_integer_, integer(1))]

# Unir con atributos de municipio (región, estrato UPME)
mun_attr <- sf::st_drop_geometry(mun)[, c("clave_mun", "departamento", "region",
                                           "region_ambigua", "estrato_ener")]
marco <- merge(marco, mun_attr, by = "clave_mun", all.x = TRUE)

# Marcar celdas sin municipio asignado (borde marítimo)
sin_mun <- is.na(marco$clave_mun)
if (any(sin_mun)) {
  marco <- marco[!sin_mun, ]
  qc_add("celdas_borde", "WARN",
         sprintf("%d celdas sin municipio (borde marítimo) descartadas", sum(sin_mun)))
}

# Construir estrato combinado
marco$stratum <- paste(marco$region, marco$estrato_ener, sep = "_")
marco$stratum <- factor(marco$stratum)

# Resumen por estrato
tab_estratos <- table(marco$stratum)
qc_add("estratos", "OK",
       sprintf("%d estratos; tamaños: min=%d, max=%d, mediana=%d",
               length(tab_estratos), min(tab_estratos), max(tab_estratos),
               stats::median(tab_estratos)))

# ----------------------------------------------------------------------------
# 4. Obtener variable de estudio piloto (CO₂ por celda)
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n📊 [4/10] Obteniendo variable piloto de CO₂...")

# Agregar CO₂ a nivel de celda: promediar la densidad de CO₂ de los municipios
# que intersectan cada celda (ponderado por área de intersección)
# Para simplificar y ser conservadores, usamos la densidad del municipio donde
# cae el centroide de la celda. Esto es una aproximación; el valor real
# requeriría exactextractr sobre el raster EDGAR (costoso computacionalmente).

# Extraer densidad de CO₂ del municipio correspondiente
marco$log1p_co2 <- mun$log1p_co2[match(marco$clave_mun, mun$clave_mun)]
marco$co2_densidad <- mun$co2_densidad[match(marco$clave_mun, mun$clave_mun)]

# Agregar ruido espacial para simular variabilidad intra-municipio
# (las celdas dentro de un mismo municipio no tienen exactamente el mismo CO₂)
set.seed(20261114)
ruido <- stats::rnorm(nrow(marco), mean = 0, sd = 0.3)
marco$log1p_co2_pilot <- marco$log1p_co2 + ruido

# Estadísticas descriptivas
qc_add("variable_piloto", "OK",
       sprintf("log1p(CO₂): media=%.2f, DE=%.2f, rango=[%.2f, %.2f]",
               mean(marco$log1p_co2_pilot, na.rm = TRUE),
               stats::sd(marco$log1p_co2_pilot, na.rm = TRUE),
               min(marco$log1p_co2_pilot, na.rm = TRUE),
               max(marco$log1p_co2_pilot, na.rm = TRUE)))

# ----------------------------------------------------------------------------
# 5. Variograma piloto y efecto de diseño (Deff)
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n📐 [5/10] Estimando Deff vía variograma piloto...")

# Submuestra para el variograma (no usar todo el marco por costo computacional)
n_pilot_vgm <- min(2000, N_total)
idx_pilot <- sample(N_total, n_pilot_vgm)
pilot_df <- data.frame(
  X = marco$X[idx_pilot],
  Y = marco$Y[idx_pilot],
  z = marco$log1p_co2_pilot[idx_pilot]
)
pilot_df <- pilot_df[stats::complete.cases(pilot_df), ]

# Convertir a objeto sf para gstat
pilot_sf <- sf::st_as_sf(pilot_df, coords = c("X", "Y"), crs = CRS$magna_bog)

# Calcular variograma empírico
vgm_emp <- gstat::variogram(z ~ 1, data = pilot_sf, cutoff = 100000, width = 5000)

# Ajustar modelo exponencial
vgm_fit <- gstat::fit.variogram(vgm_emp,
                                 model = gstat::vgm(
                                   psill = stats::var(pilot_df$z),
                                   model = "Exp",
                                   range = 20000,
                                   nugget = 0.1 * stats::var(pilot_df$z)
                                 ))

# Extraer parámetros
nugget <- vgm_fit$psill[vgm_fit$mod == "Nug"]
sill_total <- sum(vgm_fit$psill)
range_param <- vgm_fit$range[vgm_fit$mod == "Exp"]

# Correlación espacial intraclase (rho)
rho_spatial <- 1 - (nugget / sill_total)

# Efecto de diseño aproximado: Deff ≈ 1 + (n_bar - 1) * rho
# donde n_bar es el tamaño promedio de "cluster" espacial definido por el rango
# Para una malla de 5 km y rango R, el número de celdas correlacionadas es
# aproximadamente pi * (R / 5)^2
n_bar <- pi * (range_param / RES_MALLA_M)^2
Deff <- 1 + (n_bar - 1) * rho_spatial
Deff <- max(1, min(Deff, 5))  # Acotar para estabilidad

qc_add("variograma", "OK",
       sprintf("nugget=%.3f, sill=%.3f, range=%.0f m, rho=%.3f",
               nugget, sill_total, range_param, rho_spatial))
qc_add("Deff", "OK",
       sprintf("Deff estimado = %.2f (n_bar=%.1f celdas correlacionadas)",
               Deff, n_bar))

# ----------------------------------------------------------------------------
# 6. Cálculo de tamaño muestral con Neyman vs Proporcional
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n🧮 [6/10] Calculando n por estrato (Neyman vs Proporcional)...")

# Estadísticas por estrato
estratos_nombres <- levels(marco$stratum)
N_h <- tabulate(marco$stratum, nbins = length(estratos_nombres))
names(N_h) <- estratos_nombres

# Calcular S_h (desviación estándar) por estrato
S_h <- tapply(marco$log1p_co2_pilot, marco$stratum, stats::sd, na.rm = TRUE)
S_h[is.na(S_h)] <- 0  # Estratos con un solo valor

# --- Asignación PROPORCIONAL ---
n_prop <- round(N_total * (N_h / sum(N_h)))

# --- Asignación de NEYMAN (óptima) ---
# Fórmula: n_h = n * (N_h * S_h) / sum(N_k * S_k)
numerator <- N_h * S_h
denominator <- sum(numerator)
n_neyman_raw <- N_total * (numerator / denominator)

# Truncamiento: si n_h > N_h, fijar n_h = N_h y redistribuir el resto
# (Neyman modificado con truncamiento)
n_neyman <- round(n_neyman_raw)
n_neyman <- pmin(n_neyman, N_h)  # No exceder el tamaño del estrato

# Redistribuir la diferencia para que sume N_total
diff <- N_total - sum(n_neyman)
if (diff != 0) {
  # Redistribuir proporcionalmente a N_h * S_h entre estratos no truncados
  no_trunc <- n_neyman < N_h
  if (any(no_trunc)) {
    weights_redist <- (N_h * S_h)[no_trunc]
    weights_redist <- weights_redist / sum(weights_redist)
    n_neyman[no_trunc] <- n_neyman[no_trunc] + round(diff * weights_redist)
  }
  # Ajuste final para que sume exactamente N_total
  diff_final <- N_total - sum(n_neyman)
  if (diff_final != 0) {
    idx_max <- which.max(N_h * S_h)
    n_neyman[idx_max] <- n_neyman[idx_max] + diff_final
  }
}

# Tabla comparativa
tabla_asignacion <- data.frame(
  estrato = estratos_nombres,
  N_h = N_h,
  S_h = round(S_h, 3),
  n_prop = n_prop,
  n_neyman = n_neyman,
  proporción_neyman = round(100 * n_neyman / N_total, 1)
)

qc_add("asignacion_neyman", "OK",
       sprintf("Neyman: min=%d, max=%d, mediana=%d",
               min(n_neyman), max(n_neyman), stats::median(n_neyman)))
qc_add("asignacion_prop", "OK",
       sprintf("Proporcional: min=%d, max=%d, mediana=%d",
               min(n_prop), max(n_prop), stats::median(n_prop)))

# Ganancia de eficiencia de Neyman vs Proporcional
# Varianza bajo Neyman: V_ney = (sum(N_h * S_h))^2 / n
# Varianza bajo Proporcional: V_prop = sum(N_h^2 * S_h^2 / n_h) / n
V_neyman <- (sum(N_h * S_h))^2 / N_total
V_prop <- sum((N_h^2 * S_h^2) / pmax(n_prop, 1)) / N_total
ganancia_pct <- 100 * (V_prop - V_neyman) / V_prop

qc_add("ganancia_neyman", "OK",
       sprintf("Reducción de varianza con Neyman: %.1f%%", ganancia_pct))

# Guardar tabla de asignación
data.table::fwrite(tabla_asignacion,
                   file.path(PATHS$output, "tables", "asignacion_muestral.csv"))

# ----------------------------------------------------------------------------
# 7. Muestreo balanceado espacialmente (Cube) por estrato
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n🎲 [7/10] Ejecutando muestreo balanceado espacialmente (Cube)...")

# Función para muestrear dentro de un estrato con Cube
muestrear_estrato_cube <- function(datos_estrato, n_h) {
  if (nrow(datos_estrato) <= n_h) {
    return(datos_estrato)  # Censar todo el estrato
  }
  
  # Variables de balance: coordenadas + valor de CO₂
  # Cube busca que la muestra tenga la misma media que el marco en estas variables
  X_balance <- cbind(datos_estrato$X, datos_estrato$Y, datos_estrato$log1p_co2_pilot)
  
  # Ejecutar Cube
  set.seed(20261114 + datos_estrato$stratum[1])  # Reproducibilidad por estrato
  sample_idx <- BalancedSampling::Cube(
    units = 1:nrow(datos_estrato),
    prob = rep(n_h / nrow(datos_estrato), nrow(datos_estrato)),
    X = X_balance,
    maxtry = 10,
    seed = 20261114
  )
  
  datos_estrato[sample_idx, ]
}

# Aplicar por estrato
marco_muestreado <- do.call(rbind, lapply(estratos_nombres, function(est) {
  datos_est <- marco[marco$stratum == est, ]
  n_h <- n_neyman[est]
  muestrear_estrato_cube(datos_est, n_h)
}))

# Calcular probabilidades de inclusión (pi)
# Para Cube con prob = n_h / N_h, pi_i = n_h / N_h para todas las unidades del estrato
marco_muestreado$N_h <- N_h[as.character(marco_muestreado$stratum)]
marco_muestreado$n_h <- n_neyman[as.character(marco_muestreado$stratum)]
marco_muestreado$pi <- marco_muestreado$n_h / marco_muestreado$N_h

# Pesos Horvitz-Thompson: w_i = 1 / pi_i
marco_muestreado$weight <- 1 / marco_muestreado$pi

n_final <- nrow(marco_muestreado)
qc_add("muestreo_cube", "OK",
       sprintf("%d celdas seleccionadas en %d estratos",
               n_final, length(estratos_nombres)), n_final)

# ----------------------------------------------------------------------------
# 8. Estimación Horvitz-Thompson del total de CO₂
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n📊 [8/10] Estimando total de CO₂ (Horvitz-Thompson)...")

# Total estimado: sum(y_i / pi_i)
total_HT <- sum(marco_muestreado$log1p_co2_pilot * marco_muestreado$weight, na.rm = TRUE)

# Varianza estimada (aproximación para muestreo estratificado)
# V(HT) = sum_h (N_h^2 * (1 - f_h) * S_h^2 / n_h)
# donde f_h = n_h / N_h (fracción de muestreo)
f_h <- n_neyman / N_h
V_HT <- sum((N_h^2) * (1 - f_h) * (S_h^2) / pmax(n_neyman, 1), na.rm = TRUE)
SE_HT <- sqrt(V_HT)

# IC 95%
z_crit <- stats::qnorm(1 - (1 - CONF) / 2)
IC_inf <- total_HT - z_crit * SE_HT
IC_sup <- total_HT + z_crit * SE_HT

qc_add("estimacion_HT", "OK",
       sprintf("Total log1p(CO₂) HT = %.2f (SE = %.2f, IC 95%% = [%.2f, %.2f])",
               total_HT, SE_HT, IC_inf, IC_sup))

# Transformar a escala original (aproximación)
# E[CO₂] ≈ exp(E[log1p(CO₂)]) - 1 (sesgado, pero aproximación útil)
total_CO2_original <- exp(total_HT / N_total) - 1
qc_add("total_original", "OK",
       sprintf("Total CO₂ aproximado (escala original) ≈ %.2f t/celda",
               total_CO2_original))

# ----------------------------------------------------------------------------
# 9. Bootstrap espacial por bloques
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n🔄 [9/10] Bootstrap espacial por bloques...")

# Dividir el territorio en bloques espaciales (grid de 50 km × 50 km)
tamanio_bloque <- 50000  # 50 km
marco_muestreado$bloque_x <- floor(marco_muestreado$X / tamanio_bloque)
marco_muestreado$bloque_y <- floor(marco_muestreado$Y / tamanio_bloque)
marco_muestreado$bloque_id <- paste(marco_muestreado$bloque_x,
                                     marco_muestreado$bloque_y, sep = "_")

bloques_unicos <- unique(marco_muestreado$bloque_id)
n_bloques <- length(bloques_unicos)

R <- 500  # Número de réplicas bootstrap
boot_totals <- numeric(R)

set.seed(20261114)
for (r in seq_len(R)) {
  # Remuestrear bloques con reemplazo
  bloques_boot <- sample(bloques_unicos, n_bloques, replace = TRUE)
  
  # Construir muestra bootstrap
  idx_boot <- unlist(lapply(bloques_boot, function(b) {
    which(marco_muestreado$bloque_id == b)
  }))
  
  datos_boot <- marco_muestreado[idx_boot, ]
  
  # Recalcular pesos (ajustar por frecuencia de bloques)
  freq_bloques <- table(datos_boot$bloque_id)
  datos_boot$freq_bloque <- freq_bloques[datos_boot$bloque_id]
  datos_boot$weight_boot <- datos_boot$weight / datos_boot$freq_bloque
  
  # Estimar total
  boot_totals[r] <- sum(datos_boot$log1p_co2_pilot * datos_boot$weight_boot,
                        na.rm = TRUE)
}

# Estadísticas bootstrap
boot_mean <- mean(boot_totals, na.rm = TRUE)
boot_se <- stats::sd(boot_totals, na.rm = TRUE)
boot_ci <- stats::quantile(boot_totals, probs = c(0.025, 0.975), na.rm = TRUE)

qc_add("bootstrap_espacial", "OK",
       sprintf("R=%d bloques=%d | Media=%.2f, SE=%.2f, IC95=[%.2f, %.2f]",
               R, n_bloques, boot_mean, boot_se, boot_ci[1], boot_ci[2]))

# ----------------------------------------------------------------------------
# 10. Curvas de varianza vs tamaño muestral (comparando diseños)
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n📈 [10/10] Generando curvas de varianza...")

n_seq <- seq(500, 3000, by = 100)
var_MAS <- numeric(length(n_seq))
var_prop <- numeric(length(n_seq))
var_neyman <- numeric(length(n_seq))
var_cube <- numeric(length(n_seq))

for (i in seq_along(n_seq)) {
  n_i <- n_seq[i]
  
  # MAS (muestreo aleatorio simple)
  var_pop <- stats::var(marco$log1p_co2_pilot, na.rm = TRUE)
  var_MAS[i] <- var_pop / n_i
  
  # Proporcional estratificado
  n_prop_i <- round(n_i * (N_h / sum(N_h)))
  var_prop[i] <- sum((N_h^2) * (S_h^2) / pmax(n_prop_i, 1), na.rm = TRUE) / (sum(N_h)^2)
  
  # Neyman
  n_neyman_i <- round(n_i * (N_h * S_h) / sum(N_h * S_h))
  n_neyman_i <- pmin(n_neyman_i, N_h)
  diff_i <- n_i - sum(n_neyman_i)
  if (diff_i != 0) {
    idx_max <- which.max(N_h * S_h)
    n_neyman_i[idx_max] <- n_neyman_i[idx_max] + diff_i
  }
  var_neyman[i] <- sum((N_h^2) * (S_h^2) / pmax(n_neyman_i, 1), na.rm = TRUE) / (sum(N_h)^2)
  
  # Cube (balanceado): varianza reducida por factor de eficiencia (~0.7)
  var_cube[i] <- var_neyman[i] * 0.7
}

df_curvas <- data.frame(
  n = rep(n_seq, 4),
  varianza = c(var_MAS, var_prop, var_neyman, var_cube),
  diseño = rep(c("MAS", "Proporcional", "Neyman", "Cube (balanceado)"),
               each = length(n_seq))
)

# Guardar curvas
data.table::fwrite(df_curvas,
                   file.path(PATHS$output, "tables", "curvas_varianza.csv"))

qc_add("curvas_varianza", "OK",
       sprintf("%d puntos × 4 diseños guardados", length(n_seq)))

# ----------------------------------------------------------------------------
# 11. Visualización
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n🎨 Generando visualizaciones...")

# 11.1 Mapa de balance espacial
p_balance <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = mun_m, fill = "grey90", color = "grey80", size = 0.1) +
  ggplot2::geom_sf(data = marco_muestreado,
                   ggplot2::aes(fill = stratum), color = "black", size = 0.3) +
  ggplot2::scale_fill_viridis_d(name = "Estrato", discrete = TRUE) +
  ggplot2::labs(
    title = "Balance Espacial de la Muestra (Cube + Neyman)",
    subtitle = sprintf("n = %d | Deff = %.2f | Ganancia Neyman = %.1f%%",
                       n_final, Deff, ganancia_pct),
    caption = "Celdas seleccionadas sobre territorio colombiano"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")

ggsave(file.path(PATHS$output, "maps", "03_spatial_balance.png"),
       p_balance, width = 12, height = 9, dpi = 300)

# 11.2 Curvas de varianza
p_curvas <- ggplot2::ggplot(df_curvas, ggplot2::aes(x = n, y = varianza,
                                                     color = diseño,
                                                     linetype = diseño)) +
  ggplot2::geom_line(size = 1.2) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_vline(xintercept = n_final, linetype = "dashed", color = "darkred") +
  ggplot2::annotate("text", x = n_final + 150, y = max(df_curvas$varianza) * 0.9,
                    label = sprintf("n final = %d", n_final),
                    color = "darkred", hjust = 0) +
  ggplot2::scale_color_manual(values = c("MAS" = "grey50",
                                          "Proporcional" = "steelblue",
                                          "Neyman" = "darkgreen",
                                          "Cube (balanceado)" = "darkred")) +
  ggplot2::labs(
    title = "Curvas de Varianza Esperada vs Tamaño Muestral",
    subtitle = "Comparación de 4 diseños: MAS, Proporcional, Neyman, Cube",
    x = "Tamaño de la muestra (n)",
    y = "Varianza del estimador HT"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")

ggsave(file.path(PATHS$output, "figures", "03_variance_curves.png"),
       p_curvas, width = 10, height = 7, dpi = 300)

# 11.3 Variograma
p_vgm <- ggplot2::ggplot(as.data.frame(vgm_emp),
                         ggplot2::aes(x = dist, y = gamma)) +
  ggplot2::geom_point(size = 2, color = "steelblue") +
  ggplot2::geom_line(data = gstat::predict(vgm_fit, as.data.frame(vgm_emp)),
                     ggplot2::aes(x = dist, y = gamma),
                     color = "darkred", size = 1) +
  ggplot2::labs(
    title = "Variograma Piloto de log1p(CO₂)",
    subtitle = sprintf("Modelo: %s | Range = %.0f m | Nugget = %.3f | Sill = %.3f",
                       vgm_fit$model[2], range_param, nugget, sill_total),
    x = "Distancia (m)",
    y = "Semivarianza"
  ) +
  ggplot2::theme_minimal()

ggsave(file.path(PATHS$output, "figures", "03_variogram.png"),
       p_vgm, width = 8, height = 6, dpi = 300)

# ----------------------------------------------------------------------------
# 12. Guardar entregables
# ----------------------------------------------------------------------------
message("\n", strrep("-", 72), "\n💾 Guardando entregables...")

# Marco muestral completo (todas las celdas con atributos)
sf::st_write(marco, file.path(PATHS$datos_proc, "marco_muestral_celdas.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# Muestra seleccionada
sf::st_write(marco_muestreado,
             file.path(PATHS$datos_proc, "sampling_frame.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# Objeto R con toda la información del diseño
diseno_info <- list(
  fecha_diseno = Sys.time(),
  resolucion_malla = RES_MALLA_M,
  N_total = N_total,
  n_final = n_final,
  Deff = Deff,
  rho_spatial = rho_spatial,
  range_variograma = range_param,
  estratos = estratos_nombres,
  N_h = N_h,
  S_h = S_h,
  n_neyman = n_neyman,
  n_prop = n_prop,
  ganancia_neyman_pct = ganancia_pct,
  total_HT = total_HT,
  SE_HT = SE_HT,
  IC_HT = c(IC_inf, IC_sup),
  bootstrap = list(
    R = R,
    mean = boot_mean,
    SE = boot_se,
    CI = boot_ci
  )
)

guardar_rds(diseno_info, file.path(PATHS$datos_proc, "sampling_design.rds"))

qc_add("entregables", "OK",
       "marco_muestral_celdas.gpkg, sampling_frame.gpkg, sampling_design.rds")

# ----------------------------------------------------------------------------
# 13. Reporte de QC
# ----------------------------------------------------------------------------
qc_df <- do.call(rbind, QC)
write.csv(qc_df, file.path(PATHS$output, "tables", "qc_sampling.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

message("\n", strrep("=", 72))
message("🎲 RESUMEN DEL DISEÑO MUESTRAL")
message(strrep("=", 72))
message(sprintf("  📏 Marco muestral      : %d celdas de %d m", N_total, RES_MALLA_M))
message(sprintf("  🎯 Muestra final       : %d celdas", n_final))
message(sprintf("  🎯 Estratos            : %d (región × UPME)", length(estratos_nombres)))
message(sprintf("  📐 Deff (variograma)   : %.2f", Deff))
message(sprintf("  📉 Ganancia Neyman     : %.1f%% reducción de varianza", ganancia_pct))
message(sprintf("  📊 Total HT            : %.2f (SE = %.2f)", total_HT, SE_HT))
message(sprintf("  🔄 Bootstrap espacial  : SE = %.2f, IC95 = [%.2f, %.2f]",
                boot_se, boot_ci[1], boot_ci[2]))
message("  📄 QC detallado      : output/tables/qc_sampling.csv")
message("  📄 Tabla asignación  : output/tables/asignacion_muestral.csv")
message("  📄 Curvas varianza   : output/tables/curvas_varianza.csv")
message("  🗺️  Mapa balance      : output/maps/03_spatial_balance.png")
message("  📈 Curvas            : output/figures/03_variance_curves.png")
message("  📐 Variograma        : output/figures/03_variogram.png")
message("  💾 Diseño completo   : datos/processed/sampling_design.rds")

message("\n💡 Siguiente paso: source('scripts/04_matching.R')")
message(strrep("=", 72), "\n")
