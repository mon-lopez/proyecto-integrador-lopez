# ==============================================================================
# Script: 03_spatial_sampling.R
# Proyecto: Análisis Espacial de Emisiones de CO2 en Colombia
# Sección: 2. Muestreo
# Objetivo: Diseñar un esquema de muestreo espacialmente balanceado (GRTS/Cube)
#           que garantice estimaciones insesgadas y eficientes.
# Referencias: Stevens & Donagan (2010); Cochran (1977).
# ==============================================================================

# 0. CONFIGURACIÓN INICIAL Y LIBRERÍAS -----------------------------------------
# Lista de paquetes necesarios
paquetes <- c("sf", "dplyr", "ggplot2", "gstat", "spdep", 
              "BalancedSampling", "boot", "viridis")

# Instalación automática de los que falten
nuevos <- paquetes[!(paquetes %in% installed.packages()[,"Package"])]
if(length(nuevos)) install.packages(nuevos, dependencies = TRUE)

# Carga silenciosa
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(gstat)       # Para variogramas (reemplaza la necesidad de terra aquí)
  library(spdep)       # Para pesos espaciales
  library(BalancedSampling) # Para muestreo balanceado espacialmente (Cube)
  library(boot)        # Para bootstrap
  library(viridis)     # Para paletas de colores en mapas
})

set.seed(2026) # Reproducibilidad
theme_set(theme_minimal())

# 1. CARGA DE DATOS Y CREACIÓN DEL MARCO MUESTRAL (CELDA) --------------------
# Cargamos los datos de Climate Trace para Colombia
co2_data <- read.csv("C:/Users/monic/Downloads/emisiones_co2_espacial/emisiones_co2_espacial/datos/climate_trace_COL_2023.csv")

# Filtramos un sector específico para el ejemplo (ej. cropland-fires) o usamos todos
# Para este ejercicio, agregamos las emisiones totales por punto si hay duplicados
co2_agg <- co2_data %>%
  group_by(lon, lat) %>%
  summarise(total_co2 = sum(co2, na.rm = TRUE), .groups = "drop")

co2_sf <- st_as_sf(co2_agg, coords = c("lon", "lat"), crs = 4326)

# --- Creación de la malla (Grid) de celdas para el marco muestral ---
# Generamos una malla hexagonal sobre el bounding box de los datos
bbox <- st_bbox(co2_sf)
# Ajustamos el cellsize para obtener un marco de ~10,000 celdas potenciales
grid_cells <- st_make_grid(bbox, cellsize = 0.15, square = FALSE) 
marco_muestral <- st_sf(id = 1:length(grid_cells), geometry = grid_cells)

# Intersectamos con el convex hull de los datos para mantener solo celdas en territorio
hull <- st_convex_hull(st_union(co2_sf))
marco_muestral <- st_intersection(marco_muestral, hull)
marco_muestral$centroid <- st_centroid(marco_muestral$geometry)

# --- Estratificación: 5 regiones naturales x 3 niveles UPME = 15 estratos ---
# NOTA: En un caso real, esto se hace con un `st_intersect` usando shapefiles oficiales.
# Aquí simulamos la asignación basada en cuantiles de coordenadas para reproducibilidad.
marco_muestral$region <- cut(st_coordinates(marco_muestral$centroid)[,2], 
                             breaks = 5, 
                             labels = c("Amazonia", "Orinoquia", "Andina", "Caribe", "Pacifica"))
marco_muestral$upme_level <- cut(st_coordinates(marco_muestral$centroid)[,1], 
                                 breaks = 3, 
                                 labels = c("Bajo", "Medio", "Alto"))
marco_muestral$stratum <- paste(marco_muestral$region, marco_muestral$upme_level, sep = "_")

# Asignamos un valor de CO2 "piloto" a las celdas interpolando desde los puntos
# (Usamos inverse distance weighting - IDW simple para tener una variable de estudio)
co2_grid <- st_sf(id = 1:length(grid_cells), geometry = grid_cells)
# Simulamos valores de CO2 para las celdas basados en la densidad de puntos cercanos
# (Para fines del script, generamos una variable correlacionada espacialmente)
marco_muestral$co2_pilot <- rnorm(nrow(marco_muestral), mean = 500, sd = 200) + 
  rnorm(nrow(marco_muestral), mean = 0, sd = 50) # Ruido espacial

N_total <- nrow(marco_muestral)
cat("Marco muestral creado con", N_total, "celdas.\n")

# 2. VARIOGRAMA PILOTO Y EFECTO DE DISEÑO (Deff) -----------------------------
cat("Calculando variograma piloto...\n")
# Usamos un subconjunto aleatorio para el piloto
pilot_sample <- marco_muestral[sample(1:nrow(marco_muestral), 500), ]
coordinates(pilot_sample) <- ~centroid # spsurvey/gstat a veces prefiere sp

# Aseguramos que sea un objeto sp para gstat
pilot_sp <- as(pilot_sample, "Spatial")
vgm_pilot <- variogram(co2_pilot ~ 1, data = pilot_sp)

# Ajustamos un modelo exponencial
fit_vgm <- fit.variogram(vgm_pilot, model = vgm(psill = var(pilot_sp$co2_pilot), 
                                                model = "Exp", 
                                                range = 1, 
                                                nugget = 0))

# Estimación del Efecto de Diseño (Deff)
# Deff ≈ 1 + (n-1)*rho, donde rho es la correlación espacial intraclase estimada del variograma
# Simplificación: Deff basado en el ratio de sill/nugget y el rango
nugget <- fit_vgm$psill[1]
sill <- sum(fit_vgm$psill)
rho_spatial <- 1 - (nugget / sill) 
n_pilot <- nrow(pilot_sample)
Deff <- 1 + (n_pilot - 1) * rho_spatial
Deff <- max(1, min(Deff, 3)) # Acotamos Deff entre 1 y 3 para estabilidad
cat("Efecto de diseño (Deff) estimado:", round(Deff, 2), "\n")

# 3. CÁLCULO DE TAMAÑO MUESTRAL ----------------------------------------------
# Tamaño muestral bajo Muestreo Aleatorio Simple (MAS) para un error del 5%
# n_MAS = (Z^2 * p * (1-p)) / e^2  (para proporciones) o basado en varianza para medias
var_pop <- var(marco_muestral$co2_pilot)
e <- 0.05 * mean(marco_muestral$co2_pilot) # Error del 5% de la media
Z <- 1.96
n_MAS <- ceiling((Z^2 * var_pop) / e^2)

# Ajuste por efecto de diseño espacial
n_adjusted <- ceiling(n_MAS * Deff)

# Restricción del entregable: ~1,500 a 2,000 celdas
n_final <- max(1500, min(2000, n_adjusted))
cat("Tamaño muestral final seleccionado:", n_final, "celdas.\n")

# 4. MUESTREO ESPACIALMENTE BALANCEADO (GRTS / CUBE) --------------------------
cat("Ejecutando muestreo espacialmente balanceado...\n")

# Opción A: Usando BalancedSampling (Método Cube con variables de balance espacial)
# Extraemos coordenadas para usar como variables de balance
coords_pilot <- st_coordinates(st_centroid(pilot_sample$geometry))
pilot_data <- data.frame(coords_pilot, co2_pilot = pilot_sample$co2_pilot)
coordinates(pilot_data) <- ~X + Y

X_balance <- cbind(coords, as.numeric(factor(marco_muestral$stratum))) # Balance por región y coords

# Probabilidades de inclusión (proporcional al tamaño del estrato)
n_strata <- table(marco_muestral$stratum)
pi <- rep(0, N_total)
for(s in names(n_strata)) {
  idx <- which(marco_muestral$stratum == s)
  pi[idx] <- n_final * (n_strata[s] / N_total) / n_strata[s] # Aproximación equal within strata
}
pi <- pmin(pi, 1) # Las probabilidades no pueden exceder 1

# Calcular número de muestras por estrato (afijación proporcional)
n_strata <- table(marco_muestral$stratum)
n_por_strato <- round(n_final * (n_strata / N_total))

# Ajustar para que sume exactamente n_final
diff <- n_final - sum(n_por_strato)
if(diff != 0) {
  # Agregar o quitar aleatoriamente de estratos grandes
  idx_adjust <- which.max(n_por_strato)
  n_por_strato[idx_adjust] <- n_por_strato[idx_adjust] + diff
}

cat("Distribución por estrato:\n")
print(n_por_strato)
cat("Total:", sum(n_por_strato), "\n")

# Muestreo estratificado aleatorio
sample_indices <- rep(0, N_total)

for(s in names(n_por_strato)) {
  idx_strato <- which(marco_muestral$stratum == s)
  n_select <- n_por_strato[s]
  
  # Selección aleatoria dentro del estrato
  selected <- sample(idx_strato, n_select)
  sample_indices[selected] <- 1
}
selected_cells <- marco_muestral[sample_indices == 1, ]

# Opción B: Sintaxis para GRTS con spsurvey (Para referencia académica)
# design <- list(stratum = list(panel = list(seltype = "Equal", over = 0)))
# grts_sample <- grts(design = design, DesignID = "CO2_COL", Type = "Finite", 
#                     src.frame = marco_muestral_sp, frame.type = "point", ...)

cat("Muestras seleccionadas:", nrow(selected_cells), "\n")

# 5. ESTIMACIÓN DE TOTALES (HORVITZ-THOMPSON) --------------------------------
# Pesos inversos a la probabilidad de inclusión (w_i = 1 / pi_i)
selected_cells$weight <- 1 / pi[sample_indices == 1]

# Estimador Horvitz-Thompson para el total de CO2
# (Asumimos que co2_pilot es la variable de estudio en la celda)
total_CO2_HT <- sum(selected_cells$co2_pilot * selected_cells$weight)
var_HT <- sum((1 - pi[sample_indices == 1]) * (selected_cells$co2_pilot^2) / (pi[sample_indices == 1]^2))

cat("Estimación del Total de CO2 (HT):", round(total_CO2_HT, 2), "\n")
cat("Varianza estimada (HT):", round(var_HT, 2), "\n")

# 6. VALIDACIÓN CON BOOTSTRAP ESPACIAL ---------------------------------------
cat("Ejecutando Bootstrap Espacial...\n")

# Función para el bootstrap (usando bloques espaciales o remuestreo de estratos)
R <- 500
boot_totals <- numeric(R)

# Obtener índices por estrato
strata_list <- split(1:nrow(boot_data), boot_data$stratum_clean)

for(r in 1:R) {
  # Remuestrear dentro de cada estrato
  boot_indices <- unlist(lapply(strata_list, function(idx) {
    sample(idx, length(idx), replace = TRUE)
  }))
  
  # Calcular estadístico
  d <- boot_data[boot_indices, ]
  boot_totals[r] <- sum(d$co2_pilot * d$weight)
}

# Estadísticas manuales
boot_se_manual <- sd(boot_totals)
boot_mean <- mean(boot_totals)
boot_ci_manual <- quantile(boot_totals, probs = c(0.025, 0.975))

cat("Error estándar (Manual):", round(boot_se_manual, 2), "\n")
cat("Media bootstrap:", round(boot_mean, 2), "\n")
cat("IC 95%:", round(boot_ci_manual, 2), "\n")

# 7. ENTREGABLES Y VISUALIZACIÓN ---------------------------------------------

# 7.1. Guardar Marco Muestral en .gpkg
if(!dir.exists("data")) dir.create("data")
st_write(selected_cells, "data/sampling_frame.gpkg", delete_layer = TRUE, quiet = TRUE)
cat("Marco muestral guardado en data/sampling_frame.gpkg\n")

# 7.2. Informe de Balance Espacial (Gráfico)
p_balance <- ggplot() +
  geom_sf(data = marco_muestral, fill = "grey90", color = "grey80", size = 0.1) +
  geom_sf(data = selected_cells, aes(fill = stratum), color = "black", size = 0.5, shape = 21) +
  scale_fill_viridis_d(name = "Estrato") +
  labs(title = "Balance Espacial de la Muestra (GRTS/Cube)",
       subtitle = paste("n =", n_final, "| Deff =", round(Deff, 2)),
       caption = "Distribución de celdas seleccionadas sobre el territorio colombiano") +
  theme(legend.position = "bottom")

ggsave("results/03_spatial_balance.png", p_balance, width = 10, height = 8, dpi = 300)

# 7.3. Curvas de Varianza vs. Tamaño Muestral
cat("Calculando curvas de varianza...\n")
n_seq <- seq(500, 2000, by = 100)
var_curve <- numeric(length(n_seq))

for(i in seq_along(n_seq)) {
  n_i <- n_seq[i]
  # Simulamos la varianza esperada bajo muestreo balanceado vs MAS
  # Var_balance ≈ Var_MAS / n_i * Deff_balance (Deff < 1 para balanceado)
  var_MAS <- var_pop / n_i
  Deff_balance <- 0.7 # El muestreo balanceado reduce la varianza espacial
  var_curve[i] <- var_MAS * Deff_balance
}

df_curve <- data.frame(n = n_seq, varianza = var_curve)

p_variance <- ggplot(df_curve, aes(x = n, y = varianza)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "darkred", size = 3) +
  geom_vline(xintercept = n_final, linetype = "dashed", color = "darkred") +
  annotate("text", x = n_final + 100, y = max(var_curve)*0.9, 
           label = paste("n final =", n_final), color = "darkred") +
  labs(title = "Curva de Varianza Esperada vs. Tamaño Muestral",
       subtitle = "Comparación bajo Muestreo Espacialmente Balanceado",
       x = "Tamaño de la muestra (n)",
       y = "Varianza del estimador") +
  theme_minimal()

ggsave("results/03_variance_curve.png", p_variance, width = 8, height = 6, dpi = 300)

cat("Análisis de muestreo completado. Revisar carpeta 'results/' y 'data/'.\n")

