#!/usr/bin/env Rscript
# ============================================================================
# 01b_EDA.R — Análisis Exploratorio de Datos (EDA)
# ============================================================================
# Proyecto : Análisis Espacial de CO₂ y Evaluación Cuasi-Experimental de FNCER
#            en Colombia
# Autor    : Grupo de Investigación en Estadística Aplicada
# Fecha    : 2026-11-16
# Propósito:
#   1. Cargar los 11 datasets procesados desde datos/processed/
#   2. Evaluar calidad de datos (NA, tipos, rangos, outliers)
#   3. Triangular fuentes de emisiones (EDGAR vs OWID vs INGEI)
#   4. Explorar la matriz energética colombiana (UPME)
#   5. Analizar la penetración FNCER por tecnología y departamento
#   6. Estudiar emisiones evitadas por energías renovables (Ember)
#   7. Examinar correlaciones y estructura de covarianza
#   8. Visualizar el panel nacional con tratamientos T0–T4
#   9. Generar un dashboard consolidado
#  10. Exportar tablas resumen y figuras para el manuscrito
#
# Uso      : Rscript scripts/01b_EDA.R
#            source("scripts/01b_EDA.R")
#
# Entradas : datos/processed/*.rds
# Salidas  : output/figures/eda_*.png
#            output/tables/eda_*.html
#            datos/processed/eda_resumen.rds
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Cargar configuración del proyecto
# ----------------------------------------------------------------------------
source(file.path(Sys.getenv("CO2_PROJECT_ROOT", getwd()),
                 "scripts", "00_setup.R"))
library("purrr")
library("tidyverse")
log_msg("Iniciando ANÁLISIS EXPLORATORIO DE DATOS (EDA)", "INFO")
inicio <- Sys.time()

# ----------------------------------------------------------------------------
# 0.1. Detectar carpeta de datos procesados (flexible)
# ----------------------------------------------------------------------------
# Soporta tanto "datos/processed" como "datos/datos_procesados"
carpeta_proc <- c(
  file.path(RAIZ_PROYECTO, "datos", "processed"),
  file.path(RAIZ_PROYECTO, "datos", "datos_procesados")
)
carpeta_proc <- carpeta_proc[dir.exists(carpeta_proc)][1]

if (is.na(carpeta_proc)) {
  stop("No se encontró carpeta de datos procesados. Ejecuta 01_preprocessing.R.")
}
log_msg(paste("Carpeta de datos procesados:", carpeta_proc), "INFO")

# ============================================================================
# MÓDULO A — CARGA DE DATASETS PROCESADOS
# ============================================================================
log_msg("MÓDULO A: Cargando datasets procesados", "INFO")

cargar_rds_seguro <- function(nombre) {
  path <- file.path(carpeta_proc, paste0(nombre, ".rds"))
  if (!file.exists(path)) {
    log_msg(paste("Dataset no encontrado:", nombre), "WARN")
    return(NULL)
  }
  obj <- readRDS(path)
  log_msg(sprintf("  ✓ %-32s : %d filas × %d cols",
                  nombre, nrow(obj), ncol(obj)), "OK")
  return(obj)
}

# --- Fuentes de emisiones ---------------------------------------------------
edgar_co2_nacional       <- cargar_rds_seguro("edgar_co2_nacional")
owid_colombia            <- cargar_rds_seguro("owid_colombia")
ingei_energia_serie      <- cargar_rds_seguro("ingei_energia_serie")
climate_trace_colombia   <- cargar_rds_seguro("climate_trace_colombia")

# --- Matriz energética ------------------------------------------------------
upme_generacion_grupo    <- cargar_rds_seguro("upme_generacion_grupo")
upme_generacion_fuente   <- cargar_rds_seguro("upme_generacion_fuente")
ember_indicadores        <- cargar_rds_seguro("ember_indicadores_colombia")

# --- FNCER ------------------------------------------------------------------
minenergia_proyectos     <- cargar_rds_seguro("minenergia_meta_fncer_proyectos")
meta_fncer_tecnologia    <- cargar_rds_seguro("meta_fncer_tecnologia")

# --- Emisiones evitadas -----------------------------------------------------
ember_evitadas_estimadas <- cargar_rds_seguro("ember_evitadas_estimadas")

# --- Consolidados -----------------------------------------------------------
covariables_nacionales   <- cargar_rds_seguro("covariables_nacionales")
panel_nacional           <- cargar_rds_seguro("panel_nacional_completo")

# ============================================================================
# MÓDULO B — RESUMEN GENERAL Y CALIDAD DE DATOS
# ============================================================================
log_msg("MÓDULO B: Resumen de calidad de datos", "INFO")

# --- B1. Tabla resumen de todos los datasets --------------------------------
datasets_lista <- list(
  "EDGAR CO2 nacional"          = edgar_co2_nacional,
  "OWID Colombia"               = owid_colombia,
  "INGEI energía"               = ingei_energia_serie,
  "Climate TRACE activos"       = climate_trace_colombia,
  "UPME generación grupo"       = upme_generacion_grupo,
  "UPME generación fuente"      = upme_generacion_fuente,
  "Ember indicadores"           = ember_indicadores,
  "MinEnergía proyectos FNCER"  = minenergia_proyectos,
  "FNCER por tecnología"        = meta_fncer_tecnologia,
  "Ember evitadas estimadas"    = ember_evitadas_estimadas,
  "Covariables nacionales"      = covariables_nacionales,
  "Panel nacional"              = panel_nacional
)

resumen_datasets <- map_dfr(names(datasets_lista), function(nombre) {
  df <- datasets_lista[[nombre]]
  if (is.null(df)) {
    return(tibble(
      dataset = nombre, n_filas = NA_integer_, n_cols = NA_integer_,
      pct_na = NA_real_, rango_anios = NA_character_
    ))
  }
  # Detectar columna de año
  col_anio <- intersect(c("year", "anio"), names(df))[1]
  rango <- if (!is.na(col_anio)) {
    sprintf("%d-%d", min(df[[col_anio]], na.rm = TRUE),
            max(df[[col_anio]], na.rm = TRUE))
  } else NA_character_
  
  tibble(
    dataset   = nombre,
    n_filas   = nrow(df),
    n_cols    = ncol(df),
    pct_na    = round(100 * mean(is.na(df)), 2),
    rango_anios = rango
  )
})

log_msg("B1: Resumen de datasets cargados:", "INFO")
print(resumen_datasets)

# --- B2. Exportar tabla resumen ---------------------------------------------
path_tabla_resumen <- file.path(PATHS$output, "tables",
                                "eda_resumen_datasets.html")
kableExtra::kable(resumen_datasets, format = "html",
                  caption = "Resumen de datasets procesados") %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"),
                            full_width = FALSE) %>%
  kableExtra::save_kable(path_tabla_resumen)

# --- B3. Matriz de NA por dataset -------------------------------------------
log_msg("B3: Matriz de datos faltantes por columna...", "INFO")

matriz_na <- map_dfr(names(datasets_lista), function(nombre) {
  df <- datasets_lista[[nombre]]
  if (is.null(df)) return(NULL)
  tibble(
    dataset  = nombre,
    columna  = names(df),
    n_na     = map_int(df, ~sum(is.na(.x))),
    pct_na   = map_dbl(df, ~round(100 * mean(is.na(.x)), 2))
  ) %>%
    filter(pct_na > 0)
})

if (nrow(matriz_na) > 0) {
  log_msg(sprintf("B3: %d columnas con NA > 0%%", nrow(matriz_na)), "WARN")
  top_na <- matriz_na %>% arrange(desc(pct_na)) %>% slice_head(n = 15)
  print(top_na)
}

# ============================================================================
# MÓDULO C — TRIANGULACIÓN DE FUENTES DE EMISIONES
# ============================================================================
log_msg("MÓDULO C: Triangulación EDGAR vs OWID vs INGEI", "INFO")

# --- C1. Unir series de emisiones nacionales --------------------------------
# EDGAR: co2_total_t (toneladas) | OWID: co2 (Mt) | INGEI: energia_total_ktco2eq (kt)
emisiones_compuestas <- tibble(year = 1990:2025) %>%
  left_join(edgar_co2_nacional %>%
              mutate(co2_edgar_mt = co2_total_t / 1e6) %>%
              select(year, co2_edgar_mt),
            by = "year") %>%
  left_join(owid_colombia %>%
              select(year, co2_owid_mt = co2),
            by = "year") %>%
  left_join(ingei_energia_serie %>%
              mutate(co2_ingei_mt = energia_total_ktco2eq / 1e3) %>%
              select(year, co2_ingei_mt),
            by = "year")

# --- C2. Gráfico de triangulación -------------------------------------------
emisiones_larga <- emisiones_compuestas %>%
  pivot_longer(starts_with("co2_"),
               names_to = "fuente", values_to = "co2_mt") %>%
  mutate(fuente = str_replace(fuente, "co2_", "") %>%
           str_replace("_mt", "") %>%
           str_to_upper())

p_emisiones <- ggplot(emisiones_larga, aes(x = year, y = co2_mt,
                                           color = fuente, group = fuente)) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(shape = fuente), size = 2.5) +
  scale_color_manual(values = c("EDGAR" = "#d73027",
                                 "OWID"  = "#3288bd",
                                 "INGEI" = "#66c2a5")) +
  labs(
    title = "Emisiones de CO₂ en Colombia — Triangulación de fuentes",
    subtitle = "Comparación EDGAR (UE) vs OWID vs INGEI (Colombia)",
    x = "Año", y = "Mt CO₂", color = "Fuente", shape = "Fuente"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

path_fig_emisiones <- file.path(PATHS$output, "figures",
                                "eda_emisiones_triangulacion.png")
ggsave(path_fig_emisiones, p_emisiones, width = 12, height = 7, dpi = 300)
log_msg(sprintf("C2: Figura guardada: %s", basename(path_fig_emisiones)), "OK")

# --- C3. Correlación entre fuentes ------------------------------------------
cor_emisiones <- emisiones_compuestas %>%
  select(starts_with("co2_")) %>%
  cor(use = "pairwise.complete.obs")

log_msg("C3: Matriz de correlación entre fuentes:", "INFO")
print(round(cor_emisiones, 3))

# ============================================================================
# MÓDULO D — MATRIZ ENERGÉTICA (UPME)
# ============================================================================
log_msg("MÓDULO D: Explorando la matriz energética (UPME)", "INFO")

# --- D1. Evolución de la generación por grupo --------------------------------
p_grupo <- ggplot(upme_generacion_grupo,
                  aes(x = year, y = total_gwh / 1000, fill = "total")) +
  geom_col(color = "black", alpha = 0.3) +
  geom_line(aes(y = hidro_gwh / 1000), color = PALETA_PROYECTO["hidro"],
            linewidth = 1.2) +
  geom_line(aes(y = fosil_gwh / 1000), color = PALETA_PROYECTO["fossil"],
            linewidth = 1.2) +
  geom_line(aes(y = solar_gwh / 1000), color = PALETA_PROYECTO["solar"],
            linewidth = 1.2) +
  geom_line(aes(y = eolica_gwh / 1000), color = PALETA_PROYECTO["eolica"],
            linewidth = 1.2) +
  geom_line(aes(y = biomasa_gwh / 1000), color = PALETA_PROYECTO["biomasa"],
            linewidth = 1.2) +
  labs(
    title = "Generación eléctrica por fuente — Colombia (UPME)",
    subtitle = "TWh por año (2010–2025)",
    x = "Año", y = "TWh"
  ) +
  theme_minimal(base_size = 12)

path_fig_grupo <- file.path(PATHS$output, "figures",
                            "eda_upme_grupo.png")
ggsave(path_fig_grupo, p_grupo, width = 12, height = 7, dpi = 300)

# --- D2. Participación de renovables (%) -------------------------------------
p_pct <- ggplot(upme_generacion_grupo, aes(x = year)) +
  geom_line(aes(y = pct_renovable, color = "Renovable"), linewidth = 1.2) +
  geom_line(aes(y = pct_fosil, color = "Fósil"), linewidth = 1.2) +
  geom_hline(yintercept = 30, linetype = "dashed", color = "gray40") +
  annotate("text", x = Inf, y = 31, label = "Meta PDES 30%",
           hjust = 1.05, size = 3.5, color = "gray40") +
  scale_color_manual(values = c("Renovable" = "#66c2a5", "Fósil" = "#d73027")) +
  labs(
    title = "Participación renovable vs fósil en la generación eléctrica",
    x = "Año", y = "%", color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

path_fig_pct <- file.path(PATHS$output, "figures", "eda_upme_pct_renovable.png")
ggsave(path_fig_pct, p_pct, width = 11, height = 6, dpi = 300)

# --- D3. Descomposición por fuente (UPME mensual agregado) -------------------
if (!is.null(upme_generacion_fuente)) {
  p_fuente <- ggplot(upme_generacion_fuente,
                     aes(x = year, y = generacion_anual_gwh,
                         fill = fuente)) +
    geom_col(width = 0.85) +
    scale_fill_brewer(palette = "Set2", name = "Fuente") +
    labs(
      title = "Generación anual por fuente (UPME mensual agregado)",
      x = "Año", y = "GWh"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  path_fig_fuente <- file.path(PATHS$output, "figures",
                               "eda_upme_fuente.png")
  ggsave(path_fig_fuente, p_fuente, width = 12, height = 7, dpi = 300)
}

# ============================================================================
# MÓDULO E — FNCER: PROYECTOS Y TECNOLOGÍAS
# ============================================================================
log_msg("MÓDULO E: Explorando FNCER (MinEnergía)", "INFO")

# --- E1. Proyectos por tipo de tecnología ------------------------------------
if (!is.null(minenergia_proyectos)) {
  p_tipo <- ggplot(minenergia_proyectos, aes(x = reorder(tipo, -capacidad_mw))) +
    geom_bar(aes(fill = tipo), show.legend = FALSE) +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      title = "Proyectos FNCER por tipo de tecnología",
      subtitle = sprintf("Total: %d proyectos | Capacidad: %.1f MW",
                         nrow(minenergia_proyectos),
                         sum(minenergia_proyectos$capacidad_mw, na.rm = TRUE)),
      x = "Tecnología", y = "Número de proyectos"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  path_fig_tipo <- file.path(PATHS$output, "figures", "eda_fncer_tipo.png")
  ggsave(path_fig_tipo, p_tipo, width = 11, height = 6, dpi = 300)
  
  # --- E2. Distribución geográfica por departamento --------------------------
  depto_resumen <- minenergia_proyectos %>%
    group_by(departamento) %>%
    summarise(
      n_proyectos = n(),
      capacidad_mw = sum(capacidad_mw, na.rm = TRUE),
      inversion_cop = sum(inversion_cop, na.rm = TRUE),
      empleos = sum(empleos, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_proyectos)) %>%
    slice_head(n = 15)
  
  p_depto <- ggplot(depto_resumen,
                    aes(x = reorder(departamento, -n_proyectos),
                        y = n_proyectos, fill = capacidad_mw)) +
    geom_col() +
    scale_fill_viridis_c(name = "Capacidad (MW)") +
    labs(
      title = "Top 15 departamentos por proyectos FNCER",
      x = "Departamento", y = "Número de proyectos"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  path_fig_depto <- file.path(PATHS$output, "figures", "eda_fncer_departamento.png")
  ggsave(path_fig_depto, p_depto, width = 12, height = 7, dpi = 300)
  
  # --- E3. Relación inversión vs capacidad -----------------------------------
  p_inv <- ggplot(minenergia_proyectos,
                  aes(x = capacidad_mw, y = inversion_cop / 1e9,
                      color = tipo)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linetype = "dashed") +
    scale_color_brewer(palette = "Dark2", name = "Tecnología") +
    scale_x_log10() +
    labs(
      title = "Inversión vs capacidad instalada en proyectos FNCER",
      x = "Capacidad (MW, escala log)",
      y = "Inversión estimada (miles de millones COP)"
    ) +
    theme_minimal(base_size = 12)
  
  path_fig_inv <- file.path(PATHS$output, "figures", "eda_fncer_inversion.png")
  ggsave(path_fig_inv, p_inv, width = 11, height = 7, dpi = 300)
}

# --- E4. Evolución temporal de FNCER por tecnología --------------------------
if (!is.null(meta_fncer_tecnologia)) {
  p_evolutivo <- ggplot(meta_fncer_tecnologia,
                        aes(x = year, y = capacidad_mw, color = tipo)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Dark2", name = "Tecnología") +
    labs(
      title = "Evolución de capacidad FNCER por tecnología",
      subtitle = "Año de fecha estimada FPO (First Operation)",
      x = "Año", y = "Capacidad (MW)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  path_fig_evo <- file.path(PATHS$output, "figures", "eda_fncer_evolutivo.png")
  ggsave(path_fig_evo, p_evolutivo, width = 12, height = 7, dpi = 300)
}

# ============================================================================
# MÓDULO F — EMISIONES EVITADAS (EMBER)
# ============================================================================
log_msg("MÓDULO F: Emisiones evitadas por renovables (Ember)", "INFO")

if (!is.null(ember_evitadas_estimadas)) {
  # --- F1. Emisiones evitadas totales -----------------------------------------
  p_evitadas <- ggplot(ember_evitadas_estimadas,
                       aes(x = year, y = co2_evitado_mt)) +
    geom_col(fill = "#66c2a5", color = "white") +
    geom_smooth(method = "loess", se = FALSE, color = "#3288bd",
                linewidth = 1, linetype = "dashed") +
    labs(
      title = "Emisiones de CO₂ evitadas por energías renovables no convencionales",
      subtitle = "Escenario A (inventario) + Escenario B (modelo) — Ember",
      x = "Año", y = "Mt CO₂e evitadas"
    ) +
    theme_minimal(base_size = 12)
  
  path_fig_evit <- file.path(PATHS$output, "figures", "eda_emisiones_evitadas.png")
  ggsave(path_fig_evit, p_evitadas, width = 11, height = 6, dpi = 300)
  
  # --- F2. Descomposición por tecnología -------------------------------------
  evitadas_larga <- ember_evitadas_estimadas %>%
    select(year, solar_gwh, eolica_gwh, biomasa_gwh, fncer_no_hidro_gwh) %>%
    pivot_longer(-year, names_to = "tecnologia", values_to = "gwh")
  
  p_descomp <- ggplot(evitadas_larga,
                      aes(x = year, y = gwh, fill = tecnologia)) +
    geom_col() +
    scale_fill_manual(values = c("solar_gwh" = "#fee08b",
                                  "eolica_gwh" = "#66c2a5",
                                  "biomasa_gwh" = "#8c510a",
                                  "fncer_no_hidro_gwh" = "#d73027")) +
    labs(
      title = "Generación FNCER no-hidráulica por tecnología (GWh)",
      x = "Año", y = "GWh", fill = "Tecnología"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  path_fig_descomp <- file.path(PATHS$output, "figures",
                                "eda_evitadas_descomposicion.png")
  ggsave(path_fig_descomp, p_descomp, width = 11, height = 6, dpi = 300)
}

# ============================================================================
# MÓDULO G — CORRELACIONES Y ESTRUCTURA DE COVARIANZA
# ============================================================================
log_msg("MÓDULO G: Matriz de correlaciones", "INFO")

if (!is.null(covariables_nacionales)) {
  # --- G1. Seleccionar variables numéricas clave -----------------------------
  vars_corr <- c("co2_total_t", "co2", "co2_per_capita", "co2_per_gdp",
                 "energy_per_capita", "primary_energy_consumption",
                 "total_gwh", "pct_renovable", "pct_fosil",
                 "fncer_no_hidro_gwh", "solar_gwh", "eolica_gwh",
                 "biomasa_gwh", "hidro_gwh", "fosil_gwh")
  
  vars_disponibles <- intersect(vars_corr, names(covariables_nacionales))
  
  if (length(vars_disponibles) >= 3) {
    matriz_cor <- covariables_nacionales %>%
      select(all_of(vars_disponibles)) %>%
      cor(use = "pairwise.complete.obs")
    
    # --- G2. Heatmap de correlaciones ----------------------------------------
    cor_df <- as.data.frame(matriz_cor) %>%
      rownames_to_column("variable") %>%
      pivot_longer(-variable, names_to = "pares", values_to = "rho")
    
    p_cor <- ggplot(cor_df, aes(x = variable, y = pares, fill = rho)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
      scale_fill_gradient2(low = "#d73027", mid = "white", high = "#3288bd",
                           midpoint = 0, limits = c(-1, 1), name = "ρ") +
      labs(
        title = "Matriz de correlaciones — Covariables nacionales",
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank())
    
    path_fig_cor <- file.path(PATHS$output, "figures", "eda_correlaciones.png")
    ggsave(path_fig_cor, p_cor, width = 12, height = 10, dpi = 300)
    
    log_msg("G2: Heatmap de correlaciones guardado", "OK")
  }
}

# ============================================================================
# MÓDULO H — PANEL NACIONAL CON TRATAMIENTOS
# ============================================================================
log_msg("MÓDULO H: Panel nacional y tratamientos T0–T4", "INFO")

if (!is.null(panel_nacional)) {
  # --- H1. Distribución de tratamientos por año ------------------------------
  p_trat <- ggplot(panel_nacional, aes(x = year, y = co2_total_t / 1e6,
                                        color = tratamiento)) +
    geom_line(aes(group = 1), color = "gray50", linewidth = 0.8) +
    geom_point(size = 3.5) +
    scale_color_manual(values = c("T0" = PALETA_PROYECTO["fossil"],
                                   "T1" = PALETA_PROYECTO["solar"],
                                   "T2" = PALETA_PROYECTO["eolica"],
                                   "T3" = PALETA_PROYECTO["hidro"],
                                   "T4" = PALETA_PROYECTO["biomasa"]),
                       name = "Tratamiento dominante") +
    labs(
      title = "Panel nacional: emisiones CO₂ y tratamiento energético dominante",
      subtitle = "Tratamiento asignado por fuente con mayor generación en el año",
      x = "Año", y = "Mt CO₂"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  path_fig_trat <- file.path(PATHS$output, "figures", "eda_panel_tratamientos.png")
  ggsave(path_fig_trat, p_trat, width = 12, height = 7, dpi = 300)
  
  # --- H2. Tabla resumen de tratamientos -------------------------------------
  resumen_trat <- panel_nacional %>%
    group_by(tratamiento, nombre_tratamiento) %>%
    summarise(
      n_anios = n(),
      co2_media_mt = mean(co2_total_t, na.rm = TRUE) / 1e6,
      pct_renovable_media = mean(pct_renovable, na.rm = TRUE),
      pct_fosil_media = mean(pct_fosil, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(tratamiento)
  
  log_msg("H2: Resumen de tratamientos:", "INFO")
  print(as.data.frame(resumen_trat))
  
  path_tabla_trat <- file.path(PATHS$output, "tables",
                               "eda_resumen_tratamientos.html")
  kableExtra::kable(resumen_trat, format = "html",
                    caption = "Resumen por tratamiento dominante") %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"),
                              full_width = FALSE) %>%
    kableExtra::save_kable(path_tabla_trat)
}

# ============================================================================
# MÓDULO I — CLIMATE TRACE (si disponible)
# ============================================================================
log_msg("MÓDULO I: Climate TRACE — activos emisores", "INFO")

if (!is.null(climate_trace_colombia)) {
  
  # --- I1. Distribución por sector -------------------------------------------
  # Convertir sf a data.frame puro para ggplot
  ct_df <- if (inherits(climate_trace_colombia, "sf")) {
    sf::st_drop_geometry(climate_trace_colombia)
  } else {
    as.data.frame(climate_trace_colombia)
  }
  
  # Verificar que las columnas necesarias existen
  if (all(c("sector", "emissions_quantity") %in% names(ct_df))) {
    
    p_sector <- ggplot(ct_df,
                       aes(x = reorder(sector, -emissions_quantity),
                           y = emissions_quantity, fill = sector)) +
      geom_col(show.legend = FALSE) +
      scale_fill_brewer(palette = "Spectral") +
      labs(
        title = "Emisiones por sector — Climate TRACE Colombia",
        x = "Sector", y = "Emisiones (t CO₂e/año)"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    path_fig_sector <- file.path(PATHS$output, "figures",
                                 "eda_climate_trace_sector.png")
    ggsave(path_fig_sector, p_sector, width = 12, height = 7, dpi = 300)
    log_msg(sprintf("I1: Gráfico de sectores guardado: %s", path_fig_sector), "OK")
    
  } else {
    log_msg("I1: Columnas 'sector' o 'emissions_quantity' no encontradas", "WARN")
  }
  
  # --- I2. Mapa de puntos si tiene geometría ---------------------------------
  if (inherits(climate_trace_colombia, "sf")) {
    
    ct_pts <- climate_trace_colombia
    
    p_mapa <- ggplot(data = ct_pts) +
      geom_sf(aes(color = sector), size = 1.5, alpha = 0.7) +
      scale_color_brewer(palette = "Set1", name = "Sector") +
      labs(
        title = "Distribución espacial de activos emisores — Climate TRACE",
        subtitle = sprintf("%d activos en Colombia", nrow(ct_pts))
      ) +
      theme_minimal()
    
    path_fig_mapa <- file.path(PATHS$output, "maps",
                               "eda_climate_trace_mapa.png")
    ggsave(path_fig_mapa, p_mapa, width = 12, height = 8, dpi = 300)
    log_msg(sprintf("I2: Mapa de Climate TRACE guardado: %s", path_fig_mapa), "OK")
    
  } else {
    log_msg("I2: Climate TRACE no tiene geometría, omitiendo mapa", "INFO")
  }
  
} else {
  log_msg("MÓDULO I: Climate TRACE no disponible, omitiendo", "WARN")
}
# ============================================================================
# MÓDULO J — DASHBOARD CONSOLIDADO
# ============================================================================
log_msg("MÓDULO J: Generando dashboard EDA consolidado", "INFO")

# --- J1. Componer las figuras principales en un dashboard --------------------
figuras_existentes <- list.files(file.path(PATHS$output, "figures"),
                                  pattern = "eda_.*\\.png$", full.names = TRUE)

if (length(figuras_existentes) >= 2) {
  # Crear un PDF multi-página con todas las figuras
  path_dashboard <- file.path(PATHS$output, "figures", "eda_dashboard.pdf")
  
  pdf(path_dashboard, width = 12, height = 8)
  for (fig in figuras_existentes) {
    img <- png::readPNG(fig)
    grid::grid.raster(img)
  }
  dev.off()
  
  log_msg(sprintf("J1: Dashboard PDF con %d figuras: %s",
                  length(figuras_existentes), path_dashboard), "OK")
}

# ============================================================================
# MÓDULO K — EXPORTACIÓN DE RESULTADOS EDA
# ============================================================================
log_msg("MÓDULO K: Exportando resultados del EDA", "INFO")

resultados_eda <- list(
  fecha_eda = Sys.time(),
  resumen_datasets = resumen_datasets,
  matriz_na = matriz_na,
  emisiones_compuestas = emisiones_compuestas,
  cor_emisiones = cor_emisiones,
  resumen_tratamientos = if (!is.null(panel_nacional)) resumen_trat else NULL,
  depto_resumen = if (!is.null(minenergia_proyectos)) depto_resumen else NULL,
  figuras_generadas = figuras_existentes,
  tiempo_total_min = round(as.numeric(difftime(Sys.time(), inicio,
                                                units = "mins")), 2)
)

path_resultados <- file.path(carpeta_proc, "eda_resumen.rds")
saveRDS(resultados_eda, path_resultados, compress = "gzip")
log_msg(sprintf("K: Resultados EDA exportados: %s", path_resultados), "OK")

# ============================================================================
# MÓDULO L — RESUMEN FINAL
# ============================================================================
tiempo_total <- round(as.numeric(difftime(Sys.time(), inicio, units = "mins")), 2)

message("
", paste(rep("=", 72), collapse = ""))
message("ANÁLISIS EXPLORATORIO COMPLETADO")
message(paste(rep("=", 72), collapse = ""))
message("  Tiempo total        : ", tiempo_total, " min")
message("  Datasets analizados : ", sum(!sapply(datasets_lista, is.null)))
message("  Figuras generadas   : ", length(figuras_existentes))
message("  Tablas exportadas   : ",
        length(list.files(file.path(PATHS$output, "tables"),
                          pattern = "eda_.*\\.html$")))
message(paste(rep("=", 72), collapse = ""))

message("
Hallazgos clave:")
message("  1. Triangulación de emisiones: EDGAR vs OWID vs INGEI (ver figura C2)")
message("  2. Matriz energética: transición renovable vs fósil (ver figura D2)")
message("  3. FNCER: distribución geográfica y por tecnología (figuras E1–E4)")
message("  4. Emisiones evitadas: impacto cuantificado de renovables (figura F1)")
message("  5. Correlaciones: estructura multivariada (figura G2)")
message("  6. Tratamientos nacionales: asignación T0–T4 por año (figura H1)")

message("
Siguiente paso sugerido:")
message("   source('scripts/02_spatial_sampling.R')  # Diseño muestral GRTS")
message("   source('scripts/03_gwr.R')               # Regresión geográficamente ponderada")
message("   source('scripts/04_matching.R')          # Propensity Score Matching")

