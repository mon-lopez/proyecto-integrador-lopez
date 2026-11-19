#!/usr/bin/env Rscript
# ============================================================================
# 01_preprocessing.R — Limpieza, transformación y uniones espaciales
# ============================================================================
# Proyecto : Análisis Espacial de CO₂ y Evaluación Cuasi-Experimental de FNCER
#            en Colombia
# Autor    : Grupo de Investigación en Estadística Aplicada
# Fecha    : 2026-11-16
# Propósito:
#   1. Normalizar identificadores geográficos (códigos DANE)
#   2. Limpiar cada dataset según sus columnas reales
#   3. Ejecutar uniones espaciales (point-in-polygon, zonal stats, joins)
#   4. Construir el panel municipal (municipio × año) con covariables
#   5. Asignar tratamientos T0–T4 según la matriz energética municipal
#   6. Exportar datasets procesados a datos/processed/
#
# Uso      : Rscript scripts/01_preprocessing.R
#            source("scripts/01_preprocessing.R")
#
# Entradas : datos/raw/**  (archivos .xlsx, .csv, .gpkg, .rds)
# Salidas  : datos/processed/*.rds, *.gpkg
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Cargar configuración del proyecto
# ----------------------------------------------------------------------------
source(file.path(Sys.getenv("CO2_PROJECT_ROOT", getwd()),
                 "scripts", "00_setup.R"))
library(tidyverse)

log_msg("Iniciando PREPROCESSING", "INFO")
inicio <- Sys.time()

# ============================================================================
# MÓDULO A — FUNCIONES DE LIMPIEZA GENÉRICAS
# ============================================================================

# --- Normalizar código DANE de municipio (10 dígitos: DDMMNNNNNN) -----------
normalizar_dane <- function(codigo) {
  if (is.null(codigo)) return(NA_character_)
  x <- as.character(codigo)
  x <- gsub("[^0-9]", "", x)
  x <- str_pad(x, width = 10, side = "left", pad = "0")
  x[nchar(x) != 10] <- NA_character_
  return(x)
}

# --- Extraer código de departamento (primeros 2 dígitos) --------------------
extraer_depto <- function(dane) {
  str_sub(dane, 1, 2)
}

# --- Limpiar nombres de columnas --------------------------------------------
limpiar_nombres <- function(df) {
  names(df) <- names(df) %>%
    str_to_lower() %>%
    str_replace_all("[áà]", "a") %>%
    str_replace_all("[éè]", "e") %>%
    str_replace_all("[íì]", "i") %>%
    str_replace_all("[óò]", "o") %>%
    str_replace_all("[úù]", "u") %>%
    str_replace_all("[ñ]", "n") %>%
    str_replace_all("[^a-z0-9_]", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
  df
}

# --- Convertir columnas numéricas con formato colombiano (1.234,56) ---------
convertir_numero_co <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- x %>%
    str_replace_all("\\.", "") %>%
    str_replace(",", ".") %>%
    str_replace_all("[^0-9.\\-]", "")
  as.numeric(x)
}

# --- Eliminar duplicados por clave -------------------------------------------
quitar_duplicados <- function(df, clave, nombre = "dataset") {
  n0 <- nrow(df)
  df <- df %>% distinct(across(all_of(clave)), .keep_all = TRUE)
  n1 <- nrow(df)
  if (n0 != n1) {
    log_msg(paste(nombre, ": eliminados", n0 - n1, "duplicados por",
                  paste(clave, collapse = "+")), "WARN")
  }
  return(df)
}

# --- Reportar calidad de datos -----------------------------------------------
reportar_calidad <- function(df, nombre) {
  n_na <- sum(is.na(df))
  n_total <- nrow(df) * ncol(df)
  pct <- round(100 * n_na / max(n_total, 1), 2)
  log_msg(sprintf("%s: %d filas × %d cols | NA: %d (%.2f%%)",
                  nombre, nrow(df), ncol(df), n_na, pct), "OK")
  invisible(df)
}

# ============================================================================
# MÓDULO B — LIMPIEZA DE DATOS DE EMISIONES
# ============================================================================
log_msg("MÓDULO B: Limpieza de datos de emisiones", "INFO")

# --- B1. EDGAR — Emisiones nacionales de CO2 --------------------------------
# Columnas reales: iso3, pais, anio, co2_gg (gigagramos = kilotoneladas)
edgar_co2_limpio <- tryCatch({
  df <- edgar_co2_pais %>% limpiar_nombres()
  
  # Filtrar solo Colombia
  if ("pais" %in% names(df)) {
    df <- df %>% filter(str_detect(tolower(pais), "colombia"))
  } else if ("iso3" %in% names(df)) {
    df <- df %>% filter(iso3 == "COL")
  }
  
  # co2_gg está en gigagramos = kilotoneladas; convertir a toneladas
  df <- df %>%
    mutate(
      year    = as.integer(anio),
      co2_gg  = as.numeric(co2_gg),
      co2_t   = co2_gg * 1e3   # de Gg a toneladas (1 Gg = 1000 t)
    ) %>%
    select(year, co2_gg, co2_t) %>%
    group_by(year) %>%
    summarise(co2_total_t = sum(co2_t, na.rm = TRUE), .groups = "drop") %>%
    arrange(year)
  
  reportar_calidad(df, "EDGAR CO2 nacional")
  df
}, error = function(e) {
  log_msg(paste("EDGAR: error de limpieza —", e$message), "ERROR")
  NULL
})

# --- B2. Climate TRACE — Activos emisores en Colombia ------------------------
# Columnas reales: id, name, sector, subsector, country, assetType, sourceType,
#                  gas, emissionsQuantity, emissionsFactor, emissionsFactorUnits,
#                  activity, activityUnits, capacity, capacityUnits,
#                  capacityFactor, year, lon, lat
climate_trace_limpio <- tryCatch({
  df <- climate_trace_col %>% limpiar_nombres()
  
  cols_clave <- c("id", "name", "sector", "subsector",
                  "emissionsquantity", "capacity", "lon", "lat")
  validar_columnas(df, cols_clave, "Climate TRACE")
  
  df <- df %>%
    mutate(
      lat = as.numeric(lat),
      lon = as.numeric(lon),
      emissions_quantity = as.numeric(emissionsquantity),
      capacity = as.numeric(capacity),
      year = if ("year" %in% names(df)) as.integer(year) else 2023L,
      sector = str_to_title(sector),
      asset_type = str_to_title(assettype)
    ) %>%
    filter(
      !is.na(lat), !is.na(lon),
      lat >= -4.5, lat <= 13.5,
      lon >= -82,  lon <= -66
    ) %>%
    quitar_duplicados("id", "Climate TRACE") %>%
    select(id, name, sector, subsector, asset_type,
           lat, lon, emissions_quantity, capacity, year)
  
  # Crear objeto sf de puntos
  ct_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = CRS$wgs84)
  
  reportar_calidad(df, "Climate TRACE activos")
  list(tib = df, sf = ct_sf)
}, error = function(e) {
  log_msg(paste("Climate TRACE: error —", e$message), "ERROR")
  NULL
})

# --- B3. OWID — Series históricas de CO2 de Colombia -------------------------
# Columnas reales: country, year, iso_code, population, gdp, cement_co2, co2,
# co2_growth_abs, co2_growth_prct, co2_including_luc, co2_per_capita,
# co2_per_gdp, co2_per_unit_energy, coal_co2, consumption_co2, cumulative_*,
# energy_per_capita, energy_per_gdp, flaring_co2, gas_co2, ghg_*, land_use_*,
# methane, nitrous_oxide, oil_co2, other_*, primary_energy_consumption,
# share_global_*, temperature_change_*, total_ghg, trade_co2, trade_co2_share
#
# Reducimos a las columnas relevantes para el análisis
owid_col_limpio <- tryCatch({
  df <- owid_co2_data %>% limpiar_nombres()
  
  # Seleccionar solo columnas clave para Colombia
  cols_relevantes <- c("country", "year", "iso_code", "population", "gdp",
                       "co2", "co2_per_capita", "co2_per_gdp",
                       "co2_including_luc", "co2_including_luc_per_capita",
                       "energy_per_capita", "energy_per_gdp",
                       "primary_energy_consumption",
                       "coal_co2", "oil_co2", "gas_co2", "cement_co2",
                       "consumption_co2", "trade_co2",
                       "total_ghg", "methane", "nitrous_oxide")
  
  cols_disponibles <- intersect(cols_relevantes, names(df))
  
  df <- df %>%
    filter(country == "Colombia") %>%
    select(all_of(cols_disponibles)) %>%
    mutate(
      year = as.integer(year),
      population = as.numeric(population),
      gdp = as.numeric(gdp),
      co2 = as.numeric(co2),                    # Mt CO2
      co2_per_capita = as.numeric(co2_per_capita),
      co2_per_gdp = as.numeric(co2_per_gdp),
      energy_per_capita = as.numeric(energy_per_capita),
      energy_per_gdp = as.numeric(energy_per_gdp),
      primary_energy_consumption = as.numeric(primary_energy_consumption),
      total_ghg = as.numeric(total_ghg)
    ) %>%
    arrange(year)
  
  reportar_calidad(df, "OWID Colombia (columnas reducidas)")
  df
}, error = function(e) {
  log_msg(paste("OWID: error —", e$message), "ERROR")
  NULL
})

# --- B4. Ember — Indicadores Colombia ----------------------------------------
# Columnas reales: anio, fuente, generacion_twh, capacidad_gw,
#                  emisiones_mtco2e, intensidad_gco2e_kwh
ember_indicadores_limpio <- tryCatch({
  df <- ember_indicadores %>% limpiar_nombres()
  
  df <- df %>%
    mutate(
      year = as.integer(anio),
      fuente = str_to_title(fuente),
      generacion_twh = as.numeric(generacion_twh),
      capacidad_gw = as.numeric(capacidad_gw),
      emisiones_mtco2e = as.numeric(emisiones_mtco2e),
      intensidad_gco2e_kwh = as.numeric(intensidad_gco2e_kwh)
    ) %>%
    select(year, fuente, generacion_twh, capacidad_gw,
           emisiones_mtco2e, intensidad_gco2e_kwh) %>%
    arrange(year, fuente)
  
  reportar_calidad(df, "Ember indicadores")
  df
}, error = function(e) {
  log_msg(paste("Ember indicadores: error —", e$message), "ERROR")
  NULL
})

# --- B5. Ember — Emisiones evitadas (reportadas) -----------------------------
# Columnas reales: anio, solar_gwh, eolica_gwh, biomasa_gwh, fncer_no_hidro_gwh,
#                  escA_mtco2eq, escB_mtco2eq, total_generacion_twh,
#                  intensidad_gco2e_kwh
ember_evitadas_reportadas <- tryCatch({
  df <- emisiones_evitadas_ember %>% limpiar_nombres()
  
  df <- df %>%
    mutate(
      year = as.integer(anio),
      solar_gwh = as.numeric(solar_gwh),
      eolica_gwh = as.numeric(eolica_gwh),
      biomasa_gwh = as.numeric(biomasa_gwh),
      fncer_no_hidro_gwh = as.numeric(fncer_no_hidro_gwh),
      escA_mtco2eq = as.numeric(esca_mtco2eq),
      escB_mtco2eq = as.numeric(escb_mtco2eq),
      total_generacion_twh = as.numeric(total_generacion_twh),
      co2_evitado_mt = escA_mtco2eq + escB_mtco2eq
    ) %>%
    select(year, solar_gwh, eolica_gwh, biomasa_gwh, fncer_no_hidro_gwh,
           escA_mtco2eq, escB_mtco2eq, co2_evitado_mt,
           total_generacion_twh) %>%
    arrange(year)
  
  reportar_calidad(df, "Ember emisiones evitadas (reportadas)")
  df
}, error = function(e) {
  log_msg(paste("Ember evitadas reportadas: error —", e$message), "ERROR")
  NULL
})

# --- B6. Ember — Emisiones evitadas estimadas --------------------------------
# Columnas reales: anio, anio_parcial, solar_gwh, eolica_gwh, biomasa_gwh,
#                  fncer_no_hidro_gwh, escA_inventario_tco2eq, escB_mdl_tco2eq,
#                  escA_mtco2eq, escB_mtco2eq

ember_evitadas_estimadas <- tryCatch({
  df <- ember_emisiones_evitadas %>% limpiar_nombres()
  
  df <- df %>%
    mutate(
      year = as.integer(anio),
      solar_gwh = as.numeric(solar_gwh),
      eolica_gwh = as.numeric(eolica_gwh),
      biomasa_gwh = as.numeric(biomasa_gwh),
      fncer_no_hidro_gwh = as.numeric(fncer_no_hidro_gwh),
      escA_mtco2eq = as.numeric(esca_mtco2eq),
      escB_mtco2eq = as.numeric(escb_mtco2eq),
      co2_evitado_mt = esca_mtco2eq + escb_mtco2eq
    ) %>%
    select(year, solar_gwh, eolica_gwh, biomasa_gwh, fncer_no_hidro_gwh,
           escA_mtco2eq, escB_mtco2eq, co2_evitado_mt) %>%
    arrange(year)
  
  reportar_calidad(df, "Ember emisiones evitadas (estimadas)")
  df
}, error = function(e) {
  log_msg(paste("Ember evitadas estimadas: error —", e$message), "ERROR")
  NULL
})

# ============================================================================
# MÓDULO C — LIMPIEZA DE DATOS FNCER Y ENERGÉTICOS
# ============================================================================
log_msg("MÓDULO C: Limpieza de datos FNCER y energéticos", "INFO")

# --- C1. UPME — Generación anual por grupo -----------------------------------
# Columnas reales: anio, meses_reportados, anio_parcial,
#                  Hidroeléctrica (convencional + PAH), Biomasa, Solar,
#                  Eólica, Fósil, Otros/no identificado, total_gwh,
#                  fncer_no_hidro_gwh, pct_fncer_no_hidro, pct_renovable_total,
#                  pct_fosil
upme_anual_limpio <- tryCatch({
  df <- upme_generacion_anual %>% limpiar_nombres()
  
  # Los nombres con tildes y espacios ya fueron limpiados por limpiar_nombres()
  # Hidroeléctrica (convencional + PAH) -> hidroelectrica_convencional__pah_
  # Eólica -> eolica
  # Fósil -> fosil
  
  df <- df %>%
    mutate(
      year = as.integer(anio),
      meses_reportados = as.integer(meses_reportados),
      # Buscar columnas dinámicamente (los nombres pueden variar tras limpieza)
      hidro_gwh = .data[[grep("hidro", names(df), value = TRUE)[1]]] %>% as.numeric(),
      biomasa_gwh = as.numeric(biomasa),
      solar_gwh = as.numeric(solar),
      eolica_gwh = as.numeric(eolica),
      fosil_gwh = as.numeric(fosil),
      otros_gwh = .data[[grep("otros", names(df), value = TRUE)[1]]] %>% as.numeric(),
      total_gwh = as.numeric(total_gwh),
      fncer_no_hidro_gwh = as.numeric(fncer_no_hidro_gwh),
      pct_renovable = as.numeric(pct_renovable_total),
      pct_fosil = as.numeric(pct_fosil)
    ) %>%
    select(year, meses_reportados, hidro_gwh, biomasa_gwh, solar_gwh,
           eolica_gwh, fosil_gwh, otros_gwh, total_gwh,
           fncer_no_hidro_gwh, pct_renovable, pct_fosil) %>%
    arrange(year)
  
  reportar_calidad(df, "UPME generación anual por grupo")
  df
}, error = function(e) {
  log_msg(paste("UPME anual: error —", e$message), "ERROR")
  NULL
})

# --- C2. UPME — Generación por fuente (mensual) ------------------------------
# Columnas reales: ano, mes, fuente, promedio_mensual_gwh, maximo_gwh,
#                  minimo_gwh, generacion_total_gwh, porcentaje_participacion
upme_fuente_limpio <- tryCatch({
  df <- upme_generacion_fuente %>% limpiar_nombres()
  
  df <- df %>%
    mutate(
      year = as.integer(ano),
      mes = as.integer(mes),
      fuente = str_to_title(fuente),
      promedio_mensual_gwh = as.numeric(promedio_mensual_gwh),
      maximo_gwh = as.numeric(maximo_gwh),
      minimo_gwh = as.numeric(minimo_gwh),
      generacion_total_gwh = as.numeric(generacion_total_gwh),
      porcentaje_participacion = as.numeric(porcentaje_participacion)
    ) %>%
    # Agregar a nivel anual por fuente
    group_by(year, fuente) %>%
    summarise(
      generacion_anual_gwh = sum(generacion_total_gwh, na.rm = TRUE),
      promedio_mensual_gwh = mean(promedio_mensual_gwh, na.rm = TRUE),
      maximo_gwh = max(maximo_gwh, na.rm = TRUE),
      minimo_gwh = min(minimo_gwh, na.rm = TRUE),
      porcentaje_participacion = mean(porcentaje_participacion, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year, fuente)
  
  reportar_calidad(df, "UPME generación por fuente (anual)")
  df
}, error = function(e) {
  log_msg(paste("UPME fuente: error —", e$message), "ERROR")
  NULL
})

# --- C3. INGEI — Series de energía -------------------------------------------
# Columnas reales: anio, energia_total_ktco2eq, electricidad_publica_ktco2eq,
#                  transporte_terrestre_ktco2eq, industrias_manuf_ktco2eq,
#                  fugitivas_ktco2eq, generacion_upme_gwh,
#                  fe_implicito_ingei_sobre_generacion_t_mwh
ingei_energia_limpio <- tryCatch({
  df <- ingei_energia %>% limpiar_nombres()
  
  df <- df %>%
    mutate(
      year = as.integer(anio),
      energia_total_ktco2eq = as.numeric(energia_total_ktco2eq),
      electricidad_publica_ktco2eq = as.numeric(electricidad_publica_ktco2eq),
      transporte_terrestre_ktco2eq = as.numeric(transporte_terrestre_ktco2eq),
      industrias_manuf_ktco2eq = as.numeric(industrias_manuf_ktco2eq),
      fugitivas_ktco2eq = as.numeric(fugitivas_ktco2eq),
      generacion_upme_gwh = as.numeric(generacion_upme_gwh),
      fe_implicito = as.numeric(fe_implicito_ingei_sobre_generacion_t_mwh)
    ) %>%
    select(year, energia_total_ktco2eq, electricidad_publica_ktco2eq,
           transporte_terrestre_ktco2eq, industrias_manuf_ktco2eq,
           fugitivas_ktco2eq, generacion_upme_gwh, fe_implicito) %>%
    arrange(year)
  
  reportar_calidad(df, "INGEI energía")
  df
}, error = function(e) {
  log_msg(paste("INGEI: error —", e$message), "ERROR")
  NULL
})

# --- C4. MinEnergía — Meta FNCER (nivel proyecto) ----------------------------
# Columnas reales: proyecto, tipo, capacidad, departamento, municipio,
#                  codigo_departamento, codigo_municipio, fecha_estimada_fpo,
#                  energia_kwh_d_a, usuarios, inversi_n_estimada_cop,
#                  empleos_estimados, emisiones_co2_ton_a_o
minenergia_meta_limpio <- tryCatch({
  df <- minenergia_meta %>% limpiar_nombres()
  
  df <- df %>%
    mutate(
      codigo_municipio = normalizar_dane(codigo_municipio),
      codigo_departamento = str_pad(as.character(codigo_departamento),
                                    width = 2, side = "left", pad = "0"),
      capacidad_mw = as.numeric(capacidad) / 1000,  # asumir kW → MW
      energia_kwh_d_a = as.numeric(energia_kwh_d_a),
      energia_anual_mwh = energia_kwh_d_a * 365 / 1000,
      inversion_cop = as.numeric(inversi_n_estimada_cop),
      empleos = as.integer(empleos_estimados),
      emisiones_co2_ton = as.numeric(emisiones_co2_ton_a_o),
      tipo = str_to_title(tipo),
      year = if ("fecha_estimada_fpo" %in% names(df)) {
        year(fecha_estimada_fpo)
      } else {
        NA_integer_
      }
    ) %>%
    select(proyecto, tipo, codigo_departamento, codigo_municipio,
           departamento, municipio, year, capacidad_mw, energia_anual_mwh,
           inversion_cop, empleos, emisiones_co2_ton)
  
  reportar_calidad(df, "MinEnergía meta FNCER (proyectos)")
  df
}, error = function(e) {
  log_msg(paste("MinEnergía: error —", e$message), "ERROR")
  NULL
})

# --- C5. Agregar FNCER por departamento y año --------------------------------
meta_fncer_depto <- tryCatch({
  df <- minenergia_meta_limpio %>%
    filter(!is.na(codigo_departamento), !is.na(year)) %>%
    group_by(year, codigo_departamento) %>%
    summarise(
      n_proyectos = n(),
      capacidad_total_mw = sum(capacidad_mw, na.rm = TRUE),
      energia_total_mwh = sum(energia_anual_mwh, na.rm = TRUE),
      inversion_total_cop = sum(inversion_cop, na.rm = TRUE),
      empleos_total = sum(empleos, na.rm = TRUE),
      emisiones_evitadas_ton = sum(emisiones_co2_ton, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year, codigo_departamento)
  
  reportar_calidad(df, "FNCER agregado por departamento")
  df
}, error = function(e) {
  log_msg(paste("FNCER depto: error —", e$message), "ERROR")
  NULL
})

# --- C6. Agregar FNCER por tecnología ----------------------------------------
meta_fncer_tecnologia_limpio <- tryCatch({
  df <- minenergia_meta_limpio %>%
    filter(!is.na(year)) %>%
    group_by(year, tipo) %>%
    summarise(
      n_proyectos = n(),
      capacidad_mw = sum(capacidad_mw, na.rm = TRUE),
      energia_mwh = sum(energia_anual_mwh, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year, tipo)
  
  reportar_calidad(df, "FNCER por tecnología")
  df
}, error = function(e) {
  log_msg(paste("FNCER tecnología: error —", e$message), "ERROR")
  NULL
})

# ============================================================================
# MÓDULO D — PREPARACIÓN DE CLIMATE TRACE COMO PUNTOS SF
# ============================================================================
log_msg("MÓDULO D: Preparación de Climate TRACE como puntos sf", "INFO")

# --- D1. Consolidar activos por sector y año ---------------------------------
if (!is.null(climate_trace_col)) {
  ct_por_sector <- climate_trace_col %>%
    group_by(year, sector) %>%
    summarise(
      n_activos = n(),
      emisiones_sector_t = sum(emissionsQuantity, na.rm = TRUE),
      capacidad_sector = sum(capacity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = sector,
                values_from = c(n_activos, emisiones_sector_t, capacidad_sector),
                names_sep = "_",
                values_fill = 0) %>%
    mutate(
      emisiones_total_t = rowSums(select(., starts_with("emisiones_sector_t_")),
                                  na.rm = TRUE),
      n_activos_total = rowSums(select(., starts_with("n_activos_")),
                                na.rm = TRUE)
    )
  
  log_msg(sprintf("D1: Climate TRACE consolidado: %d años × %d sectores",
                  n_distinct(ct_por_sector$year),
                  n_distinct(climate_trace_limpio$tib$sector)), "OK")
} else {
  ct_por_sector <- NULL
  log_msg("D1: Climate TRACE no disponible", "WARN")
}

# --- D2. Resumen nacional de emisiones por año -------------------------------
if (!is.null(climate_trace_col)) {
  ct_nacional <- climate_trace_col %>%
    group_by(year) %>%
    summarise(
      n_activos_total = n(),
      emisiones_total_t = sum(emissionsQuantity, na.rm = TRUE),
      capacidad_total = sum(capacity, na.rm = TRUE),
      .groups = "drop"
    )
  
  log_msg(sprintf("D2: Climate TRACE nacional: %d años", nrow(ct_nacional)), "OK")
} else {
  ct_nacional <- NULL
}

# ============================================================================
# MÓDULO E — CONSOLIDACIÓN DE COVARIABLES NACIONALES
# ============================================================================

# --- E1. Crear grid de años completo -----------------------------------------
years_panel <- 2010:2025

covariables_nacionales <- tibble(year = years_panel)

# Helper: hacer left_join solo si el dataset no es NULL
join_si_existe <- function(df_base, df_a_unir, by = "year", nombre = "dataset") {
  if (is.null(df_a_unir)) {
    log_msg(paste("E1: Omitiendo", nombre, "(NULL)"), "WARN")
    return(df_base)
  }
  if (!"year" %in% names(df_a_unir)) {
    log_msg(paste("E1: Omitiendo", nombre, "(no tiene columna 'year')"), "WARN")
    return(df_base)
  }
  log_msg(paste("E1: Uniendo", nombre), "INFO")
  left_join(df_base, df_a_unir, by = by)
}

covariables_nacionales <- covariables_nacionales %>%
  join_si_existe(
    edgar_co2_limpio %>% select(year, co2_total_t),
    nombre = "EDGAR CO2"
  ) %>%
  join_si_existe(
    owid_col_limpio %>% select(year, co2, co2_per_capita, co2_per_gdp,
                               energy_per_capita, energy_per_gdp,
                               primary_energy_consumption, total_ghg,
                               population, gdp),
    nombre = "OWID Colombia"
  ) %>%
  join_si_existe(
    ember_indicadores_limpio %>%
      filter(fuente == "Ember") %>%
      select(year, generacion_twh, capacidad_gw, emisiones_mtco2e),
    nombre = "Ember indicadores"
  ) %>%
  join_si_existe(
    upme_anual_limpio %>%
      select(year, total_gwh, fncer_no_hidro_gwh,
             pct_renovable, pct_fosil,
             solar_gwh, eolica_gwh, biomasa_gwh, hidro_gwh, fosil_gwh),
    nombre = "UPME anual"
  ) %>%
  join_si_existe(
    ingei_energia_limpio %>%
      select(year, energia_total_ktco2eq, generacion_upme_gwh, fe_implicito),
    nombre = "INGEI energía"
  ) %>%
  join_si_existe(
    ct_nacional %>% select(year, n_activos_total, emisiones_total_t),
    nombre = "Climate TRACE nacional"
  ) %>%
  join_si_existe(
    ember_evitadas_estimadas %>% select(year, co2_evitado_mt),
    nombre = "Ember evitadas reportadas"
  )

log_msg(sprintf("E1: Covariables nacionales: %d años × %d cols",
                nrow(covariables_nacionales), ncol(covariables_nacionales)), "OK")

# ============================================================================
# MÓDULO F — CLASIFICACIÓN DE TRATAMIENTOS NACIONALES POR AÑO
# ============================================================================
log_msg("MÓDULO F: Clasificando tratamientos nacionales por año", "INFO")

# --- F1. Determinar fuente dominante por año ---------------------------------
# Criterio: la fuente con mayor generación determina el tratamiento del año
asignar_tratamiento_anio <- function(solar, eolica, hidro, biomasa, fosil) {
  caps <- c(
    "T1" = solar,
    "T2" = eolica,
    "T3" = hidro,
    "T4" = biomasa,
    "T0" = fosil
  )
  caps <- caps[!is.na(caps)]
  if (length(caps) == 0 || sum(caps) == 0) return("T0")
  names(caps)[which.max(caps)]
}

panel_nacional <- covariables_nacionales %>%
  mutate(
    tratamiento = mapply(
      asignar_tratamiento_anio,
      solar = solar_gwh,
      eolica = eolica_gwh,
      hidro = hidro_gwh,
      biomasa = biomasa_gwh,
      fosil = fosil_gwh
    ),
    tratamiento = factor(tratamiento, levels = c("T0", "T1", "T2", "T3", "T4")),
    nombre_tratamiento = factor(
      case_when(
        tratamiento == "T0" ~ "Control (fósil)",
        tratamiento == "T1" ~ "Solar FV",
        tratamiento == "T2" ~ "Eólica",
        tratamiento == "T3" ~ "Pequeña Hidro",
        tratamiento == "T4" ~ "Biomasa/Biogás"
      ),
      levels = TRATAMIENTOS$nombre
    )
  )

# --- F2. Verificar distribución de tratamientos ------------------------------
dist_trat <- panel_nacional %>%
  count(tratamiento, name = "n_anios") %>%
  mutate(pct = round(100 * n_anios / sum(n_anios), 1))

log_msg("Distribución de tratamientos por año:", "INFO")
print(as.data.frame(dist_trat))

# ============================================================================
# MÓDULO G — VARIABLES DERIVADAS Y TRANSFORMACIONES
# ============================================================================
log_msg("MÓDULO G: Variables derivadas", "INFO")

panel_nacional <- panel_nacional %>%
  mutate(
    # Logaritmo de emisiones (para modelos lineales)
    log_co2_edgar = log1p(co2_total_t),
    log_co2_owid = log1p(co2 * 1e6),  # OWID en Mt → toneladas
    
    # Intensidad de carbono
    intensidad_carbono = co2_total_t / pmax(generacion_upme_gwh * 1e3, 1),  # t CO2 / MWh
    
    # Cambio porcentual año a año
    co2_lag1 = dplyr::lag(co2_total_t, 1),
    co2_pct_change = (co2_total_t - co2_lag1) / pmax(co2_lag1, 1) * 100,
    
    # Indicador binario de transición energética
    es_renovable = ifelse(tratamiento != "T0", 1L, 0L),
    
    # Año desde el inicio del panel
    anio_rel = year - min(year, na.rm = TRUE),
    
    # Participación de renovables no hidro
    pct_fncer_no_hidro = fncer_no_hidro_gwh / pmax(total_gwh, 1) * 100
  )

log_msg(sprintf("G: Panel nacional con %d columnas derivadas", ncol(panel_nacional)), "OK")

# ============================================================================
# MÓDULO H — VALIDACIÓN FINAL
# ============================================================================
log_msg("MÓDULO H: Validación final del panel", "INFO")

validaciones <- list(
  n_anios = n_distinct(panel_nacional$year),
  n_filas = nrow(panel_nacional),
  pct_na_co2_edgar = round(100 * mean(is.na(panel_nacional$co2_total_t)), 2),
  pct_na_co2_owid = round(100 * mean(is.na(panel_nacional$co2)), 2),
  pct_na_trat = round(100 * mean(is.na(panel_nacional$tratamiento)), 2),
  duplicados = nrow(panel_nacional) - n_distinct(panel_nacional$year)
)

for (v in names(validaciones)) {
  log_msg(sprintf("H1: %s = %s", v, validaciones[[v]]), "INFO")
}

if (validaciones$duplicados > 0) {
  log_msg("H1: DUPLICADOS detectados — revisando", "WARN")
  panel_nacional <- panel_nacional %>%
    distinct(year, .keep_all = TRUE)
}

# ============================================================================
# MÓDULO I — EXPORTACIÓN DE DATOS PROCESADOS
# ============================================================================
log_msg("MÓDULO I: Exportando datos procesados", "INFO")

exportar <- function(obj, nombre) {
  if (is.null(obj)) return(invisible(NULL))
  path <- file.path(PATHS$datos, paste0(nombre, ".rds"))
  guardar_rds(obj, path)
}

exportar(edgar_co2_limpio, "edgar_co2_nacional")
exportar(if (!is.null(climate_trace_limpio)) climate_trace_limpio else NULL,
         "climate_trace_colombia")
exportar(owid_col_limpio, "owid_colombia")
exportar(ember_indicadores_limpio, "ember_indicadores_colombia")
exportar(ember_evitadas_reportadas, "ember_evitadas_reportadas")
exportar(ember_evitadas_estimadas, "ember_evitadas_estimadas")
exportar(upme_anual_limpio, "upme_generacion_grupo")
exportar(upme_fuente_limpio, "upme_generacion_fuente")
exportar(ingei_energia_limpio, "ingei_energia_serie")
exportar(minenergia_meta_limpio, "minenergia_meta_fncer_proyectos")
exportar(meta_fncer_tecnologia_limpio, "meta_fncer_tecnologia")
exportar(covariables_nacionales, "covariables_nacionales")
exportar(panel_nacional, "panel_nacional_completo")

# --- Guardar Climate TRACE como sf -------------------------------------------
if (!is.null(climate_trace_limpio)) {
  path_ct <- file.path(PATHS$datos_proc, "climate_trace_colombia_sf.gpkg")
  st_write(climate_trace_limpio$sf, path_ct, driver = "GPKG", quiet = TRUE)
  log_msg(sprintf("I2: Climate TRACE sf guardado: %d activos",
                  nrow(climate_trace_limpio$sf)), "OK")
}

# ============================================================================
# MÓDULO J — RESUMEN FINAL Y LOG
# ============================================================================
tiempo_total <- round(as.numeric(difftime(Sys.time(), inicio, units = "mins")), 2)

message("
", paste(rep("=", 72), collapse = ""))
message("PREPROCESSING COMPLETADO")
message(paste(rep("=", 72), collapse = ""))
message("  Tiempo total        : ", tiempo_total, " min")
message("  Años en panel       : ", validaciones$n_anios)
message("  Filas del panel     : ", validaciones$n_filas)
message("  Tratamientos        : ",
        paste(levels(panel_nacional$tratamiento), collapse = ", "))
message("  NA en CO2 EDGAR     : ", validaciones$pct_na_co2_edgar, "%")
message("  NA en CO2 OWID      : ", validaciones$pct_na_co2_owid, "%")
message("  NA en tratamiento   : ", validaciones$pct_na_trat, "%")
message("  Duplicados          : ", validaciones$duplicados)
message(paste(rep("=", 72), collapse = ""))

log_path <- file.path(PATHS$datos_proc, "preprocessing_log.txt")
writeLines(c(
  paste("Fecha:", Sys.time()),
  paste("Tiempo:", tiempo_total, "min"),
  paste("Años:", validaciones$n_anios),
  paste("Filas panel:", validaciones$n_filas),
  paste("NA CO2 EDGAR%:", validaciones$pct_na_co2_edgar),
  paste("NA CO2 OWID%:", validaciones$pct_na_co2_owid),
  paste("NA Trat%:", validaciones$pct_na_trat),
  paste("Duplicados:", validaciones$duplicados),
  "",
  "Datasets exportados:",
  list.files(PATHS$datos, pattern = "\\.rds$|\\.gpkg$")
), log_path)

message("Log guardado en: ", log_path)

message("
Siguiente paso sugerido:")
message("   source('scripts/02_spatial_sampling.R')  # Diseño muestral GRTS")

