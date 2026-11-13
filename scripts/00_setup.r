#!/usr/bin/env Rscript
# ============================================================================
# 00_setup.R — Configuración inicial del proyecto
# ============================================================================
# Proyecto : Análisis Espacial de CO₂ y Evaluación Cuasi-Experimental de FNCER
#            en Colombia
# Autor    : Grupo de Investigación en Estadística Aplicada
# Fecha    : 2026-11-14
# Propósito:
#   1. Definir la estructura de directorios del proyecto
#   2. Verificar y activar el entorno reproducible (renv)
#   3. Cargar paquetes organizados por asignatura metodológica
#   4. Definir constantes globales (URLs, paths, versiones de datos)
#   5. Configurar opciones de R y de paquetes espaciales
#   6. Definir funciones utilitarias reutilizadas en todo el pipeline
#   7. Verificar la disponibilidad de Python para scripts de cosecha
# Uso      : source("scripts/00_setup.R")
# ============================================================================

# ----------------------------------------------------------------------------
# 0. Mensaje de bienvenida
# ----------------------------------------------------------------------------
message("\n", paste(rep("=", 72), collapse = ""))
message("🌎  PROYECTO CO₂ COLOMBIA — Configuración inicial")
message(paste(rep("=", 72), collapse = ""), "\n")

# ----------------------------------------------------------------------------
# 1. Detectar la raíz del proyecto (funciona en RStudio, Rscript y Quarto)
# ----------------------------------------------------------------------------
detectar_raiz <- function() {
  # Prioridad: variable de entorno > RStudio project > archivo actual
  if (nzchar(Sys.getenv("CO2_PROJECT_ROOT"))) {
    return(normalizePath(Sys.getenv("CO2_PROJECT_ROOT")))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() &&
      !is.null(rstudioapi::getActiveProject())) {
    return(normalizePath(rstudioapi::getActiveProject()))
  }
  # Fallback: subir desde el script hasta encontrar README.md
  p <- normalizePath(dirname(sys.frame(1)$ofile %||% "."))
  while (p != dirname(p)) {
    if (file.exists(file.path(p, "README.md"))) return(p)
    p <- dirname(p)
  }
  stop("No se encontró la raíz del proyecto (falta README.md). ",
       "Define la variable de entorno CO2_PROJECT_ROOT.")
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

RAIZ_PROYECTO <- detectar_raiz()
Sys.setenv(CO2_PROJECT_ROOT = RAIZ_PROYECTO)
message("📁 Raíz del proyecto : ", RAIZ_PROYECTO)

# ----------------------------------------------------------------------------
# 2. Crear estructura de directorios si no existe
# ----------------------------------------------------------------------------
dirs <- c(
  "datos/raw/edgar",
  "datos/raw/climate_trace",
  "datos/raw/owid",
  "datos/raw/gadm",
  "datos/raw/datos_gov",
  "datos/processed",
  "scripts",
  "output/tables",
  "output/maps",
  "output/figures",
  "docs",
  "manuscript",
  "obsidian-vault/01_Protocolo",
  "obsidian-vault/02_Diario",
  "obsidian-vault/03_Literatura",
  "obsidian-vault/04_Glosario"
)

dir.create(file.path(RAIZ_PROYECTO, dirs),
           recursive = TRUE, showWarnings = FALSE)
message("✅ Estructura de directorios verificada (", length(dirs), " carpetas)")

# ----------------------------------------------------------------------------
# 3. Verificar y activar renv (reproducibilidad)
# ----------------------------------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  message("⚠️  renv no está instalado. Instalando...")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

if (file.exists(file.path(RAIZ_PROYECTO, "renv.lock"))) {
  renv::activate(project = RAIZ_PROYECTO, quiet = TRUE)
  message("🔒 Entorno renv activado. Usa renv::restore() si faltan paquetes.")
} else {
  message("⚠️  renv.lock no encontrado. Ejecuta renv::snapshot() al finalizar.")
}

# ----------------------------------------------------------------------------
# 4. Función auxiliar para cargar paquetes (instala si es necesario)
# ----------------------------------------------------------------------------
cargar_paquete <- function(pkg, repos = "https://cloud.r-project.org") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("  📦 Instalando ", pkg, "...")
    install.packages(pkg, repos = repos, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  return(invisible(TRUE))
}

# Paquetes de Bioconductor (instalación diferente)
cargar_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  return(invisible(TRUE))
}

# ----------------------------------------------------------------------------
# 5. Cargar paquetes organizados por asignatura metodológica
# ----------------------------------------------------------------------------
message("\n📚 Cargando paquetes por asignatura...")

# --- Base y utilidades ------------------------------------------------------
message("\n  🔧 Base y utilidades:")
pkgs_base <- c("tidyverse", "data.table", "jsonlite", "here", "fs",
               "glue", "cli", "purrr", "stringr", "lubridate")
walk(pkgs_base, cargar_paquete)

# --- 🗺️ Estadística Espacial -----------------------------------------------
message("\n  🗺️  Estadística Espacial:")
pkgs_espacial <- c("sf", "terra", "stars", "spdep", "sfdep", "gstat",
                   "spatialreg", "gwrr", "exactextractr", "tmap",
                   "leaflet", "mapview")
walk(pkgs_espacial, cargar_paquete)

# --- 🎲 Teoría del Muestreo -------------------------------------------------
message("\n  🎲 Teoría del Muestreo:")
pkgs_muestreo <- c("spsurvey", "BalancedSampling", "survey", "sampling",
                   "boot")
walk(pkgs_muestreo, cargar_paquete)

# --- 🤖 Machine Learning ----------------------------------------------------
message("\n  🤖 Machine Learning:")
pkgs_ml <- c("grf", "xgboost", "lightgbm", "randomForest", "missForest",
             "blockCV", "shapviz", "vip", "tidymodels", "rsample",
             "yardstick", "DoubleML")
walk(pkgs_ml, cargar_paquete)

# --- 🧪 Diseño de Experimentos / Inferencia Causal --------------------------
message("\n  🧪 Diseño de Experimentos:")
pkgs_exp <- c("MatchIt", "cobalt", "lme4", "lmerTest", "did",
              "rdrobust", "AER", "ivreg", "dagitty", "ggdag",
              "fixest", "sandwich", "lmtest")
walk(pkgs_exp, cargar_paquete)

# --- 📊 Visualización y Reporte ---------------------------------------------
message("\n  📊 Visualización y Reporte:")
pkgs_viz <- c("ggplot2", "patchwork", "scales", "viridis", "RColorBrewer",
              "kableExtra", "gt", "knitr", "quarto")
walk(pkgs_viz, cargar_paquete)

# ----------------------------------------------------------------------------
# 6. Configuración de opciones globales
# ----------------------------------------------------------------------------
message("\n⚙️  Configurando opciones globales...")

# Opciones de R
options(
  scipen        = 999,          # Evitar notación científica
  digits        = 4,
  timeout       = 600,          # Timeout para descargas (segundos)
  mc.cores      = parallel::detectCores() - 1,  # Paralelización
  future.plan   = "multisession"
)

# Opciones de sf (GDAL)
sf::sf_use_s2(TRUE)             # Geometría esférica (más precisa)
Sys.setenv(GDAL_MAX_DATASET_POOL_SIZE = "200")

# Opciones de terra
terra::terraOptions(verbosity = 0)

# Tema ggplot2 consistente
theme_set(theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold"),
                  strip.text   = element_text(face = "bold")))

# Paleta de colores del proyecto (accesible, daltonismo-friendly)
PALETA_PROYECTO <- c(
  "fossil"    = "#d73027",  # Rojo — fósil
  "solar"     = "#fee08b",  # Amarillo — solar
  "eolica"    = "#66c2a5",  # Verde-azul — eólica
  "hidro"     = "#3288bd",  # Azul — hidro
  "biomasa"   = "#8c510a",  # Marrón — biomasa
  "control"   = "#999999"   # Gris — control
)

message("✅ Opciones globales configuradas")

# ----------------------------------------------------------------------------
# 7. Constantes del proyecto
# ----------------------------------------------------------------------------
message("\n📌 Definiendo constantes del proyecto...")

# --- Paths ------------------------------------------------------------------
PATHS <- list(
  raiz         = RAIZ_PROYECTO,
  datos_raw    = file.path(RAIZ_PROYECTO, "datos/raw"),
  datos_proc   = file.path(RAIZ_PROYECTO, "datos/processed"),
  output       = file.path(RAIZ_PROYECTO, "output"),
  scripts      = file.path(RAIZ_PROYECTO, "scripts"),
  manuscript   = file.path(RAIZ_PROYECTO, "manuscript"),
  catalogo_gov = file.path(RAIZ_PROYECTO, "datos/catalogo_datos_gov.csv"),
  referencias  = file.path(RAIZ_PROYECTO, "datos/referencias_verificadas.json")
)

# --- URLs de fuentes globales -----------------------------------------------
URLS <- list(
  # EDGAR v2026
  edgar_base   = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/EDGAR_2026_GHG/",
  edgar_01km   = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/EDGAR_2026_GHG/CO2/TOTALS/emi_nc/",
  edgar_1km    = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/EDGAR_2026_GHG/EDGAR_2025_1km/",
  
  # OWID
  owid_co2     = "https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv",
  
  # Climate TRACE (Zenodo)
  climate_trace_zenodo = "https://doi.org/10.5281/zenodo.20830133",
  climate_trace_api    = "https://api.climatetrace.org/v7/",
  
  # GADM
  gadm_base    = "https://geodata.ucdavis.edu/gadm/gadm4.1/json/",
  gadm_col_0   = "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_COL_0.json",
  gadm_col_1   = "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_COL_1.json",
  gadm_col_2   = "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_COL_2.json",
  
  # Natural Earth
  ne_10m       = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip",
  
  # datos.gov.co (Socrata API)
  datos_gov_api = "https://www.datos.gov.co/resource/"
)

# --- Parámetros del diseño muestral -----------------------------------------
DISENO_MUESTRAL <- list(
  resolucion_malla   = 5000,       # metros (5 km × 5 km)
  n_estratos_region  = 5,          # Andina, Caribe, Pacífica, Orinoquía, Amazonía
  n_estratos_upme    = 3,          # Alto, Medio, Bajo desarrollo energético
  confianza          = 0.95,
  margen_error       = 0.05,
  muestra_objetivo   = 1750,       # celdas
  crs_proyectado     = "EPSG:9377" # MAGNA-SIRGAS / Colombia Bogota grid
)

# --- Tratamientos energéticos -----------------------------------------------
TRATAMIENTOS <- tibble::tribble(
  ~codigo, ~nombre,           ~descripcion,
  "T0",    "Control (fósil)", "Matriz energética fósil predominante (gas/carbón)",
  "T1",    "Solar FV",        "Energía solar fotovoltaica",
  "T2",    "Eólica",          "Energía eólica terrestre",
  "T3",    "Pequeña Hidro",   "Hidroeléctrica de pequeña escala (<10 MW)",
  "T4",    "Biomasa/Biogás",  "Residuos agrícolas y urbanos"
)

# --- Regiones naturales de Colombia -----------------------------------------
REGIONES <- c("Andina", "Caribe", "Pacífica", "Orinoquía", "Amazonía")

# --- CRS (Sistemas de Referencia de Coordenadas) ----------------------------
CRS <- list(
  wgs84       = 4326,   # WGS84 (lat/lon) — intercambio
  magna_orig  = 3116,   # MAGNA-SIRGAS origen
  magna_bog   = 9377,   # MAGNA-SIRGAS Bogotá grid (proyectado, metros)
  utm_18n     = 32618   # UTM 18N (backup)
)

message("✅ Constantes definidas")

# ----------------------------------------------------------------------------
# 8. Funciones utilitarias reutilizadas en el pipeline
# ----------------------------------------------------------------------------
message("\n🛠️  Definiendo funciones utilitarias...")

# --- Guardar/cargar objetos de R de forma eficiente -------------------------
guardar_rds <- function(x, path, compress = "gz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = path, compress = compress)
  message("  💾 Guardado: ", path, " (", 
          format(file.size(path) / 1e6, digits = 2), " MB)")
  invisible(x)
}

cargar_rds <- function(path) {
  if (!file.exists(path)) {
    stop("Archivo no encontrado: ", path)
  }
  readRDS(path)
}

# --- Descarga con reintentos y cache ----------------------------------------
descargar_archivo <- function(url, destino, reintentos = 3, pausa = 5) {
  if (file.exists(destino)) {
    message("  ⏩ Ya existe: ", basename(destino))
    return(destino)
  }
  dir.create(dirname(destino), recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(reintentos)) {
    tryCatch({
      download.file(url, destino, mode = "wb", quiet = TRUE)
      message("  ✅ Descargado: ", basename(destino))
      return(destino)
    }, error = function(e) {
      message("  ⚠️  Intento ", i, "/", reintentos, " falló: ", e$message)
      if (i < reintentos) Sys.sleep(pausa)
    })
  }
  stop("No se pudo descargar: ", url)
}

# --- Validar que un data.frame tiene columnas esperadas ---------------------
validar_columnas <- function(df, cols_esperadas, nombre = "data.frame") {
  faltantes <- setdiff(cols_esperadas, names(df))
  if (length(faltantes) > 0) {
    stop("Columnas faltantes en ", nombre, ": ",
         paste(faltantes, collapse = ", "))
  }
  invisible(TRUE)
}

# --- Log con timestamp ------------------------------------------------------
log_msg <- function(msg, nivel = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  icon <- switch(nivel,
                 "INFO"  = "ℹ️",
                 "OK"    = "✅",
                 "WARN"  = "⚠️",
                 "ERROR" = "❌",
                 "ℹ️")
  message(paste0(icon, " [", timestamp, "] ", msg))
}

# --- Conversión de unidades de CO₂ ------------------------------------------
# EDGAR reporta en toneladas de CO₂; algunos modelos usan toneladas de C
tC_a_tCO2 <- function(tC) tC * 44 / 12  # = tC * 3.667
tCO2_a_tC <- function(tCO2) tCO2 * 12 / 44

# --- Resumen rápido de objeto sf --------------------------------------------
resumen_sf <- function(obj) {
  stopifnot(inherits(obj, "sf"))
  list(
    n_filas    = nrow(obj),
    n_columnas = ncol(obj),
    crs        = sf::st_crs(obj)$input,
    bbox       = sf::st_bbox(obj),
    tipo_geom  = unique(sf::st_geometry_type(obj))
  )
}

message("✅ Funciones utilitarias cargadas")

# ----------------------------------------------------------------------------
# 9. Verificar disponibilidad de Python (para scripts de cosecha)
# ----------------------------------------------------------------------------
message("\n🐍 Verificando Python...")
python_path <- Sys.which("python3")
if (nzchar(python_path)) {
  py_version <- system2(python_path, "--version", stdout = TRUE, stderr = TRUE)
  message("  ✅ Python disponible: ", py_version, " (", python_path, ")")
  PYTHON_OK <- TRUE
} else {
  python_path <- Sys.which("python")
  if (nzchar(python_path)) {
    py_version <- system2(python_path, "--version", stdout = TRUE, stderr = TRUE)
    message("  ⚠️  Python (no python3) disponible: ", py_version)
    PYTHON_OK <- TRUE
  } else {
    message("  ❌ Python no encontrado. Los scripts 01_catalogo_datos_gov.py ",
            "y 05_referencias.py no podrán ejecutarse.")
    PYTHON_OK <- FALSE
  }
}

# ----------------------------------------------------------------------------
# 10. Resumen final de inicio
# ----------------------------------------------------------------------------
message("\n", paste(rep("-", 72), collapse = ""))
message("📋 RESUMEN DE CONFIGURACIÓN")
message(paste(rep("-", 72), collapse = ""))
message("  📁 Raíz del proyecto  : ", RAIZ_PROYECTO)
message("  📦 Paquetes cargados  : ",
        length(c(pkgs_base, pkgs_espacial, pkgs_muestreo,
                 pkgs_ml, pkgs_exp, pkgs_viz)))
message("  🎯 Muestra objetivo   : ", DISENO_MUESTRAL$muestra_objetivo, " celdas")
message("  🗺️  Resolución malla   : ", DISENO_MUESTRAL$resolucion_malla, " m")
message("  🧪 Tratamientos       : ", nrow(TRATAMIENTOS), " (T0–T4)")
message("  🌎 Regiones           : ", length(REGIONES))
message("  🐍 Python disponible  : ", PYTHON_OK)
message("  🔢 CRS proyectado     : EPSG:", CRS$magna_bog, " (MAGNA-SIRGAS Bogotá)")
message(paste(rep("-", 72), collapse = ""))
message("✅ SETUP COMPLETADO — Listo para ejecutar el pipeline\n")

# ----------------------------------------------------------------------------
# 11. Exportar objetos al entorno global (para que otros scripts los usen)
# ----------------------------------------------------------------------------
# Los objetos RAIZ_PROYECTO, PATHS, URLS, DISENO_MUESTRAL, TRATAMIENTOS,
# REGIONES, CRS, PALETA_PROYECTO y las funciones utilitarias quedan
# disponibles en el entorno global al hacer source() de este script.

# Guardar un "snapshot" del setup para referencia
SETUP_INFO <- list(
  fecha_setup     = Sys.time(),
  r_version       = R.version.string,
  raiz_proyecto   = RAIZ_PROYECTO,
  paquetes_cargados = c(pkgs_base, pkgs_espacial, pkgs_muestreo,
                        pkgs_ml, pkgs_exp, pkgs_viz),
  python_ok       = PYTHON_OK,
  diseño_muestral = DISENO_MUESTRAL
)

message("💡 Siguiente paso sugerido:")
message("   source('scripts/01b_data_acquisition.R')  # Descargar fuentes globales")
message("   O bien: python3 scripts/01_catalogo_datos_gov.py  # Cosechar datos.gov.co\n")