# ============================================================================
# 00_setup.R — Configuración del proyecto (RUTA, PAQUETES, FUNCIONES)
# Proyecto: Emisiones de CO2 y transición energética en Colombia 2015-2030
# ============================================================================
# Este script NO descarga datos ni calcula resultados: solo prepara el entorno.
# Se ejecuta con:  Rscript scripts/00_setup.R
# ============================================================================

# --- 1. Raíz del proyecto ---------------------------------------------------
detectar_raiz <- function() {
  if (nzchar(Sys.getenv("CO2_ROOT"))) return(normalizePath(Sys.getenv("CO2_ROOT")))
  p <- normalizePath(getwd(), winslash = "/")
  for (i in 1:6) {
    if (file.exists(file.path(p, "co2-emisions-col.Rproj")) ||
        (file.exists(file.path(p, "README.md")) && dir.exists(file.path(p, "scripts")))) {
      return(p)
    }
    p <- dirname(p)
  }
  stop("No se encontró la raíz del proyecto. Ejecuta desde la carpeta del proyecto.")
}
RAIZ <- detectar_raiz()
Sys.setenv(CO2_ROOT = RAIZ)

# --- 2. Rutas ---------------------------------------------------------------
RUTAS <- list(
  datos      = file.path(RAIZ, "datos"),
  brutos     = file.path(RAIZ, "datos", "brutos"),
  geometrias = file.path(RAIZ, "datos", "geometrias"),
  procesados = file.path(RAIZ, "datos", "procesados"),
  salida     = file.path(RAIZ, "output"),
  figuras    = file.path(RAIZ, "output", "figures"),
  mapas      = file.path(RAIZ, "output", "maps"),
  tablas     = file.path(RAIZ, "output", "tables"),
  manuscrito = file.path(RAIZ, "manuscript")
)
invisible(lapply(RUTAS, dir.create, recursive = TRUE, showWarnings = FALSE))

# --- 3. Paquetes por asignatura --------------------------------------------
paquetes <- list(
  base = c("dplyr", "tidyr", "purrr", "tibble", "readr", "stringr", "forcats",
           "ggplot2", "scales", "patchwork", "viridis", "gridExtra"),
  espacial = c("sf", "terra", "spdep", "gstat", "spsurvey", "spatialreg"),
  doe = c("car"),
  ml = c("randomForest", "ranger", "xgboost", "caret"),
  reporte = c("knitr", "rmarkdown", "jsonlite")
)
faltantes <- character(0)
for (grupo in names(paquetes)) {
  for (p in paquetes[[grupo]]) {
    ok <- suppressWarnings(suppressPackageStartupMessages(
      requireNamespace(p, quietly = TRUE)))
    if (!ok) faltantes <- c(faltantes, p)
  }
}
if (length(faltantes)) {
  warning("Paquetes no disponibles: ", paste(unique(faltantes), collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(sf); library(terra)
})
sf_use_s2(FALSE)
terraOptions(progress = 0, tempdir = tempdir())

# --- 4. Constantes del estudio ---------------------------------------------
CONST <- list(
  anio_base      = 2015,   # línea base observacional (EDGAR)
  anio_ref       = 2023,   # último año EDGAR observado en la grilla
  anio_horizonte = 2030,   # meta Ley 2169 de 2021
  # Ley 2169 de 2021, art. 5 / NDC 2020: -51% de GEI en 2030 frente al
  # escenario de referencia (BAU) de 345,8 MtCO2eq -> techo de 169,44 MtCO2eq.
  bau_2030_mtco2e    = 345.8,
  meta_2030_mtco2e   = 169.44,
  reduccion_ley      = 0.51,
  crs_geo            = 4326,
  crs_proj           = 9377,  # MAGNA-SIRGAS / Origen Nacional
  semilla            = 2026
)
set.seed(CONST$semilla)

# --- 5. Paletas y tema gráfico ---------------------------------------------
PAL_REGION <- c(
  "Caribe" = "#E5862B", "Andina" = "#2E6F9E", "Pacifica" = "#1B9E77",
  "Orinoquia" = "#8C6BB1", "Amazonia" = "#4D9221"
)
PAL_TEC <- c(
  "Solar" = "#F4C430", "Eolica" = "#4A90D9", "Biomasa" = "#7A9E3B",
  "Hidro" = "#1F78B4", "Fosil" = "#5B5B5B"
)
PAL_ESC <- c("BAU" = "#B2182B", "Moderado" = "#EF8A62", "Acelerado" = "#2166AC")

tema <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(colour = "grey30"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.title = element_text(colour = "grey20")
    )
}
theme_set(tema())

# --- 6. Utilidades ----------------------------------------------------------
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                 paste0(...)))
paso <- function(titulo) {
  cat("\n", strrep("=", 72), "\n", titulo, "\n", strrep("=", 72), "\n", sep = "")
}

fmt_num <- function(x, dec = 1) formatC(x, format = "f", big.mark = ",", decimal.mark = ".", digits = dec)
fmt_p   <- function(p) ifelse(p < 0.001, "< 0,001", formatC(p, format = "f", digits = 3, decimal.mark = ","))

guardar_fig <- function(p, nombre, ancho = 9, alto = 6, dir = RUTAS$figuras, dpi = 300) {
  f <- file.path(dir, paste0(nombre, ".png"))
  ggsave(f, p, width = ancho, height = alto, dpi = dpi, bg = "white")
  msg("figura -> ", basename(f))
  invisible(f)
}

guardar_tabla <- function(df, nombre, dir = RUTAS$tablas) {
  f <- file.path(dir, paste0(nombre, ".csv"))
  utils::write.csv(df, f, row.names = FALSE, fileEncoding = "UTF-8")
  msg("tabla  -> ", basename(f))
  invisible(f)
}

# IC 95% de un vector (bootstrap percentil, reproducible)
ic95 <- function(x, B = 2000) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(c(li = NA, ls = NA))
  bs <- replicate(B, mean(sample(x, replace = TRUE)))
  stats::quantile(bs, c(0.025, 0.975), names = FALSE) |> setNames(c("li", "ls"))
}

if (sys.nframe() == 0) {
  paso("ENTORNO LISTO")
  cat("Raíz:", RAIZ, "\n")
  cat("R:", R.version.string, "\n")
  cat("Paquetes faltantes:", ifelse(length(faltantes), paste(faltantes, collapse = ", "), "ninguno"), "\n")
  cat("Salidas:", RUTAS$salida, "\n")
}
