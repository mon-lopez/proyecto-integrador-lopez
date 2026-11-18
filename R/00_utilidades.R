# =============================================================================
# 00_utilidades.R
# Funciones base para trabajar con emisiones de CO2 por pais y por coordenada.
# Dependencias: solo jsonlite (viene con R). No requiere sf ni terra.
# =============================================================================

# Si existe la carpeta .Rlib del proyecto (sf/terra ya instalados), se usa
if (dir.exists(".Rlib")) .libPaths(c(".Rlib", .libPaths()))

suppressPackageStartupMessages(library(jsonlite))

# --- Rutas -------------------------------------------------------------------
# Se asume que este script vive en <proyecto>/R/.
.ruta_proyecto <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) == 1) return(normalizePath(file.path(dirname(f), "..")))
  normalizePath(".")   # fallback: directorio actual
}
RUTA_DATOS <- function(...) file.path(.ruta_proyecto(), "datos", ...)

# =============================================================================
# 1. LEER POLIGONOS DE PAISES (GeoJSON -> lista indexable)
# =============================================================================
# Devuelve una lista con un elemento por pais:
#   iso3, pais, continente, bbox = c(xmin, ymin, xmax, ymax),
#   poligonos = lista de poligonos; cada uno es lista de anillos (1 = exterior, 2+ = huecos)
leer_paises <- function(archivo = RUTA_DATOS("paises_ne_10m.geojson"),
                        campo_iso = "ISO_A3_EH",
                        campo_nombre = "ADMIN") {
  gj <- jsonlite::fromJSON(archivo, simplifyVector = FALSE)
  out <- vector("list", length(gj$features))

  for (i in seq_along(gj$features)) {
    ft   <- gj$features[[i]]
    pr   <- ft$properties
    geom <- ft$geometry
    if (is.null(geom)) next

    # Normalizar Polygon y MultiPolygon a una lista de poligonos
    if (identical(geom$type, "Polygon")) {
      poligonos <- list(geom$coordinates)
    } else if (identical(geom$type, "MultiPolygon")) {
      poligonos <- geom$coordinates
    } else next

    # Convertir el GeoJSON anidado en matrices conservando la jerarquia
    pol_mats <- list()
    for (pol in poligonos) {
      anillos <- list()
      for (anillo in pol) {
        m <- do.call(rbind, lapply(anillo, function(p) c(p[[1]], p[[2]])))
        if (!is.null(m) && nrow(m) >= 3) anillos[[length(anillos) + 1]] <- m
      }
      if (length(anillos)) pol_mats[[length(pol_mats) + 1]] <- anillos
    }
    if (!length(pol_mats)) next

    xy <- do.call(rbind, unlist(pol_mats, recursive = FALSE))
    out[[i]] <- list(
      iso3       = if (is.null(pr[[campo_iso]])) NA_character_ else pr[[campo_iso]],
      pais       = if (is.null(pr[[campo_nombre]])) NA_character_ else pr[[campo_nombre]],
      continente = if (is.null(pr[["CONTINENT"]])) NA_character_ else pr[["CONTINENT"]],
      bbox       = c(xmin = min(xy[, 1]), ymin = min(xy[, 2]),
                     xmax = max(xy[, 1]), ymax = max(xy[, 2])),
      poligonos  = pol_mats
    )
  }
  out <- out[!vapply(out, is.null, logical(1))]
  out
}

# =============================================================================
# 2. PUNTO DENTRO DE POLIGONO (ray casting vectorizado)
# =============================================================================
# x, y: vectores de coordenadas. anillo: matriz de 2 columnas (lon, lat).
# Solo se evaluan los segmentos cuya latitud puede cruzar la semirrecta del
# punto; con lotes de puntos en una misma celda eso reduce el costo enormemente.
.en_anillo <- function(x, y, anillo) {
  n <- nrow(anillo)
  if (n < 3) return(rep(FALSE, length(x)))
  xi <- anillo[, 1]; yi <- anillo[, 2]
  xj <- c(xi[-1], xi[1]); yj <- c(yi[-1], yi[1])

  # Descartar segmentos fuera del rango de latitudes del lote de puntos
  if (length(y) > 1) {
    ymin <- min(y); ymax <- max(y)
  } else {
    ymin <- ymax <- y
  }
  keep <- (pmax(yi, yj) >= ymin) & (pmin(yi, yj) <= ymax)
  if (!any(keep)) return(rep(FALSE, length(x)))
  xi <- xi[keep]; yi <- yi[keep]; xj <- xj[keep]; yj <- yj[keep]

  dentro <- rep(FALSE, length(x))
  for (i in seq_along(xi)) {
    cond <- (yi[i] > y) != (yj[i] > y)          # el segmento cruza la latitud?
    if (!any(cond)) next
    xcruce <- (xj[i] - xi[i]) * (y[cond] - yi[i]) / (yj[i] - yi[i]) + xi[i]
    idx <- which(cond)
    cruza <- idx[x[idx] < xcruce]
    dentro[cruza] <- !dentro[cruza]
  }
  dentro
}

# =============================================================================
# 2b. INDICE ESPACIAL por celdas (para miles de puntos)
# =============================================================================
# Aplana todos los poligonos con su bbox y dice que poligonos pueden tocar
# cada celda de una rejilla gruesa. Evita recorrer los 250 paises por punto.
construir_indice <- function(paises, res = 5) {
  pols <- list()
  for (i in seq_along(paises)) {
    for (pol in paises[[i]]$poligonos) {
      xy <- do.call(rbind, pol)
      pols[[length(pols) + 1]] <- list(
        pais_idx = i,
        anillos  = pol,
        bbox     = c(xmin = min(xy[, 1]), ymin = min(xy[, 2]),
                     xmax = max(xy[, 1]), ymax = max(xy[, 2]))
      )
    }
  }
  ncol <- ceiling(360 / res); nrow <- ceiling(180 / res)
  celdas <- vector("list", ncol * nrow)

  for (k in seq_along(pols)) {
    bb <- pols[[k]]$bbox
    i0 <- max(1,     floor((bb["xmin"] + 180) / res) + 1)
    i1 <- min(ncol,  floor((bb["xmax"] + 180) / res) + 1)
    j0 <- max(1,     floor((bb["ymin"] +  90) / res) + 1)
    j1 <- min(nrow,  floor((bb["ymax"] +  90) / res) + 1)
    for (j in j0:j1) for (i in i0:i1) {
      idx <- (j - 1) * ncol + i
      celdas[[idx]] <- c(celdas[[idx]], k)
    }
  }
  list(poligonos = pols, celdas = celdas, res = res, ncol = ncol, nrow = nrow)
}

# Devuelve el iso3 del pais que contiene cada punto (NA si cae en el mar).
pais_de_coordenada <- function(lon, lat, paises = NULL, indice = NULL) {
  if (is.null(paises)) paises <- leer_paises()
  stopifnot(length(lon) == length(lat))
  if (is.null(indice)) indice <- construir_indice(paises)

  res <- rep(NA_character_, length(lon))
  ok  <- !is.na(lon) & !is.na(lat) & lon >= -180 & lon <= 180 & lat >= -90 & lat <= 90
  if (!any(ok)) return(res)

  celda <- (pmin(pmax(floor((lat[ok] + 90) / indice$res) + 1, 1), indice$nrow) - 1) *
             indice$ncol +
           pmin(pmax(floor((lon[ok] + 180) / indice$res) + 1, 1), indice$ncol)
  pos   <- which(ok)

  for (cl in unique(celda)) {
    sel <- pos[celda == cl]
    for (k in indice$celdas[[cl]]) {
      sel <- sel[is.na(res[sel])]
      if (!length(sel)) break
      pol <- indice$poligonos[[k]]
      x <- lon[sel]; y <- lat[sel]

      dentro <- .en_anillo(x, y, pol$anillos[[1]])
      if (length(pol$anillos) > 1) {
        for (h in 2:length(pol$anillos)) {
          hueco <- .en_anillo(x, y, pol$anillos[[h]])
          dentro[hueco] <- FALSE
        }
      }
      if (any(dentro)) res[sel[dentro]] <- paises[[pol$pais_idx]]$iso3
    }
  }
  res
}

# =============================================================================
# 3. CARGAR EMISIONES POR PAIS
# =============================================================================
# EDGAR entrega CO2 en Gg (= kilotoneladas). 1 Gg = 1 kt = 1000 t.
cargar_emisiones_pais <- function(archivo = RUTA_DATOS("emisiones_co2_pais_edgar.csv")) {
  d <- utils::read.csv(archivo, stringsAsFactors = FALSE)
  d$co2_mt <- d$co2_gg / 1000          # megatoneladas de CO2
  d
}

# OWID: pais-anio con decenas de indicadores (incluye LULUCF, per capita, PIB)
cargar_owid <- function(archivo = RUTA_DATOS("owid-co2-data.csv")) {
  utils::read.csv(archivo, stringsAsFactors = FALSE)
}

# Serie de un pais concreto por ISO3 o por nombre
serie_pais <- function(pais, datos = NULL) {
  if (is.null(datos)) datos <- cargar_emisiones_pais()
  if (nchar(pais) == 3) datos[datos$iso3 == toupper(pais), ]
  else datos[grepl(pais, datos$pais, ignore.case = TRUE), ]
}

# =============================================================================
# 4. LEER LA MALLA (RASTER) DE EMISIONES
# =============================================================================
# Usa GDAL por linea de comandos (viene con QGIS). Si tienes terra, es mas
# comodo: terra::rast("datos/co2_grid_edgar_2023.tif")
# Devuelve data.frame lon, lat, co2 (toneladas por celda de 0.1 grados).
malla_a_tabla <- function(anio = 2023, archivo = NULL, verbose = TRUE) {
  if (is.null(archivo)) archivo <- RUTA_DATOS(sprintf("co2_grid_edgar_%s.tif", anio))
  stopifnot(file.exists(archivo))

  tmp <- tempfile(fileext = ".xyz")
  status <- system2("gdal_translate",
                    c("-q", "-of", "XYZ", shQuote(archivo), shQuote(tmp)),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("gdal_translate fallo. Verifica que GDAL este en el PATH.")

  d <- utils::read.table(tmp, col.names = c("lon", "lat", "co2"),
                         colClasses = c("numeric", "numeric", "numeric"))
  unlink(tmp)
  d <- d[d$co2 > 0, ]            # celdas sin emision
  if (verbose) message(sprintf("%d celdas con emision > 0 (%s)", nrow(d), anio))
  d
}

# =============================================================================
# 5. AGREGAR LA MALLA POR PAIS
# =============================================================================
# Suma las celdas de 0.1 grados por pais (asignacion por centroide de celda).
# Para resultados exactos por area, use exactextractr (ver 03_analisis_espacial.R).
agregar_malla_por_pais <- function(malla, paises = NULL) {
  if (is.null(paises)) paises <- leer_paises()
  malla$iso3 <- pais_de_coordenada(malla$lon, malla$lat, paises)
  agregado <- stats::aggregate(co2 ~ iso3, data = malla, FUN = sum, na.rm = TRUE)
  agregado[!is.na(agregado$iso3), ]
}
