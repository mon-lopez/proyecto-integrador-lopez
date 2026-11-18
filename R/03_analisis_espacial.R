# =============================================================================
# 03_analisis_espacial.R
# Analisis de la malla 0.1 grados: mapa, agregacion por pais y autocorrelacion
# espacial (Moran's I global y local) sin materializar la matriz de pesos.
#
# Requiere:  install.packages(c("terra", "sf"))
# Ejecutar:  Rscript R/03_analisis_espacial.R
# =============================================================================
source("R/00_utilidades.R")

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Falta el paquete 'terra'. Instalalo con: install.packages('terra')")
}
suppressPackageStartupMessages(library(terra))

# --- 1. Leer la malla --------------------------------------------------------
r <- rast(RUTA_DATOS("co2_grid_edgar_2023.tif"))
cat("Malla:", ncell(r), "celdas |", paste(res(r), collapse = " x "), "grados\n")
cat("Total global:", round(global(r, "sum", na.rm = TRUE)[[1]] / 1e9, 2), "Gt\n")

# --- 2. Mapa global (escala log: la distribucion es muy asimetrica) ----------
png(RUTA_DATOS("fig_mapa_global.png"), width = 1400, height = 700, res = 110)
op <- par(mar = c(2, 2, 3, 4))
plot(log10(r + 1), col = hcl.colors(50, "inferno", rev = TRUE),
     main = "Emisiones de CO2 2023 (log10 toneladas por celda de 0.1°)",
     axes = FALSE, box = FALSE, plg = list(cex = 0.7))
par(op)
invisible(dev.off())
cat("Mapa guardado en datos/fig_mapa_global.png\n")

# --- 3. Agregacion por pais --------------------------------------------------
if (requireNamespace("sf", quietly = TRUE)) {
  suppressPackageStartupMessages(library(sf))
  paises <- st_make_valid(st_read(RUTA_DATOS("paises_ne_10m.geojson"), quiet = TRUE))
  val <- terra::extract(r, vect(paises), fun = sum, na.rm = TRUE)
  paises$co2_mt <- val[, 2] / 1e6
  tabla <- st_drop_geometry(paises[, c("ISO_A3_EH", "ADMIN", "co2_mt")])
  tabla <- tabla[order(-tabla$co2_mt), ]
  cat("\nTop 10 por agregacion de la malla (Mt):\n")
  print(head(tabla, 10), row.names = FALSE)

  ofic <- cargar_emisiones_pais()
  ofic <- ofic[ofic$anio == 2023, c("iso3", "co2_mt")]
  cmp <- merge(tabla, ofic, by.x = "ISO_A3_EH", by.y = "iso3",
               suffixes = c("_malla", "_oficial"))
  cmp$dif_pct <- 100 * (cmp$co2_mt_malla - cmp$co2_mt_oficial) / cmp$co2_mt_oficial
  cat("\nDiferencia malla vs tabla oficial:",
      round(median(abs(cmp$dif_pct), na.rm = TRUE), 1), "% (mediana absoluta)\n")
  cat("(la malla deja fuera emisiones asignadas a celdas oceanicas;\n",
      " ver la advertencia en FUENTES_CO2_ESPACIAL.md)\n")
}

# =============================================================================
# 4. AUTOCORRELACION ESPACIAL
# =============================================================================
# Se agrega a 2 grados para que el calculo sea tratable y se usan solo celdas
# con emision > 0 (los ceros son estructurales: oceano, desierto, hielo).
ra <- terra::aggregate(r, fact = 20, fun = sum, na.rm = TRUE)
nr <- nrow(ra); nc <- ncol(ra)
M  <- matrix(terra::values(ra), nrow = nr, ncol = nc, byrow = TRUE)
z  <- log10(M + 1)
z[M == 0 | is.na(M)] <- NA
n  <- sum(!is.na(z))
cat(sprintf("\nRetícula de trabajo: %d x %d | celdas usadas: %d\n", nr, nc, n))

# Vecinos tipo rook (4) sin construir la matriz de pesos: puro desplazamiento
shift <- function(Z, di, dj) {
  Zs <- matrix(NA_real_, nr, nc)
  ir <- seq_len(nr); jc <- seq_len(nc)
  rd <- ir[(ir + di) >= 1 & (ir + di) <= nr]
  cd <- jc[(jc + dj) >= 1 & (jc + dj) <= nc]
  if (!length(rd) || !length(cd)) return(Zs)
  Zs[rd, cd] <- Z[rd + di, cd + dj]
  Zs
}

# Rezago espacial = promedio de vecinos validos (pesos fila-estandarizados)
rezago <- function(Z) {
  S <- matrix(0, nr, nc); C <- matrix(0, nr, nc)
  for (d in list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))) {
    Zs <- shift(Z, d[1], d[2])
    S  <- S + ifelse(is.na(Zs), 0, Zs)
    C  <- C + (!is.na(Zs))
  }
  out <- S / C
  out[is.na(Z)] <- NA
  out
}

Znn <- z
Znn[is.na(Znn)] <- 0
Zc  <- Znn - mean(z, na.rm = TRUE)

moran_I <- function(V) sum(V * rezago(V), na.rm = TRUE) / sum(V^2)
I_obs <- moran_I(Zc)
cat(sprintf("Moran's I global = %.4f   (valor esperado bajo azar: %.5f)\n",
            I_obs, -1 / (n - 1)))

set.seed(1)
per <- replicate(199, moran_I(matrix(sample(Zc), nr, nc)))
pval <- (1 + sum(per >= I_obs)) / (1 + length(per))
cat(sprintf("p-valor por permutaciones (199): %.3f\n", pval))

# --- 5. Moran local (LISA) ---------------------------------------------------
m2    <- sum(Zc^2) / n
I_loc <- Zc * rezago(Zc) / m2

p <- rast(nrows = nr, ncols = nc, xmin = xmin(ra), xmax = xmax(ra),
          ymin = ymin(ra), ymax = ymax(ra), crs = crs(ra))
values(p) <- as.vector(t(I_loc))
names(p) <- "I_local"

png(RUTA_DATOS("fig_lisa.png"), width = 1400, height = 700, res = 110)
op <- par(mar = c(2, 2, 3, 4))
plot(p, col = hcl.colors(50, "RdBu", rev = TRUE),
     main = "Indice de Moran local (2023, escala 2°)",
     axes = FALSE, box = FALSE, plg = list(cex = 0.7))
par(op)
invisible(dev.off())
cat("Mapa LISA guardado en datos/fig_lisa.png\n")
