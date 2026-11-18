# =============================================================================
# 04_climate_trace.R
# Descarga los activos georreferenciados de Climate TRACE para un pais.
# Cada activo trae lat/lon (centroid) y su emision -> datos de punto reales.
#
# Uso:  Rscript R/04_climate_trace.R [ISO3] [ANIO]
#       Rscript R/04_climate_trace.R COL 2023
# =============================================================================
source("R/00_utilidades.R")
suppressPackageStartupMessages({library(jsonlite); library(curl)})

args <- commandArgs(trailingOnly = TRUE)
pais <- if (length(args) >= 1) toupper(args[1]) else "COL"
anio <- if (length(args) >= 2) as.integer(args[2]) else 2023L
gas  <- "co2"          # "co2", "ch4", "co2e_100yr", ...

# --- 1. Paginar la API -------------------------------------------------------
pagina <- function(pais, anio, gas, limite = 100, offset = 0) {
  url <- sprintf(paste0("https://api.climatetrace.org/v7/sources?",
                        "gadmId=%s&year=%d&gas=%s&limit=%d&offset=%d"),
                 pais, anio, gas, limite, offset)
  h <- new_handle()
  res <- curl_fetch_memory(url, handle = h)
  if (res$status_code != 200) stop("HTTP ", res$status_code, " en ", url)
  jsonlite::fromJSON(rawToChar(res$content))
}

cat("Descargando activos de", pais, anio, "gas =", gas, "\n")
todo <- list(); offset <- 0; limite <- 100
repeat {
  d <- pagina(pais, anio, gas, limite, offset)
  if (length(d) == 0 || is.null(nrow(d))) break
  todo[[length(todo) + 1]] <- d
  if ((offset + nrow(d)) %% 1000 == 0) cat("  ", offset + nrow(d), "registros\n")
  if (nrow(d) < limite) break
  offset <- offset + limite
  if (offset > 200000) break          # salvaguarda
}
if (!length(todo)) stop("La API no devolvio registros.")

activos <- jsonlite::rbind_pages(todo)   # une paginas sin problemas de row.names
rownames(activos) <- NULL
# Centroid viene como data.frame anidado
activos$lon <- activos$centroid$longitude
activos$lat <- activos$centroid$latitude
activos$centroid <- NULL

cat(sprintf("\n%d activos | emision total: %.2f Mt\n",
            nrow(activos), sum(activos$emissionsQuantity, na.rm = TRUE) / 1e6))
print(sort(table(activos$sector), decreasing = TRUE))

# --- 2. Guardar --------------------------------------------------------------
salida <- RUTA_DATOS(sprintf("climate_trace_%s_%d.csv", pais, anio))
write.csv(activos, salida, row.names = FALSE, fileEncoding = "UTF-8")
cat("Guardado en", salida, "\n")

# --- 3. Comprobar que las coordenadas caen en el pais esperado ---------------
paises <- leer_paises(); idx <- construir_indice(paises)
activos$iso3_calc <- pais_de_coordenada(activos$lon, activos$lat, paises, idx)
tab <- table(activos$iso3_calc, useNA = "ifany")
cat("\nPais asignado a los centroides (top 5):\n")
print(head(sort(tab, decreasing = TRUE), 5))
cat(sprintf("Coinciden con %s: %.1f %%\n", pais,
            100 * mean(activos$iso3_calc == pais, na.rm = TRUE)))

# --- 4. Mapa ----------------------------------------------------------------
bb <- c(min(activos$lon, na.rm = TRUE), max(activos$lon, na.rm = TRUE),
        min(activos$lat, na.rm = TRUE), max(activos$lat, na.rm = TRUE))

png(RUTA_DATOS(sprintf("fig_puntos_%s.png", pais)), width = 1000, height = 1000, res = 120)
op <- par(mar = c(3, 3, 3, 1))
plot(NA, xlim = bb[1:2], ylim = bb[3:4], xlab = "lon", ylab = "lat", bty = "l",
     main = sprintf("Climate TRACE: activos de %s (%d)", pais, anio))
for (p in paises) {
  if (is.na(p$iso3)) next
  b <- p$bbox
  if (b["xmax"] < bb[1] || b["xmin"] > bb[2] ||
      b["ymax"] < bb[3] || b["ymin"] > bb[4]) next
  for (pol in p$poligonos) {
    xy <- pol[[1]]
    polygon(xy[, 1], xy[, 2], col = "grey94", border = "grey75", lwd = 0.4)
  }
}
q <- activos$emissionsQuantity
cex <- 0.3 + 2.2 * sqrt(pmax(q, 0) / max(q, na.rm = TRUE))
points(activos$lon, activos$lat, pch = 21, cex = cex,
       bg = adjustcolor("#c0392b", 0.55), col = NA)
par(op)
invisible(dev.off())
cat("Mapa guardado en", RUTA_DATOS(sprintf("fig_puntos_%s.png", pais)), "\n")
