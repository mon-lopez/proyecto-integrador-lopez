# =============================================================================
# 02_ejemplo_coordenadas.R
# De coordenadas a pais: asignacion de puntos y union con las emisiones.
# Ejecutar desde la carpeta del proyecto:  Rscript R/02_ejemplo_coordenadas.R
# =============================================================================
source("R/00_utilidades.R")

# --- 1. Poligonos e indice espacial -----------------------------------------
paises <- leer_paises()                 # Natural Earth 1:10m (258 entidades)
indice <- construir_indice(paises)      # rejilla de 5 grados para acelerar

# --- 2. Un punto cualquiera --------------------------------------------------
pais_de_coordenada(-74.07, 4.71, paises, indice)   # "COL"

# --- 3. Un conjunto de puntos ------------------------------------------------
sitios <- data.frame(
  sitio = c("Bogota", "Medellin", "Cali", "Barranquilla", "Santa Marta",
            "Quito", "Lima", "Sao Paulo"),
  lon = c(-74.07, -75.56, -76.53, -74.78, -74.20,
          -78.51, -77.04, -46.63),
  lat = c( 4.71,   6.24,   3.45,  10.96,  11.24,
           0.22, -12.05, -23.55)
)
sitios$iso3 <- pais_de_coordenada(sitios$lon, sitios$lat, paises, indice)
print(sitios, row.names = FALSE)

# --- 4. Unir con las emisiones del pais --------------------------------------
edgar <- cargar_emisiones_pais()
sitios$co2_pais_mt <- edgar$co2_mt[match(
  paste(sitios$iso3, 2023), paste(edgar$iso3, edgar$anio))]
print(sitios[, c("sitio", "iso3", "co2_pais_mt")], row.names = FALSE)

# --- 5. Limite de precision: Leticia (frontera triple) -----------------------
# Natural Earth 1:10m suaviza el saliente de Leticia (~20 km de ancho) y el
# centro del pueblo cae del lado brasilero.
leticia <- c(-69.9406, -4.2153)
cat("\nNatural Earth 10m :",
    pais_de_coordenada(leticia[1], leticia[2], paises, indice), "\n")

gadm_col <- leer_paises(RUTA_DATOS("gadm41_COL_0.json"),
                        campo_iso = "GID_0", campo_nombre = "COUNTRY")
cat("GADM nivel 0      :",
    pais_de_coordenada(leticia[1], leticia[2], gadm_col,
                       construir_indice(gadm_col)), "\n")
cat("(GADM lo resuelve porque esta a escala 1:1m; ver FUENTES_CO2_ESPACIAL.md)\n")

# --- 6. Puntos generados al azar (prueba de escala) --------------------------
set.seed(1)
N <- 200000
lon <- runif(N, -180, 180); lat <- runif(N, -90, 90)
t0 <- Sys.time()
aleatorio <- pais_de_coordenada(lon, lat, paises, indice)
cat(sprintf("\n%d puntos aleatorios en %.1f s | %% sobre tierra: %.1f%%\n",
            N, as.numeric(difftime(Sys.time(), t0, units = "secs")),
            100 * mean(!is.na(aleatorio))))

# --- 7. Mapa de los sitios ---------------------------------------------------
png("datos/fig_sitios.png", width = 1100, height = 800, res = 120)
op <- par(mar = c(2, 2, 3, 1))
plot(NA, xlim = c(-85, -35), ylim = c(-25, 15), xlab = "", ylab = "",
     bty = "n", main = "Sitios georreferenciados")
for (p in paises) {
  if (is.na(p$iso3)) next
  bb <- p$bbox
  if (bb["xmax"] < -85 || bb["xmin"] > -35 ||
      bb["ymax"] < -25 || bb["ymin"] > 15) next
  for (pol in p$poligonos) {
    xy <- pol[[1]]
    polygon(xy[, 1], xy[, 2], col = "grey92", border = "grey70", lwd = 0.5)
  }
}
points(sitios$lon, sitios$lat, pch = 21, bg = "#c0392b", cex = 1.4)
text(sitios$lon, sitios$lat, sitios$sitio, pos = 3, cex = 0.8, col = "#333333")
par(op)
invisible(dev.off())
cat("Mapa guardado en datos/fig_sitios.png\n")
