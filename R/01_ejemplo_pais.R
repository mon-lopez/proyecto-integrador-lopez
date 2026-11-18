# =============================================================================
# 01_ejemplo_pais.R
# Serie de emisiones por pais: descarga, filtro y grafico.
# Ejecutar desde la carpeta del proyecto:  Rscript R/01_ejemplo_pais.R
# =============================================================================
source("R/00_utilidades.R")

# --- 1. Cargar el panel pais-anio de EDGAR -----------------------------------
edgar <- cargar_emisiones_pais()
cat("Filas:", nrow(edgar), "| paises:", length(unique(edgar$iso3)),
    "| anios:", min(edgar$anio), "-", max(edgar$anio), "\n\n")

# --- 2. Filtrar por pais -----------------------------------------------------
pais <- "COL"                      # cambia por cualquier ISO3
col  <- serie_pais(pais, edgar)
cat("Serie de", col$pais[1], "- ultimos 10 anios:\n")
print(tail(col[, c("anio", "co2_mt")], 10), row.names = FALSE)

# --- 3. Comparar varios paises -----------------------------------------------
seleccion <- c("COL", "MEX", "ARG", "CHL", "PER", "BRA")
comp <- edgar[edgar$iso3 %in% seleccion & edgar$anio >= 1990, ]
comp$co2_mt <- comp$co2_gg / 1000

# --- 4. Ranking del ultimo anio ----------------------------------------------
ultimo <- max(edgar$anio)
top <- edgar[edgar$anio == ultimo, ]
top <- top[order(-top$co2_mt), ]
cat("\nTop 10 emisores", ultimo, "(Mt CO2):\n")
print(head(top[, c("pais", "co2_mt")], 10), row.names = FALSE)

# --- 5. Grafico --------------------------------------------------------------
png("datos/fig_serie_paises.png", width = 1200, height = 700, res = 120)
op <- par(mar = c(4, 4.5, 3, 1), las = 1)
plot(NA, xlim = range(comp$anio), ylim = c(0, max(comp$co2_mt, na.rm = TRUE)),
     xlab = "Anio", ylab = "Mt CO2 (fosil)", bty = "l",
     main = "Emisiones de CO2 - America Latina")
cols <- setNames(c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"),
                 seleccion)
for (p in seleccion) {
  d <- comp[comp$iso3 == p, ]
  d <- d[order(d$anio), ]
  lines(d$anio, d$co2_mt, col = cols[p], lwd = 2)
  text(max(d$anio), tail(d$co2_mt, 1), p, col = cols[p], pos = 4, cex = 0.8)
}
par(op)
invisible(dev.off())
cat("\nGrafico guardado en datos/fig_serie_paises.png\n")
