# ============================================================================
# 02b_diseno_final.R — Diseño muestral final con estrato de certeza
# ----------------------------------------------------------------------------
# El estudio Monte Carlo de 02_sampling.R muestra que, sin estrato de certeza,
# el CV del total es ~20 % con n = 400 porque el 1 % de las celdas concentra
# ~64 % de las emisiones. La solución estándar en poblaciones tan asimétricas
# es un estrato de certeza (inclusión forzosa de las unidades dominantes).
#
# Diseño final:
#   - Estrato de certeza: celdas con emisión 2015 >= percentil 98 (censo).
#   - Resto: cuatro estratos por cuartiles de emisión 2015, asignación Neyman,
#     selección GRTS dentro de cada estrato.
#   - n total = 400.
# Estimador: censo del estrato de certeza + Horvitz-Thompson del resto.
# ============================================================================
source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages(library(spsurvey))

grilla <- readRDS(file.path(RUTAS$procesados, "grilla_co2.rds"))
ruta_res <- file.path(RUTAS$procesados, "resultados_muestreo.rds")
res <- if (file.exists(ruta_res)) readRDS(ruta_res) else list()

N <- nrow(grilla); Y <- grilla$co2_2023; TOTAL_REAL <- sum(Y)
n_total <- 400
set.seed(CONST$semilla + 11)

# ---- Estrato de certeza ----------------------------------------------------
umbral_certeza <- quantile(grilla$co2_2015, 0.98)
certeza <- grilla$co2_2015 >= umbral_certeza
n_certeza <- sum(certeza)
msg("Estrato de certeza: ", n_certeza, " celdas (", round(100 * n_certeza / N, 2),
    "%) con emisión 2015 >= ", fmt_num(umbral_certeza, 0), " t; concentran ",
    round(100 * sum(Y[certeza]) / TOTAL_REAL, 1), "% del CO2 de 2023")

# ---- Estratos del resto ----------------------------------------------------
resto <- grilla[!certeza, ]
resto$estrato <- factor(dplyr::ntile(resto$co2_2015, 4), levels = 1:4,
                        labels = c("R1", "R2", "R3", "R4"))
N_h <- as.numeric(table(resto$estrato))
S_h <- as.numeric(tapply(resto$co2_2023, resto$estrato, sd))
n_resto <- n_total - n_certeza
a <- n_resto * (N_h * S_h) / sum(N_h * S_h)
a <- pmax(round(a), 2)
while (sum(a) > n_resto) a[which.max(a)] <- a[which.max(a)] - 1
while (sum(a) < n_resto) a[which.max(N_h - a)] <- a[which.max(N_h - a)] + 1
names(a) <- levels(resto$estrato)
msg("Asignación Neyman del resto (n = ", n_resto, "): ",
    paste(names(a), a, sep = "=", collapse = " | "))

# ---- Selección GRTS dentro de cada estrato del resto -----------------------
idx_resto <- unlist(lapply(levels(resto$estrato), function(h) {
  dd <- resto[resto$estrato == h, ]
  s <- spsurvey::grts(sframe = dd, n_base = a[[h]])$sites_base
  match(s$celda_edgar, grilla$celda_edgar)
}))
idx_final <- c(which(certeza), idx_resto)
muestra <- grilla[idx_final, ]
muestra$estrato_final <- ifelse(muestra$celda_edgar %in% grilla$celda_edgar[certeza],
                                "Certeza", as.character(resto$estrato[
                                  match(muestra$celda_edgar, resto$celda_edgar)]))
muestra$N_h <- ifelse(muestra$estrato_final == "Certeza", n_certeza,
                      N_h[match(muestra$estrato_final, levels(resto$estrato))])
muestra$n_h <- ifelse(muestra$estrato_final == "Certeza", n_certeza,
                      a[match(muestra$estrato_final, levels(resto$estrato))])
muestra$peso <- muestra$N_h / muestra$n_h

# ---- Estimación -----------------------------------------------------------
est <- sum(Y[certeza])
vv <- 0
for (h in levels(resto$estrato)) {
  ih <- idx_resto[grilla$celda_edgar[idx_resto] %in% resto$celda_edgar[resto$estrato == h]]
  Nh <- N_h[match(h, levels(resto$estrato))]; nh <- length(ih)
  est <- est + Nh * mean(Y[ih])
  vv <- vv + Nh^2 * (1 - nh / Nh) * var(Y[ih]) / nh
}
se <- sqrt(vv)
ic <- est + c(-1.96, 1.96) * se
msg(sprintf("Diseño con certeza: n = %d (%d certeza + %d muestreadas) | HT = %.2f Mt (EE %.2f) | IC95 [%.2f; %.2f] | real = %.2f | error = %.1f%%",
            nrow(muestra), n_certeza, length(idx_resto), est / 1e6, se / 1e6,
            ic[1] / 1e6, ic[2] / 1e6, TOTAL_REAL / 1e6,
            100 * (est - TOTAL_REAL) / TOTAL_REAL))

# ---- Estimador de razón con la emisión de 2015 como auxiliar --------------
# El inventario de 2015 está disponible para TODAS las celdas del marco y
# explica la mayor parte de la variación de 2023; usarlo como auxiliar reduce
# drásticamente la varianza del total del resto.
X2015 <- grilla$co2_2015
X_resto_total <- sum(X2015[!certeza])
w_resto <- 1 / (muestra$n_h[muestra$estrato_final != "Certeza"] /
                  muestra$N_h[muestra$estrato_final != "Certeza"])
y_r <- Y[idx_resto]; x_r <- X2015[idx_resto]
R_razon <- sum(w_resto * y_r) / sum(w_resto * x_r)
est_razon <- sum(Y[certeza]) + X_resto_total * R_razon
e_r <- y_r - R_razon * x_r
se_razon <- X_resto_total * sqrt(sum((w_resto * e_r)^2)) / sum(w_resto * x_r)
ic_razon <- est_razon + c(-1.96, 1.96) * se_razon
msg(sprintf("Estimador de razón (auxiliar 2015): %.2f Mt (EE %.2f) | IC95 [%.2f; %.2f] | error = %.1f%% | R = %.4f",
            est_razon / 1e6, se_razon / 1e6, ic_razon[1] / 1e6, ic_razon[2] / 1e6,
            100 * (est_razon - TOTAL_REAL) / TOTAL_REAL, R_razon))

# CV esperado del diseño (analítico, sin Monte Carlo)
cv_esperado <- 100 * se / TOTAL_REAL
msg("CV esperado del diseño con certeza (HT): ", round(cv_esperado, 2), "% | razón: ",
    round(100 * se_razon / TOTAL_REAL, 2), "%")

# ---- Figuras --------------------------------------------------------------
p <- ggplot(grilla) +
  geom_sf(aes(fill = region), colour = NA, linewidth = 0, alpha = .45) +
  geom_sf(data = muestra |> filter(estrato_final != "Certeza"),
          colour = "black", size = .5, shape = 21, fill = "white") +
  geom_sf(data = muestra |> filter(estrato_final == "Certeza"),
          colour = "#B2182B", size = .8, shape = 21, fill = "#B2182B", alpha = .8) +
  scale_fill_manual(values = PAL_REGION, name = "Región") +
  labs(title = "Diseño muestral final con estrato de certeza (n = 400)",
       subtitle = sprintf("%d celdas de inclusión forzosa (rojo) + %d celdas GRTS/Neyman (blanco)",
                          n_certeza, length(idx_resto)),
       caption = "Marco: 9.657 celdas EDGAR de 0,1°") +
  tema() + theme(axis.text = element_blank(), axis.title = element_blank(),
                 panel.grid = element_blank())
guardar_fig(p, "muestreo_05_diseno_certeza", 8.5, 9, dir = RUTAS$mapas)

# ---- Guardado -------------------------------------------------------------
diseno_certeza <- list(
  umbral_certeza = umbral_certeza, n_certeza = n_certeza,
  n_resto = n_resto, asignacion_resto = a,
  estimacion = c(est_mt = est / 1e6, se_mt = se / 1e6, li_mt = ic[1] / 1e6,
                 ls_mt = ic[2] / 1e6, real_mt = TOTAL_REAL / 1e6,
                 error_rel_pct = 100 * (est - TOTAL_REAL) / TOTAL_REAL),
  estimacion_razon = c(est_mt = est_razon / 1e6, se_mt = se_razon / 1e6,
                       li_mt = ic_razon[1] / 1e6, ls_mt = ic_razon[2] / 1e6,
                       error_rel_pct = 100 * (est_razon - TOTAL_REAL) / TOTAL_REAL,
                       R = R_razon),
  cv_esperado_pct = cv_esperado,
  muestra = data.frame(n = nrow(muestra), peso_medio = mean(muestra$peso),
                       n_certeza = n_certeza)
)
res$diseno_certeza <- diseno_certeza
saveRDS(res, file.path(RUTAS$procesados, "resultados_muestreo.rds"))
saveRDS(muestra, file.path(RUTAS$procesados, "muestra_final.rds"))

tabla <- data.frame(
  metrica = c("n total", "Celdas de certeza", "n para el resto",
              "HT: estimación (Mt)", "HT: error estándar (Mt)", "HT: IC95 inferior",
              "HT: IC95 superior", "Razón: estimación (Mt)", "Razón: error estándar (Mt)",
              "Razón: IC95 inferior", "Razón: IC95 superior", "Valor real EDGAR (Mt)",
              "HT: error relativo (%)", "Razón: error relativo (%)",
              "CV esperado HT (%)", "CV esperado razón (%)"),
  valor = round(c(nrow(muestra), n_certeza, length(idx_resto), est / 1e6, se / 1e6,
                  ic[1] / 1e6, ic[2] / 1e6, est_razon / 1e6, se_razon / 1e6,
                  ic_razon[1] / 1e6, ic_razon[2] / 1e6, TOTAL_REAL / 1e6,
                  100 * (est - TOTAL_REAL) / TOTAL_REAL,
                  100 * (est_razon - TOTAL_REAL) / TOTAL_REAL,
                  cv_esperado, 100 * se_razon / TOTAL_REAL), 3)
)
guardar_tabla(tabla, "02b_muestreo_diseno_final")
print(tabla, row.names = FALSE)
paso("DISEÑO FINAL COMPLETADO")
