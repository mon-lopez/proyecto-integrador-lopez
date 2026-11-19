# ============================================================================
# 02_sampling.R — Teoría del muestreo: diseño y estimación del total de CO2
# ----------------------------------------------------------------------------
# Marco muestral: 9.657 celdas EDGAR de 0,1° que cubren Colombia (N = 9.657).
# Población objetivo: total de emisiones de CO2 de 2023 (parámetro = 102,14 Mt).
#
# Estratificación: cinco estratos definidos con INFORMACIÓN AUXILIAR DISPONIBLE
# ANTES del operativo (emisiones EDGAR de 2015), no con la variable de interés
# de 2023. Es lo que haría un diseño real: el marco se construye con la
# cartografía de emisiones del año base.
#
# Diseños comparados por Monte Carlo:
#   (1) MAS          : aleatorio simple, igual probabilidad
#   (2) ESTRAT-PROP  : estratificado con asignación proporcional
#   (3) ESTRAT-NEYMAN: estratificado con asignación óptima de Neyman
#   (4) GRTS         : espacialmente balanceado, igual probabilidad
#   (5) GRTS-NEYMAN  : GRTS dentro de cada estrato con asignación Neyman
#                      (diseño final del estudio)
# Estimadores: Horvitz-Thompson y de razón con auxiliar (Climate TRACE 2023).
# Salidas: datos/procesados/resultados_muestreo.rds, muestra_grts.rds,
#          output/figures/muestreo_*, output/maps/muestreo_*, tables/02_*
# ============================================================================
source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages(library(spsurvey))

grilla <- readRDS(file.path(RUTAS$procesados, "grilla_co2.rds"))
grilla$reg <- as.character(grilla$region)
grilla$estrato_aux <- factor(dplyr::ntile(grilla$co2_2015, 5), levels = 1:5,
                             labels = c("E1 (muy baja)", "E2", "E3", "E4", "E5 (muy alta)"))
grilla$estrato_reg <- grilla$reg
# Auxiliar para el estimador de razón: total Climate TRACE 2023 de la celda
cols_ct <- grep("^ct_", names(grilla), value = TRUE)
grilla$ct_total <- rowSums(st_drop_geometry(grilla[, cols_ct]))

N <- nrow(grilla); Y <- grilla$co2_2023; X <- grilla$ct_total
TOTAL_REAL <- sum(Y); TOTAL_X <- sum(X)
paso(sprintf("MARCO MUESTRAL: N = %d celdas | total EDGAR real = %.2f Mt | total CT = %.2f Mt",
             N, TOTAL_REAL / 1e6, TOTAL_X / 1e6))

E <- grilla$estrato_aux
resumen_estratos <- data.frame(
  estrato = levels(E),
  N_h = as.numeric(table(E)),
  media_2015 = as.numeric(tapply(grilla$co2_2015, E, mean)),
  media_2023 = as.numeric(tapply(Y, E, mean)),
  S_h = as.numeric(tapply(Y, E, sd)),
  total_2023_mt = as.numeric(tapply(Y, E, sum)) / 1e6,
  share_pct = 100 * as.numeric(tapply(Y, E, sum)) / TOTAL_REAL
)
print(resumen_estratos, row.names = FALSE, digits = 3)
guardar_tabla(resumen_estratos |> mutate(across(where(is.numeric), ~ round(.x, 3))),
              "02_estratos")

N_h <- as.numeric(table(E)); S_h <- as.numeric(tapply(Y, E, sd)); S_h[is.na(S_h)] <- 0
suaviza <- function(a, n, Nh) {
  a <- pmax(round(a), 2)
  while (sum(a) > n) {
    cand <- which(a > 2 & a == max(a)); a[cand[1]] <- a[cand[1]] - 1
  }
  while (sum(a) < n) a[which.max(Nh - a)] <- a[which.max(Nh - a)] + 1
  as.numeric(a)
}
asignacion_prop <- function(n) suaviza(n * N_h / sum(N_h), n, N_h)
asignacion_neyman <- function(n) suaviza(n * (N_h * S_h) / sum(N_h * S_h), n, N_h)
n_para_cv <- function(cv) ceiling((sum(N_h * S_h) / (cv * TOTAL_REAL))^2)
msg("n para CV = 10% (Neyman): ", n_para_cv(0.10), " | CV = 5%: ", n_para_cv(0.05))

# ---- Selección -------------------------------------------------------------
selecciona_mas <- function(n) sample.int(N, n)
selecciona_estrat <- function(n_h) unlist(lapply(seq_along(N_h), function(h) {
  idx <- which(E == levels(E)[h]); sample(idx, n_h[h])
}))
selecciona_grts <- function(n) {
  s <- spsurvey::grts(sframe = grilla, n_base = n)$sites_base
  match(s$celda_edgar, grilla$celda_edgar)
}
selecciona_grts_estrat <- function(n_h) unlist(lapply(seq_along(N_h), function(h) {
  idx <- which(E == levels(E)[h]); dd <- grilla[idx, ]
  s <- spsurvey::grts(sframe = dd, n_base = n_h[h])$sites_base
  idx[match(s$celda_edgar, dd$celda_edgar)]
}))

# ---- Estimadores -----------------------------------------------------------
estima_mas <- function(idx) {
  n <- length(idx)
  c(est = N * mean(Y[idx]), se = sqrt(N^2 * (1 - n / N) * var(Y[idx]) / n))
}
estima_estrat <- function(idx) {
  est <- 0; vv <- 0
  for (h in seq_along(N_h)) {
    ih <- idx[E[idx] == levels(E)[h]]
    if (!length(ih)) next
    nh <- length(ih); Nh <- N_h[h]
    est <- est + Nh * mean(Y[ih])
    vv <- vv + Nh^2 * (1 - nh / Nh) * var(Y[ih]) / nh
  }
  c(est = est, se = sqrt(vv))
}
estima_razon <- function(idx) {
  w <- 1 / (as.numeric(table(E)[E[idx]]) / N_h[E[idx]] * (length(idx) / N))
  w <- ifelse(is.finite(w) & w > 0, w, 0)
  R <- sum(w * Y[idx]) / sum(w * X[idx])
  est <- TOTAL_X * R
  e <- Y[idx] - R * X[idx]
  se <- TOTAL_X * sqrt(sum((w * e)^2)) / sum(w * X[idx])
  c(est = est, se = se)
}

# ---- Estudio Monte Carlo ---------------------------------------------------
niveles_n <- c(100, 200, 400, 800)
R_aleat <- 400; R_grts <- 50
t0 <- Sys.time()
acumula <- function(etiqueta, n, reps, selector) {
  est <- se <- est_r <- se_r <- numeric(reps)
  for (r in seq_len(reps)) {
    idx <- selector()
    if (n >= N) idx <- seq_len(N)
    eh <- if (etiqueta == "MAS") estima_mas(idx) else estima_estrat(idx)
    er <- estima_razon(idx)
    est[r] <- eh["est"]; se[r] <- eh["se"]; est_r[r] <- er["est"]; se_r[r] <- er["se"]
  }
  data.frame(diseño = etiqueta, n = n, reps = reps,
             sesgo_rel = 100 * (mean(est) - TOTAL_REAL) / TOTAL_REAL,
             cv = 100 * sd(est) / TOTAL_REAL,
             rmse_rel = 100 * sqrt(mean((est - TOTAL_REAL)^2)) / TOTAL_REAL,
             cobertura95 = 100 * mean(abs(est - TOTAL_REAL) <= 1.96 * se),
             cv_razon = 100 * sd(est_r) / TOTAL_REAL,
             cobertura_razon95 = 100 * mean(abs(est_r - TOTAL_REAL) <= 1.96 * se_r))
}

filas <- list()
for (n in niveles_n) {
  filas[[length(filas) + 1]] <- acumula("MAS", n, R_aleat, function() selecciona_mas(n))
  nh_p <- asignacion_prop(n)
  filas[[length(filas) + 1]] <- acumula("Estratificado proporcional", n, R_aleat,
                                        function() selecciona_estrat(nh_p))
  nh_n <- asignacion_neyman(n)
  filas[[length(filas) + 1]] <- acumula("Estratificado Neyman", n, R_aleat,
                                        function() selecciona_estrat(nh_n))
}
filas[[length(filas) + 1]] <- acumula("GRTS", 400, R_grts, function() selecciona_grts(400))
nh_grts <- asignacion_neyman(400)
filas[[length(filas) + 1]] <- acumula("GRTS-Neyman", 400, R_grts,
                                      function() selecciona_grts_estrat(nh_grts))
comparacion <- bind_rows(filas) |> mutate(across(where(is.numeric), ~ round(.x, 2)))
comparacion$eficiencia_rel <- round(100 * (comparacion$cv[comparacion$diseño == "MAS" &
                                                            comparacion$n == 400] /
                                             comparacion$cv)^2, 1)
msg("Monte Carlo terminado en ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
guardar_tabla(comparacion, "02_muestreo_comparacion")
print(as.data.frame(comparacion[comparacion$n == 400, ]), row.names = FALSE)

# ---- Diseño final: GRTS con asignación Neyman (n = 400) --------------------
set.seed(CONST$semilla + 1)
idx_final <- selecciona_grts_estrat(nh_grts)
muestra <- grilla[idx_final, ]
muestra$N_h <- N_h[match(muestra$estrato_aux, levels(E))]
muestra$n_h <- as.numeric(table(muestra$estrato_aux)[muestra$estrato_aux])
muestra$peso <- muestra$N_h / muestra$n_h

ht <- estima_estrat(idx_final); ra <- estima_razon(idx_final)
# IC bootstrap (percentil) para el estimador de razón
set.seed(CONST$semilla + 2)
boot <- replicate(1000, {
  ii <- unlist(lapply(levels(E), function(h) {
    cand <- idx_final[E[idx_final] == h]
    if (length(cand) < 2) return(cand)
    sample(cand, length(cand), replace = TRUE)
  }))
  estima_razon(ii)["est"]
})
ic_razon <- quantile(boot, c(.025, .975), names = FALSE)
msg(sprintf("Muestra final n=%d | HT = %.2f Mt (EE %.2f) | Razón = %.2f Mt IC95 [%.2f; %.2f] | real = %.2f",
            nrow(muestra), ht["est"] / 1e6, ht["se"] / 1e6, ra["est"] / 1e6,
            ic_razon[1] / 1e6, ic_razon[2] / 1e6, TOTAL_REAL / 1e6))

est_region <- muestra |>
  st_drop_geometry() |>
  group_by(region) |>
  summarise(n = n(), celdas_region = sum(N_h / n_h),
            estimado_mt = sum(co2_2023 * peso) / 1e6,
            real_mt = sum(Y[grilla$region == first(region)]) / 1e6,
            .groups = "drop") |>
  mutate(error_rel = round(100 * (estimado_mt - real_mt) / real_mt, 2),
         across(where(is.numeric), ~ round(.x, 2)))
guardar_tabla(est_region, "02_muestreo_por_region")

# ---- Figuras ---------------------------------------------------------------
p_frame <- ggplot(grilla) +
  geom_sf(aes(fill = region), colour = NA, linewidth = 0, alpha = .45) +
  geom_sf(data = muestra, colour = "black", size = .5, shape = 21, fill = "white") +
  scale_fill_manual(values = PAL_REGION, name = "Región") +
  labs(title = "Marco muestral y muestra final (GRTS, n = 400)",
       subtitle = "Puntos: celdas seleccionadas con asignación de Neyman por estrato de emisión",
       caption = "Marco: 9.657 celdas EDGAR de 0,1° (~11 km)") +
  tema() + theme(axis.text = element_blank(), axis.title = element_blank(),
                 panel.grid = element_blank())
guardar_fig(p_frame, "muestreo_01_marco_muestra", 8.5, 9, dir = RUTAS$mapas)

p_cv <- comparacion |>
  tidyr::pivot_longer(c(cv, cv_razon), names_to = "estimador", values_to = "valor") |>
  mutate(estimador = recode(estimador, cv = "Horvitz-Thompson", cv_razon = "Razón (CT 2023)")) |>
  filter(!is.na(n)) |>
  ggplot(aes(n, valor, colour = diseño, linetype = estimador)) +
  geom_line(linewidth = .9) + geom_point(size = 2.2) +
  scale_x_continuous(breaks = niveles_n) +
  labs(title = "Error de muestreo relativo según diseño y tamaño de muestra",
       subtitle = "CV (%) del total de emisiones; 400 réplicas Monte Carlo (50 para GRTS)",
       x = "Tamaño de muestra (n)", y = "CV del total estimado (%)",
       colour = "Diseño", linetype = "Estimador") +
  tema()
guardar_fig(p_cv, "muestreo_02_cv_disenos", 9.5, 5)

p_pesos <- ggplot(muestra, aes(peso, fill = estrato_aux)) +
  geom_histogram(bins = 20, colour = "white") +
  scale_x_log10(labels = label_number(big.mark = " ")) +
  scale_fill_viridis_d(option = "viridis", name = "Estrato") +
  labs(title = "Factores de expansión de la muestra final",
       subtitle = sprintf("Peso = N_h/n_h; rango %.0f-%.0f", min(muestra$peso), max(muestra$peso)),
       x = "Factor de expansión", y = "Celdas") +
  tema()
guardar_fig(p_pesos, "muestreo_03_pesos", 8.5, 4.5)

dmat <- as.matrix(st_distance(st_transform(muestra, CONST$crs_proj))) / 1000
dnn <- apply(dmat, 2, function(x) sort(x)[2])
p_nn <- ggplot(data.frame(dist_nn = dnn), aes(dist_nn)) +
  geom_histogram(bins = 20, fill = "#1B9E77", colour = "white") +
  labs(title = "Distancia al vecino muestral más cercano",
       subtitle = sprintf("Media %.1f km; el balance espacial de GRTS evita conglomerados",
                          mean(dnn)),
       x = "Distancia al vecino más cercano (km)", y = "Celdas") +
  tema()
guardar_fig(p_nn, "muestreo_04_vecino_cercano", 8, 4.5)

# ---- Resultados ------------------------------------------------------------
resumen <- list(
  N = N, total_real_mt = TOTAL_REAL / 1e6, total_ct_mt = TOTAL_X / 1e6,
  media_celda_t = mean(Y), estratos = resumen_estratos,
  asignacion_prop_400 = asignacion_prop(400), asignacion_neyman_400 = asignacion_neyman(400),
  n_para_cv10 = n_para_cv(0.10), n_para_cv5 = n_para_cv(0.05),
  comparacion = comparacion,
  estimacion_final = c(ht_mt = ht["est"] / 1e6, ht_ee_mt = ht["se"] / 1e6,
                       razon_mt = ra["est"] / 1e6, razon_ee_mt = ra["se"] / 1e6,
                       li_mt = ic_razon[1] / 1e6, ls_mt = ic_razon[2] / 1e6,
                       real_mt = TOTAL_REAL / 1e6,
                       error_rel_ht = 100 * (ht["est"] - TOTAL_REAL) / TOTAL_REAL,
                       error_rel_razon = 100 * (ra["est"] - TOTAL_REAL) / TOTAL_REAL),
  por_region = est_region,
  balance_espacial = c(dist_nn_media_km = mean(dnn), dist_nn_sd_km = sd(dnn),
                       peso_medio = mean(muestra$peso), peso_min = min(muestra$peso),
                       peso_max = max(muestra$peso))
)
saveRDS(resumen, file.path(RUTAS$procesados, "resultados_muestreo.rds"))
saveRDS(muestra, file.path(RUTAS$procesados, "muestra_grts.rds"))

paso("MUESTREO COMPLETADO")
cat(sprintf("Real %.2f Mt | HT %.2f (EE %.2f) | Razón %.2f IC95 [%.2f; %.2f]\n",
            TOTAL_REAL / 1e6, ht["est"] / 1e6, ht["se"] / 1e6, ra["est"] / 1e6,
            ic_razon[1] / 1e6, ic_razon[2] / 1e6))
