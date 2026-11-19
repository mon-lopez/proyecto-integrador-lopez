# ============================================================================
# 05_ml.R — Machine Learning: predictores espaciales y pronóstico a 2030
# ----------------------------------------------------------------------------
# Parte A. Modelo espacial de celdas: ¿qué explica la densidad de emisiones?
#          RF / XGBoost / OLS con validación cruzada espacial por bloques.
# Parte B. Pronóstico nacional 2025-2030: modelo de tasas de crecimiento
#          (RF/XGB/OLS) con validación de origen móvil + escenarios macro.
# Parte C. Ajuste por política FNCER (modelo contable) y probabilidad de que
#          Colombia reduzca sus emisiones en 2030 frente a 2024, 2015 y la
#          meta de la Ley 2169 de 2021.
# Salidas: datos/procesados/resultados_ml.rds, output/figures/ml_*, tables/05_*
# ============================================================================
source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({
  library(ranger); library(xgboost); library(spdep)
})
HILOS <- max(1, parallel::detectCores() - 1)

grilla <- readRDS(file.path(RUTAS$procesados, "grilla_co2.rds"))
panel  <- readRDS(file.path(RUTAS$procesados, "panel_nacional.rds"))
set.seed(CONST$semilla)

metrica <- function(obs, pred) {
  e <- obs - pred
  c(RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
    R2 = 1 - sum(e^2) / sum((obs - mean(obs))^2), sesgo = mean(pred - obs))
}

# ============================================================================
# PARTE A — MODELO ESPACIAL DE CELDAS
# ============================================================================
paso("A. MODELO ESPACIAL DE DENSIDAD DE EMISIONES (n = 9.657 celdas)")
a <- grilla |>
  st_drop_geometry() |>
  mutate(
    celda_edgar = celda_edgar,
    log_dens = log_dens_2023,
    log_pob = log1p(poblacion_2020),
    log_popdens = log1p(pop_dens_2020),
    log_dist = log1p(dist_fncer_km),
    region = factor(region),
    urbano = factor(urbano)
  ) |>
  select(celda_edgar, x_9377, y_9377, log_dens, log_pob, log_popdens, urbano,
         region, log_dist, n_proy, cap_mw, area_km2, lat, lon) |>
  na.omit()
msg("Datos completos: ", nrow(a), " celdas")

form_a <- ~ log_pob + log_popdens + urbano + region + log_dist +
  n_proy + cap_mw + area_km2 + lat + lon
X <- model.matrix(form_a, data = a)[, -1]
y <- a$log_dens

# Bloques espaciales para validación cruzada (k-means sobre coordenadas)
km <- kmeans(a[, c("x_9377", "y_9377")], centers = 10, nstart = 20)
a$bloque <- km$cluster

ajusta_rf <- function(idx, ntree = 400, imp = "none") {
  ranger(x = X[idx, , drop = FALSE], y = y[idx], num.trees = ntree,
         mtry = max(2, floor(ncol(X) / 3)), importance = imp,
         num.threads = HILOS, seed = CONST$semilla)
}
pred_rf <- function(m, newX) predict(m, data = newX, num.threads = HILOS)$predictions
ajusta_xgb <- function(idx) {
  xgb.train(params = list(objective = "reg:squarederror", max_depth = 4, eta = 0.05,
                          subsample = 0.8, colsample_bytree = 0.8, nthread = HILOS),
            data = xgb.DMatrix(X[idx, , drop = FALSE], label = y[idx]),
            nrounds = 300, verbose = 0)
}
pred_xgb <- function(m, newX) predict(m, xgb.DMatrix(newX))

metricas_a <- list(); pred_cv <- list()
for (b in sort(unique(a$bloque))) {
  tr <- which(a$bloque != b); te <- which(a$bloque == b)
  rf <- ajusta_rf(tr); xg <- ajusta_xgb(tr)
  ols <- lm(y ~ ., data = data.frame(y = y[tr], X[tr, , drop = FALSE]))
  pr <- tibble(obs = y[te], rf = pred_rf(rf, X[te, , drop = FALSE]),
               xgb = pred_xgb(xg, X[te, , drop = FALSE]),
               ols = as.numeric(predict(ols, data.frame(X[te, , drop = FALSE]))),
               bloque = b)
  pred_cv[[as.character(b)]] <- pr
  metricas_a[[as.character(b)]] <- rbind(
    RF = metrica(pr$obs, pr$rf), XGB = metrica(pr$obs, pr$xgb), OLS = metrica(pr$obs, pr$ols))
}
pred_cv <- bind_rows(pred_cv)
metricas_a <- bind_rows(lapply(names(metricas_a), function(k) {
  m <- metricas_a[[k]]
  data.frame(bloque = as.integer(k), modelo = rownames(m), m, row.names = NULL)
}))
resumen_a <- metricas_a |>
  group_by(modelo) |>
  summarise(RMSE = sqrt(mean(RMSE^2)), MAE = mean(MAE), R2 = mean(R2), .groups = "drop") |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))
print(as.data.frame(resumen_a), row.names = FALSE)

# Comparación CV aleatoria vs espacial (optimismo por autocorrelación)
folds_aleat <- sample(rep(1:10, length.out = nrow(a)))
rmse_aleat <- numeric(10)
for (f in 1:10) {
  tr <- which(folds_aleat != f); te <- which(folds_aleat == f)
  rf <- ajusta_rf(tr, ntree = 300)
  rmse_aleat[f] <- metrica(y[te], pred_rf(rf, X[te, , drop = FALSE]))["RMSE"]
}
comp_cv <- tibble(
  esquema = c("Validación aleatoria (10-fold)", "Validación espacial (10 bloques)"),
  RMSE = c(mean(rmse_aleat), resumen_a$RMSE[resumen_a$modelo == "RF"])
)
print(as.data.frame(comp_cv), row.names = FALSE)

# Modelo final e importancia
rf_final <- ajusta_rf(seq_len(nrow(a)), ntree = 600, imp = "permutation")
imp <- ranger::importance(rf_final) |>
  tibble::enframe(name = "variable", value = "importancia_perm") |>
  arrange(desc(importancia_perm)) |>
  mutate(importancia_rel = 100 * importancia_perm / max(importancia_perm),
         across(where(is.numeric), ~ round(.x, 2)))
print(head(imp, 10), row.names = FALSE)

# Dependencia parcial de las tres variables numéricas más importantes
top3 <- head(imp$variable[imp$variable %in% c("log_pob", "log_popdens", "log_dist",
                                              "n_proy", "cap_mw", "area_km2", "lat", "lon")], 3)
pd_df <- bind_rows(lapply(top3, function(v) {
  g <- seq(min(a[[v]]), max(a[[v]]), length.out = 25)
  d <- a[rep(1, length(g)), , drop = FALSE]
  d[[v]] <- g
  names(d)[names(d) == "log_dens"] <- "log_dens"
  Xg <- model.matrix(form_a, data = d)[, -1]
  # alinear columnas con la matriz de entrenamiento
  Xg <- Xg[, colnames(X), drop = FALSE]
  tibble(variable = v, valor = g, pred = as.numeric(pred_rf(rf_final, Xg)))
}))

# Residuales y autocorrelación residual
a$resid <- y - pred_rf(rf_final, X)
nb <- dnearneigh(as.matrix(a[, c("x_9377", "y_9377")]), 0, 12000, use_kd_tree = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
moran_resid <- moran.test(a$resid, lw, zero.policy = TRUE)
msg("Moran's I de los residuales del RF: ", round(moran_resid$estimate[1], 3),
    " (p = ", fmt_p(moran_resid$p.value), ")")

# Triangulación con Climate TRACE
tri <- grilla |> st_drop_geometry()
cols_comb <- c("ct_power", "ct_transporte", "ct_edificaciones", "ct_industria", "ct_ffo")
tri$ct_total <- rowSums(tri[, cols_comb], na.rm = TRUE)
tri <- tri[tri$ct_total > 0 & tri$co2_2023 > 0, ]
triangulacion <- c(
  n_celdas = nrow(tri), cor_log = cor(log1p(tri$co2_2023), log1p(tri$ct_total)),
  ratio_total = sum(tri$co2_2023) / sum(tri$ct_total),
  cor_spearman = cor(tri$co2_2023, tri$ct_total, method = "spearman"))
msg("Triangulación EDGAR vs Climate TRACE: r(log) = ", round(triangulacion["cor_log"], 3),
    " | razón total = ", round(triangulacion["ratio_total"], 2))

# ============================================================================
# PARTE B — PRONÓSTICO NACIONAL 2025-2030
# ============================================================================
paso("B. PRONÓSTICO NACIONAL DE CO2 A 2030")
# Se modela la tasa de crecimiento anual del CO2. La muestra comienza en 2001
# porque la participación hidroeléctrica (principal determinante de corto plazo
# de las emisiones del sistema eléctrico) solo está disponible desde 2000.
b <- panel |>
  filter(anio >= 2001, anio <= 2024) |>
  mutate(
    g_co2 = 100 * (log(co2_owid) - log(lag(co2_owid))),
    g_pib = 100 * (log(pib_const2015) - log(lag(pib_const2015))),
    g_pob = 100 * (log(poblacion) - log(lag(poblacion))),
    pct_hidro = tidyr::replace_na(pct_hidro, 70)
  ) |>
  filter(!is.na(g_co2), !is.na(g_pib), !is.na(g_pob)) |>
  select(anio, g_co2, g_pib, g_pob, pct_hidro) |>
  na.omit()
msg("Serie de crecimiento: ", nrow(b), " años (", min(b$anio), "-", max(b$anio), ")")

db <- data.frame(g_co2 = b$g_co2, g_pib = b$g_pib, g_pob = b$g_pob,
                 pct_hidro = b$pct_hidro)
Xb <- as.matrix(db[, -1]); yb <- db$g_co2

ajusta_rf_b <- function(idx) {
  ranger(x = Xb[idx, , drop = FALSE], y = yb[idx], num.trees = 400, mtry = 2, num.threads = HILOS, seed = CONST$semilla)
}
pred_xgb_b <- function(m, nuevo) predict(m, xgb.DMatrix(nuevo))
ajusta_xgb_b <- function(idx) {
  xgb.train(params = list(objective = "reg:squarederror", max_depth = 3, eta = 0.08,
                          subsample = 0.8, colsample_bytree = 0.8, nthread = HILOS),
            data = xgb.DMatrix(Xb[idx, , drop = FALSE], label = yb[idx]),
            nrounds = 200, verbose = 0)
}

origenes <- 2014:2023
cv_b <- bind_rows(lapply(origenes, function(o) {
  tr <- which(b$anio <= o); te <- which(b$anio == o + 1)
  if (!length(te) || length(tr) < 10) return(NULL)
  rf <- ajusta_rf_b(tr); xg <- ajusta_xgb_b(tr)
  ols <- lm(g_co2 ~ ., data = db[tr, , drop = FALSE])
  pr <- tibble(anio = o + 1, obs = yb[te],
               RF = predict(rf, data = Xb[te, , drop = FALSE], num.threads = HILOS)$predictions,
               XGB = pred_xgb_b(xg, Xb[te, , drop = FALSE]),
               OLS = as.numeric(predict(ols, db[te, , drop = FALSE])),
               NAIVE = mean(yb[tr]))
  pr$ENS <- rowMeans(pr[, c("RF", "XGB", "OLS", "NAIVE")])
  pr
}))
cv_res <- bind_rows(lapply(c("RF", "XGB", "OLS", "NAIVE", "ENS"), function(m) {
  data.frame(modelo = m, t(metrica(cv_b$obs, cv_b[[m]])), row.names = NULL)
}))
cv_res[, -1] <- lapply(cv_res[, -1], function(x) round(as.numeric(x), 3))
print(cv_res, row.names = FALSE)
mejor <- cv_res$modelo[which.min(cv_res$RMSE)]
sd_cv <- cv_res$RMSE[cv_res$modelo == mejor]
msg("Mejor modelo: ", mejor, " | RMSE de validación = ", sd_cv, " p.p.")

rf_b <- ajusta_rf_b(seq_along(yb))
xgb_b <- ajusta_xgb_b(seq_along(yb))
ols_b <- lm(g_co2 ~ ., data = db)
modelo_pron <- function(nuevo) {
  switch(mejor,
         "RF" = as.numeric(predict(rf_b, data = nuevo, num.threads = HILOS)$predictions),
         "XGB" = as.numeric(pred_xgb_b(xgb_b, nuevo)),
         "OLS" = as.numeric(predict(ols_b, data.frame(nuevo))),
         "NAIVE" = rep(mean(yb), nrow(nuevo)),
         "ENS" = rowMeans(cbind(
           as.numeric(predict(rf_b, data = nuevo, num.threads = HILOS)$predictions),
           as.numeric(pred_xgb_b(xgb_b, nuevo)),
           as.numeric(predict(ols_b, data.frame(nuevo))),
           rep(mean(yb), nrow(nuevo)))))
}

# --- Escenarios -------------------------------------------------------------
escenarios <- tibble::tribble(
  ~escenario, ~g_pib, ~g_pob, ~fncer_2030,
  "BAU",        3.0,   0.9,   8,
  "Moderado",   3.0,   0.9,  20,
  "Acelerado",  3.0,   0.9,  35
) |> mutate(escenario = factor(escenario, levels = c("BAU", "Moderado", "Acelerado")))

fncer_2024 <- panel$pct_fncer_nc[panel$anio == 2024]
fncer_bau_path <- seq(fncer_2024, 8, length.out = 7)[-1]
hidro_ref <- mean(panel$pct_hidro[panel$anio >= 2015 & panel$anio <= 2024], na.rm = TRUE)
co2_2024 <- panel$co2_owid[panel$anio == 2024]
msg("CO2 2024 = ", round(co2_2024, 2), " Mt | FNCER 2024 = ", round(fncer_2024, 2),
    "% | hidro ref = ", round(hidro_ref, 1), "%")

simula_crecimiento <- function(esc, B = 2000) {
  e <- escenarios[escenarios$escenario == esc, ]
  sim <- matrix(NA_real_, B, 6)
  prev <- rep(log(co2_2024), B)
  for (i in 1:6) {
    Xn <- cbind(g_pib = e$g_pib, g_pob = e$g_pob, pct_hidro = hidro_ref)
    mu <- modelo_pron(Xn)
    g <- mu + rnorm(B, 0, sd_cv)
    prev <- prev + g / 100
    sim[, i] <- exp(prev)
  }
  colnames(sim) <- paste0("y", 2025:2030)
  sim
}

# --- Ajuste contable por FNCER (adicional respecto al BAU) ------------------
gen_2024 <- panel$gen_total_gwh[panel$anio == 2024]
crec_gen <- mean(diff(log(panel$gen_total_gwh[panel$anio >= 2015 & panel$anio <= 2024])))
ef_marginal <- 0.55
ef_rango <- c(0.45, 0.65)

ajuste_fncer <- function(esc, B = 2000) {
  e <- escenarios[escenarios$escenario == esc, ]
  path <- seq(fncer_2024, e$fncer_2030, length.out = 7)[-1]
  gen <- gen_2024 * exp(crec_gen * (1:6))
  delta <- (path - fncer_bau_path) / 100 * gen * ef_marginal / 1000   # Mt CO2
  mat <- matrix(rep(delta, each = B), B, 6) *
    matrix(runif(B * 6, ef_rango[1], ef_rango[2]) / ef_marginal, B, 6)
  mat
}

pron <- list()
for (esc in levels(escenarios$escenario)) {
  base <- simula_crecimiento(esc)
  pron[[esc]] <- base - ajuste_fncer(esc)
}
pron_long <- bind_rows(lapply(names(pron), function(esc) {
  m <- pron[[esc]]
  tibble(escenario = esc, anio = 2025:2030, media = colMeans(m),
         mediana = apply(m, 2, median),
         li = apply(m, 2, quantile, .025), ls = apply(m, 2, quantile, .975))
})) |> mutate(escenario = factor(escenario, levels = levels(escenarios$escenario)))
print(as.data.frame(pron_long |> filter(anio %in% c(2025, 2030)) |>
                      mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)

# --- Probabilidades de reducción -------------------------------------------
umbral <- c(`2024` = co2_2024, `2015` = panel$co2_owid[panel$anio == 2015])
prob <- bind_rows(lapply(names(pron), function(esc) {
  v <- pron[[esc]][, 6]
  tibble(escenario = esc, co2_2030_media = mean(v), co2_2030_mediana = median(v),
         p_menor_2024 = mean(v < umbral["2024"]),
         p_menor_2015 = mean(v < umbral["2015"]),
         reduccion_vs_2024_pct = 100 * (mean(v) / co2_2024 - 1))
})) |> mutate(escenario = factor(escenario, levels = levels(escenarios$escenario)),
              across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(prob), row.names = FALSE)

co2_bau_2030 <- median(pron[["BAU"]][, 6])
adicional <- tibble(
  escenario = c("Moderado", "Acelerado"),
  co2_2030 = c(median(pron[["Moderado"]][, 6]), median(pron[["Acelerado"]][, 6]))
) |>
  mutate(co2_bau = co2_bau_2030,
         reduccion_vs_bau_pct = 100 * (1 - co2_2030 / co2_bau),
         across(where(is.numeric), ~ round(.x, 2)))
print(as.data.frame(adicional), row.names = FALSE)

# --- Contraste con la Ley 2169 de 2021 -------------------------------------
no_co2_2024 <- panel$ghg_total_owid[panel$anio == 2024] - co2_2024
ghg_2030_bau <- co2_bau_2030 + no_co2_2024
ley <- tibble(
  concepto = c("CO2 2030 BAU", "CO2 2030 Moderado", "CO2 2030 Acelerado",
               "No-CO2 + LULUCF 2024", "GEI 2030 implícito (BAU)",
               "Techo legal 2030 (Ley 2169)", "Brecha vs techo legal (BAU)",
               "CO2 compatible con el techo si el no-CO2 no cambia"),
  mt_co2e = c(co2_bau_2030, median(pron[["Moderado"]][, 6]), median(pron[["Acelerado"]][, 6]),
              no_co2_2024, ghg_2030_bau, CONST$meta_2030_mtco2e,
              ghg_2030_bau - CONST$meta_2030_mtco2e,
              CONST$meta_2030_mtco2e - no_co2_2024)
) |> mutate(across(where(is.numeric), ~ round(.x, 2)))
print(as.data.frame(ley), row.names = FALSE)

# --- Sensibilidad: deriva temporal en la tasa de crecimiento ---------------
# Un modelo con tendencia explícita captura la reducción secular de la
# intensidad de carbono y produce una trayectoria más baja. Se reporta como
# cota inferior del escenario BAU.
db_t <- cbind(db, t = b$anio - 2000)
ols_t <- lm(g_co2 ~ ., data = db_t)
pred_t <- predict(ols_t, data.frame(g_pib = 3, g_pob = 0.9, pct_hidro = hidro_ref, t = 25:30))
co2_bau_t <- co2_2024 * exp(sum(pred_t) / 100)
msg("Sensibilidad con deriva temporal: CO2 2030 BAU = ", round(co2_bau_t, 1), " Mt")

# --- Contraste estructural OLS log-log -------------------------------------
est <- panel |>
  filter(anio >= 2001, anio <= 2024, !is.na(co2_owid), !is.na(pib_const2015)) |>
  mutate(log_co2 = log(co2_owid), log_pib = log(pib_const2015), log_pob = log(poblacion),
         t = anio - 1990)
m_est <- lm(log_co2 ~ log_pib + log_pob + t, data = est)
proy_est <- exp(predict(m_est, newdata = data.frame(
  log_pib = log(est$pib_const2015[est$anio == 2024] * (1.03)^(1:6)),
  log_pob = log(est$poblacion[est$anio == 2024] * (1.009)^(1:6)),
  t = 2024 + 1:6 - 1990)))
estructural <- tibble(anio = 2025:2030, co2_estructural = proy_est, r2 = summary(m_est)$r.squared)
print(as.data.frame(estructural |> mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)

# ============================================================================
# PARTE C — FIGURAS Y SALIDAS
# ============================================================================
paso("C. FIGURAS Y RESULTADOS")

p_imp <- imp |> filter(!grepl("^region", variable)) |>
  slice_max(importancia_rel, n = 10) |>
  ggplot(aes(reorder(variable, importancia_rel), importancia_rel)) +
  geom_col(fill = "#2E6F9E") + coord_flip() +
  labs(title = "Importancia de las variables en la densidad de emisiones",
       subtitle = "Random Forest (ranger): importancia por permutación, relativa al máximo",
       x = NULL, y = "Importancia relativa (%)") +
  tema()
guardar_fig(p_imp, "ml_01_importancia", 8.5, 5)

p_pd <- ggplot(pd_df, aes(valor, pred)) +
  geom_line(linewidth = 1, colour = "#2E6F9E") +
  facet_wrap(~variable, scales = "free_x") +
  labs(title = "Dependencia parcial de la densidad de emisiones",
       subtitle = "Efecto marginal promedio en log(t CO2/km² + 1)",
       x = "Valor de la variable", y = "Predicción (log densidad)") +
  tema()
guardar_fig(p_pd, "ml_02_dependencia_parcial", 10, 4)

p_obs <- ggplot(pred_cv, aes(obs, rf)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(alpha = .25, size = .7, colour = "#2E6F9E") +
  labs(title = "Predicción vs observación (validación espacial)",
       subtitle = sprintf("Random Forest: RMSE = %.2f, R² = %.2f en escala log",
                          resumen_a$RMSE[resumen_a$modelo == "RF"],
                          resumen_a$R2[resumen_a$modelo == "RF"]),
       x = "log densidad observada", y = "log densidad predicha") +
  tema()
guardar_fig(p_obs, "ml_03_pred_vs_obs", 7, 6)

grilla_res <- grilla |> left_join(a[, c("celda_edgar", "resid")], by = "celda_edgar")
p_res <- ggplot(grilla_res) +
  geom_sf(aes(fill = resid), colour = NA, linewidth = 0) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                       midpoint = 0, name = "Residual") +
  labs(title = "Residuales del modelo espacial de emisiones",
       subtitle = sprintf("Moran's I = %.3f (p = %s): queda estructura espacial sin explicar",
                          moran_resid$estimate[1], fmt_p(moran_resid$p.value))) +
  tema() + theme(axis.text = element_blank(), axis.title = element_blank(),
                 panel.grid = element_blank())
guardar_fig(p_res, "ml_04_residuales", 8.5, 9)

p_cv <- ggplot(comp_cv, aes(esquema, RMSE, fill = esquema)) +
  geom_col(width = .55, show.legend = FALSE) +
  geom_text(aes(label = round(RMSE, 2)), vjust = -.4) +
  scale_fill_manual(values = c("#B2182B", "#1B9E77")) +
  labs(title = "Optimismo de la validación aleatoria frente a la espacial",
       subtitle = "RMSE del mismo Random Forest (escala log)",
       x = NULL, y = "RMSE") +
  tema()
guardar_fig(p_cv, "ml_05_cv_espacial", 7.5, 4.5)

serie_graf <- panel |> filter(anio >= 2005, anio <= 2024) |> select(anio, co2 = co2_owid)
p_fc <- ggplot() +
  geom_line(data = serie_graf, aes(anio, co2), linewidth = 1, colour = "black") +
  geom_ribbon(data = pron_long, aes(anio, ymin = li, ymax = ls, fill = escenario), alpha = .18) +
  geom_line(data = pron_long, aes(anio, media, colour = escenario), linewidth = 1) +
  geom_hline(yintercept = umbral["2015"], linetype = 3, colour = "grey40") +
  annotate("text", x = 2005.5, y = umbral["2015"] - 3, size = 3.1, colour = "grey30",
           label = sprintf("Nivel de 2015: %.1f Mt", umbral["2015"])) +

  scale_fill_manual(values = PAL_ESC, name = "Escenario") +
  scale_colour_manual(values = PAL_ESC, name = "Escenario") +
  labs(title = "Pronóstico de emisiones de CO2 de Colombia a 2030",
       subtitle = "Bandas: intervalo de predicción del 95% (bootstrap de residuales)",
       x = NULL, y = "Mt CO2 / año",
       caption = "Modelo de tasas de crecimiento + ajuste contable por FNCER") +
  tema()
guardar_fig(p_fc, "ml_06_pronostico", 9.5, 5.5)

p_prob <- prob |>
  select(escenario, `vs 2024` = p_menor_2024, `vs 2015` = p_menor_2015) |>
  pivot_longer(-escenario) |>
  ggplot(aes(escenario, value, fill = name)) +
  geom_col(position = "dodge", width = .65) +
  geom_text(aes(label = percent(value, accuracy = 1)),
            position = position_dodge(.65), vjust = -.3, size = 3) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.12)) +
  scale_fill_manual(values = c("vs 2024" = "#EF8A62", "vs 2015" = "#2166AC"), name = NULL) +
  labs(title = "Probabilidad de que el CO2 de 2030 sea inferior a la referencia",
       subtitle = "2.000 simulaciones por escenario", x = NULL, y = "Probabilidad") +
  tema()
guardar_fig(p_prob, "ml_07_probabilidades", 8, 4.5)

guardar_tabla(resumen_a, "05_ml_metricas_celdas")
guardar_tabla(imp, "05_ml_importancia")
guardar_tabla(cv_res, "05_ml_cv_pronostico")
guardar_tabla(pron_long |> mutate(across(where(is.numeric), ~ round(.x, 3))), "05_ml_pronostico")
guardar_tabla(prob, "05_ml_probabilidades")
guardar_tabla(ley, "05_ml_ley2169")
guardar_tabla(comp_cv, "05_ml_cv_espacial")
guardar_tabla(estructural, "05_ml_estructural")

resultados_ml <- list(
  celdas = list(metricas = resumen_a, metricas_bloque = metricas_a, importancia = imp,
                cv_espacial = comp_cv,
                moran_residual = c(I = moran_resid$estimate[1], p = moran_resid$p.value),
                triangulacion = triangulacion, dependencia = pd_df),
  pronostico = list(cv = cv_res, mejor_modelo = mejor, sd_cv = sd_cv, co2_bau_t = co2_bau_t,
                    escenarios = escenarios, pronostico = pron_long, probabilidades = prob,
                    adicional = adicional, ley = ley, estructural = estructural,
                    co2_2024 = co2_2024, co2_2015 = umbral["2015"],
                    no_co2_2024 = no_co2_2024, ef_marginal = ef_marginal,
                    gen_2024 = gen_2024, crec_gen = crec_gen,
                    fncer_2024 = fncer_2024, hidro_ref = hidro_ref)
)
saveRDS(resultados_ml, file.path(RUTAS$procesados, "resultados_ml.rds"))

paso("ML COMPLETADO")
cat(sprintf("CO2 2030: BAU = %.1f | Moderado = %.1f | Acelerado = %.1f Mt\n",
            co2_bau_2030, median(pron[["Moderado"]][, 6]), median(pron[["Acelerado"]][, 6])))
cat(sprintf("Reducción vs BAU: Moderado %.1f%% | Acelerado %.1f%%\n",
            adicional$reduccion_vs_bau_pct[1], adicional$reduccion_vs_bau_pct[2]))
cat(sprintf("P(CO2 2030 < 2024): BAU %.2f | Moderado %.2f | Acelerado %.2f\n",
            prob$p_menor_2024[1], prob$p_menor_2024[2], prob$p_menor_2024[3]))
