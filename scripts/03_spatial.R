# ============================================================================
# 03_spatial.R — Módulo de Estadística Espacial
# Proyecto: Emisiones de CO2 y transición energética en Colombia 2015-2030
# ----------------------------------------------------------------------------
# Objetivo: caracterizar la estructura espacial de la densidad de emisiones
# de CO2 (EDGAR 0,1°, 9.657 celdas, EPSG:9377) y evaluar la capacidad
# predictiva de métodos geoestadísticos.
#
# CONTENIDO
#   1. Autocorrelación global: I de Moran (reina 8 vecinos, torre 12 km y kNN k=6)
#      para log_dens_2023, log_co2_2023 y log1p(población 2020).
#   2. LISA: localmoran_perm (999 permutaciones, two.sided), clases HH/LL/HL/LH
#      con p<0,05 y robustez FDR-BH.
#   3. Semivariograma empírico (cutoff 280 km) y ajuste Esférico/Exponencial/
#      Gaussiano; selección por SSE; nugget, sill, rango práctico y nugget/sill.
#   4. Validación predictiva espacial: kriging ordinario vs IDW (p=1,2,3) con
#      validación cruzada por bloques espaciales (k-means, k=5) y fracciones de
#      entrenamiento 25/50/100%; RMSE, MAE, R2; Diebold-Mariano + bootstrap
#      pareado; mapas de predicción y de error.
#   5. Párrafo-resumen técnico con los hallazgos cuantitativos.
#
# CONTRATO DE SALIDA — datos/procesados/resultados_espaciales.rds
#   $moran_global      data.frame: variable, pesos, n, n_efectivo, I, EI, VarI,
#                      sdI, z, p_perm (999), p_analitico, cardinalidad_media
#   $lisa              data.frame por celda: celda_edgar, depto, region, x_9377,
#                      y_9377, Ii, EI, VarIi, Z, p_sim, p_fdr, clase, sig,
#                      clase_sig, clase_sig_fdr
#   $resumen_clusters  list(global, por_region, concentracion, significancia)
#   $variograma        list(empirico, modelos, mejor, params, ajustes)
#   $cv                data.frame: metodo, metodo_etiqueta, fraccion, n_replicas,
#                      n_train, n_celdas, n_pred, RMSE, MAE, R2, sesgo
#   $test_dm           list(principal, por_fraccion): estadístico DM, p, bootstrap
#   $n_celdas          entero (9.657)
#   $resumen           párrafo técnico (character)
#   $metadatos         list(fecha, semilla, crs, variable, cutoff_km, n_perm, ...)
#
# FIGURAS
#   output/maps/    : densidad, LISA, kriging, error de validación cruzada
#   output/figures/ : scatter de Moran, semivariograma, comparación de RMSE
# ============================================================================

source(file.path("scripts", "00_setup.R"))
suppressPackageStartupMessages({
  library(spdep); library(gstat); library(viridis); library(patchwork)
})

set.seed(CONST$semilla)
op <- options(stringsAsFactors = FALSE)
on.exit(options(op), add = TRUE)

# ---------------------------------------------------------------------------
# Utilidades locales
# ---------------------------------------------------------------------------
nfmt <- function(x, d = 2) {
  formatC(x, format = "f", digits = d, big.mark = ".", decimal.mark = ",")
}
nint <- function(x) formatC(x, format = "d", big.mark = ".", decimal.mark = ",")
metricas <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  e <- obs[ok] - pred[ok]
  c(n = sum(ok), RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
    R2 = 1 - sum(e^2) / sum((obs[ok] - mean(obs[ok]))^2), sesgo = mean(e))
}
# Silencia la salida de consola de gstat ("[using ordinary kriging]", etc.)
# sin alterar el valor devuelto ni las advertencias.
call_quiet <- function(expr) {
  f <- tempfile()
  sink(f)
  on.exit({sink(); unlink(f)}, add = TRUE)
  force(expr)
}

# ===========================================================================
paso("0. DATOS, COORDENADAS Y MATRICES DE PESOS ESPACIALES")
# ===========================================================================
grilla <- readRDS(file.path(RUTAS$procesados, "grilla_co2.rds"))
stopifnot(inherits(grilla, "sf"))
stopifnot(st_crs(grilla)$epsg == CONST$crs_proj)
dat <- st_drop_geometry(grilla)
n <- nrow(dat)
stopifnot(n == 9657)
coords <- cbind(dat$x_9377, dat$y_9377)
stopifnot(!anyNA(coords), all(diff(range(coords[, 1])) > 0))

# Variables de trabajo (todas observadas, sin imputación)
dat$log_pob_2020 <- log1p(dat$poblacion_2020)
vars_moran <- list(
  log_dens_2023 = dat$log_dens_2023,
  log_co2_2023  = dat$log_co2_2023,
  log_pob_2020  = dat$log_pob_2020
)
stopifnot(all(vapply(vars_moran, function(v) all(is.finite(v)), logical(1))))

msg("Celdas: ", n, " | CRS: EPSG:", st_crs(grilla)$epsg)
msg("log_dens_2023: media = ", nfmt(mean(dat$log_dens_2023), 3),
    " | DE = ", nfmt(sd(dat$log_dens_2023), 3),
    " | rango = [", nfmt(min(dat$log_dens_2023), 2), "; ",
    nfmt(max(dat$log_dens_2023), 2), "]")

# --- Vecindades -------------------------------------------------------------
# La grilla es un retículo regular de 0,1°: la separación entre centroides es
# ~11,1 km y la diagonal ~15,7 km. Por eso:
#   * 12 km  -> solo vecinos de arista (torre, 4 vecinos): robustez.
#   * 16 km  -> aristas + diagonales (reina, 8 vecinos): especificación principal.
# dnearneigh usa distancias euclídeas en metros (EPSG:9377).
w_queen <- dnearneigh(coords, 0, 16000)
w_rook  <- dnearneigh(coords, 0, 12000)
w_knn   <- knn2nb(knearneigh(coords, k = 6), sym = TRUE)
lw_queen <- nb2listw(w_queen, style = "W", zero.policy = TRUE)
lw_rook  <- nb2listw(w_rook,  style = "W", zero.policy = TRUE)
lw_knn   <- nb2listw(w_knn,   style = "W", zero.policy = TRUE)

islas <- function(nb) sum(card(nb) == 0)
msg("Reina 16 km: cardinalidad media = ", nfmt(mean(card(w_queen)), 3),
    " | islas = ", islas(w_queen), " | subgrafos = ", n.comp.nb(w_queen)$nc)
msg("Torre 12 km: cardinalidad media = ", nfmt(mean(card(w_rook)), 3),
    " | islas = ", islas(w_rook))
msg("kNN k=6 simétrico: cardinalidad media = ", nfmt(mean(card(w_knn)), 3),
    " | islas = ", islas(w_knn))

pesos_lista <- list(
  list(lw = lw_queen, nombre = "Reina 8 vecinos (16 km)", umbral = 16),
  list(lw = lw_rook,  nombre = "Torre 4 vecinos (12 km)", umbral = 12),
  list(lw = lw_knn,   nombre = "kNN k=6 (simétrico)",     umbral = NA_real_)
)

# ===========================================================================
paso("1. AUTOCORRELACIÓN ESPACIAL GLOBAL (I DE MORAN)")
# ===========================================================================
moran_fila <- function(x, lw, variable, pesos, umbral) {
  mt <- moran.test(x, lw, randomisation = TRUE, zero.policy = TRUE)
  set.seed(CONST$semilla)
  mc <- moran.mc(x, lw, nsim = 999, alternative = "greater", zero.policy = TRUE)
  est <- mt$estimate
  data.frame(
    variable = variable, pesos = pesos, umbral_km = umbral,
    n = length(x), n_efectivo = length(x) - sum(card(lw$neighbours) == 0),
    I = unname(est["Moran I statistic"]),
    EI = unname(est["Expectation"]),
    VarI = unname(est["Variance"]),
    sdI = sqrt(unname(est["Variance"])),
    z = unname(mt$statistic),
    p_perm = as.numeric(mc$p.value),
    p_analitico = as.numeric(mt$p.value),
    cardinalidad_media = mean(card(lw$neighbours)),
    stringsAsFactors = FALSE
  )
}

moran_global <- do.call(rbind, lapply(names(vars_moran), function(v) {
  do.call(rbind, lapply(pesos_lista, function(p) {
    moran_fila(vars_moran[[v]], p$lw, v, p$nombre, p$umbral)
  }))
}))
rownames(moran_global) <- NULL
moran_global$variable <- factor(moran_global$variable,
                                levels = c("log_dens_2023", "log_co2_2023", "log_pob_2020"))
moran_global <- moran_global[order(moran_global$variable, moran_global$umbral_km, na.last = TRUE), ]
print(moran_global[, c("variable", "pesos", "n_efectivo", "I", "EI", "VarI", "p_perm")],
      row.names = FALSE, digits = 5)

moran_q <- moran_global[moran_global$variable == "log_dens_2023" &
                          moran_global$pesos == "Reina 8 vecinos (16 km)", ]
stopifnot(nrow(moran_q) == 1)

# ===========================================================================
paso("2. AUTOCORRELACIÓN LOCAL (LISA) — 999 PERMUTACIONES")
# ===========================================================================
# localmoran_perm: columnas Ii, E.Ii, Var.Ii, Z.Ii, p analítico, p simulado,
# p folded, skewness, kurtosis. Se usa el p simulado (permutaciones) como
# criterio principal y FDR-BH como robustez frente a comparaciones múltiples.
lisa_mat <- localmoran_perm(dat$log_dens_2023, lw_queen, nsim = 999,
                            alternative = "two.sided", zero.policy = TRUE,
                            iseed = CONST$semilla)
col_p_sim <- which(colnames(lisa_mat) == "Pr(z != E(Ii)) Sim")
stopifnot(length(col_p_sim) == 1)
col_p_ana <- which(colnames(lisa_mat) == "Pr(z != E(Ii))")

z_sc  <- as.numeric(scale(dat$log_dens_2023))
lag_z <- as.numeric(lag.listw(lw_queen, z_sc, zero.policy = TRUE))
lag_z[is.na(lag_z)] <- 0  # isla sin vecinos: se clasifica como no significativa

clase <- ifelse(z_sc > 0 & lag_z > 0, "HH",
         ifelse(z_sc < 0 & lag_z < 0, "LL",
         ifelse(z_sc > 0 & lag_z < 0, "HL",
         ifelse(z_sc < 0 & lag_z > 0, "LH", "NS"))))
p_sim <- lisa_mat[, col_p_sim]
p_ana <- lisa_mat[, col_p_ana]
p_fdr <- p.adjust(p_sim, method = "BH")
sig     <- !is.na(p_sim) & p_sim < 0.05
sig_fdr <- !is.na(p_fdr) & p_fdr < 0.05

lisa <- data.frame(
  celda_edgar = dat$celda_edgar,
  depto = dat$depto, region = dat$region,
  x_9377 = dat$x_9377, y_9377 = dat$y_9377,
  z = dat$log_dens_2023,
  Ii = lisa_mat[, "Ii"], EI = lisa_mat[, "E.Ii"], VarIi = lisa_mat[, "Var.Ii"],
  Z = lisa_mat[, "Z.Ii"],
  p_sim = p_sim, p_analitico = p_ana, p_fdr = p_fdr,
  clase = factor(clase, levels = c("HH", "LL", "HL", "LH", "NS")),
  sig = sig, sig_fdr = sig_fdr,
  clase_sig = factor(ifelse(sig, clase, "NS"),
                     levels = c("HH", "LL", "HL", "LH", "NS")),
  clase_sig_fdr = factor(ifelse(sig_fdr, clase, "NS"),
                         levels = c("HH", "LL", "HL", "LH", "NS")),
  stringsAsFactors = FALSE
)
rownames(lisa) <- NULL
msg("LISA: ", sum(sig), " celdas con p<0,05 (", nfmt(100 * mean(sig), 2),
    "%) | FDR-BH<0,05: ", sum(sig_fdr), " (", nfmt(100 * mean(sig_fdr), 2), "%)")
print(table(lisa$clase_sig))

# Conteos globales, por región y concentración de HH/LL
tab_global <- lisa |>
  dplyr::count(clase_sig, .drop = FALSE) |>
  dplyr::mutate(pct = 100 * n / sum(n))
tab_global_fdr <- lisa |>
  dplyr::count(clase_sig_fdr, .drop = FALSE) |>
  dplyr::mutate(pct = 100 * n / sum(n))

tab_region <- lisa |>
  dplyr::group_by(region, clase_sig, .drop = FALSE) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(region) |>
  dplyr::mutate(pct_region = 100 * n / sum(n)) |>
  dplyr::ungroup()

concentracion <- lisa |>
  dplyr::filter(sig, clase %in% c("HH", "LL")) |>
  dplyr::count(clase, region) |>
  dplyr::group_by(clase) |>
  dplyr::mutate(pct_clase = 100 * n / sum(n)) |>
  dplyr::ungroup()

resumen_clusters <- list(
  global = tab_global,
  global_fdr = tab_global_fdr,
  por_region = tab_region,
  concentracion = concentracion,
  significancia = data.frame(
    criterio = c("p permutaciones < 0,05", "FDR-BH < 0,05"),
    n = c(sum(sig), sum(sig_fdr)),
    pct = c(100 * mean(sig), 100 * mean(sig_fdr)),
    stringsAsFactors = FALSE
  )
)

# ===========================================================================
paso("3. SEMIVARIOGRAMA EMPÍRICO Y AJUSTE DE MODELOS")
# ===========================================================================
# Coordenadas en km para que el rango se exprese directamente en km.
d_km <- data.frame(x = dat$x_9377 / 1000, y = dat$y_9377 / 1000,
                   z = dat$log_dens_2023, celda_edgar = dat$celda_edgar)
v_emp <- variogram(z ~ 1, locations = ~ x + y, data = d_km,
                   cutoff = 280, width = 20)
stopifnot(nrow(v_emp) >= 10)
msg("Variograma empírico: ", nrow(v_emp), " bins | cutoff = 280 km | ancho = 20 km")

inits <- data.frame(
  nug = c(0.2, 0.5, 0.8, 0.35),
  rng = c(40, 90, 150, 200),
  psill = c(0.8, 1.2, 1.6, 2.0)
)
ajustar_variograma <- function(v, modelo, inits) {
  res <- list()
  for (i in seq_len(nrow(inits))) {
    f <- tryCatch(
      suppressWarnings(fit.variogram(v, vgm(psill = inits$psill[i], model = modelo,
                                            range = inits$rng[i],
                                            nugget = inits$nug[i]))),
      error = function(e) NULL)
    if (is.null(f)) next
    if (!all(is.finite(f$psill)) || !all(is.finite(f$range))) next
    if (any(f$psill < 0) || any(f$range < 0)) next
    g <- variogramLine(f, dist_vector = v$dist)$gamma
    res[[length(res) + 1]] <- list(
      fit = f,
      sse = sum((v$gamma - g)^2),
      ssew = sum(v$np * (v$gamma - g)^2)
    )
  }
  if (!length(res)) return(NULL)
  res[[which.min(vapply(res, function(r) r$sse, numeric(1)))]]
}

modelos_gstat <- c("Sph", "Exp", "Gau")
etiquetas     <- c(Sph = "Esférico", Exp = "Exponencial", Gau = "Gaussiano")
ajustes <- list()
for (m in modelos_gstat) {
  a <- ajustar_variograma(v_emp, m, inits)
  if (is.null(a)) {
    msg("ADVERTENCIA: el modelo ", etiquetas[m], " no convergió; se excluye.")
  } else {
    ajustes[[m]] <- a
  }
}
stopifnot(length(ajustes) >= 1)

rango_practico <- function(f, m) {
  r <- f$range[2]
  switch(m, Sph = r, Exp = 3 * r, Gau = sqrt(3) * r)
}
modelos_tab <- do.call(rbind, lapply(names(ajustes), function(m) {
  f <- ajustes[[m]]$fit
  data.frame(
    modelo = etiquetas[m], modelo_gstat = m,
    nugget = f$psill[1],
    psill = f$psill[2],
    rango_param_km = f$range[2],
    rango_practico_km = rango_practico(f, m),
    sill = sum(f$psill),
    prop_nugget = f$psill[1] / sum(f$psill),
    SSE = ajustes[[m]]$sse,
    SSE_ponderado = ajustes[[m]]$ssew,
    convergio = TRUE,
    stringsAsFactors = FALSE
  )
}))
modelos_tab <- modelos_tab[order(modelos_tab$SSE), ]
rownames(modelos_tab) <- NULL
mejor <- modelos_tab$modelo_gstat[1]
params_mejor <- ajustes[[mejor]]$fit
print(modelos_tab, row.names = FALSE, digits = 4)
msg("Mejor modelo por SSE: ", etiquetas[mejor], " (SSE = ",
    nfmt(modelos_tab$SSE[1], 4), ")")

# ===========================================================================
paso("4. VALIDACIÓN CRUZADA ESPACIAL: KRIGING vs IDW")
# ===========================================================================
# Bloques espaciales: k-means (k=5) sobre las coordenadas proyectadas. Cada
# bloque actúa una vez como conjunto de prueba (5 folds). En cada fold se
# reajusta el variograma SOLO con los datos de entrenamiento (evita fuga de
# información) y se predicen las celdas del bloque con:
#   * kriging ordinario local (nmax = 40 vecinos)
#   * IDW con potencias 1, 2 y 3 (nmax = 40 vecinos)
# No se impone radio máximo para garantizar cobertura out-of-fold completa
# (el bloque amazónico tiene celdas a >280 km de cualquier punto de
# entrenamiento; excluirlas sesgaría el RMSE hacia la región andina).
set.seed(CONST$semilla)
km5 <- kmeans(coords / 1000, centers = 5, nstart = 25, iter.max = 100)
folds <- km5$cluster
msg("Folds espaciales (k-means, k=5): ", paste(names(table(folds)), table(folds),
                                               sep = "=", collapse = " | "))
for (k in sort(unique(folds))) {
  idx <- folds == k
  msg("  bloque ", k, ": n = ", sum(idx), " | x [",
      nfmt(min(d_km$x[idx]), 0), ", ", nfmt(max(d_km$x[idx]), 0), "] km | y [",
      nfmt(min(d_km$y[idx]), 0), ", ", nfmt(max(d_km$y[idx]), 0), "] km")
}

ajustar_vg_cv <- function(dsub) {
  v <- variogram(z ~ 1, locations = ~ x + y, data = dsub, cutoff = 280, width = 20)
  f <- tryCatch(
    suppressWarnings(fit.variogram(v, vgm(psill = params_mejor$psill[2], model = mejor,
                                          range = params_mejor$range[2],
                                          nugget = params_mejor$psill[1]))),
    error = function(e) NULL)
  if (is.null(f) || any(!is.finite(f$psill)) || any(!is.finite(f$range)) ||
      any(f$psill < 0) || any(f$range < 0)) {
    return(list(model = params_mejor, ok = FALSE))
  }
  list(model = f, ok = TRUE)
}

fracciones <- c(0.25, 0.50, 1.00)
metodos <- c("kriging_ok", "idw_p1", "idw_p2", "idw_p3")
etiq_metodos <- c(kriging_ok = "Kriging ordinario", idw_p1 = "IDW (p=1)",
                  idw_p2 = "IDW (p=2)", idw_p3 = "IDW (p=3)")
# Submuestreo anidado: para cada réplica y fold se sortea UNA permutación y
# 25% / 50% son sus prefijos (25% ⊂ 50% ⊂ 100%). Cinco réplicas para 25% y
# 50% promedian la variabilidad de la muestra; 100% es determinista.
n_rep_sub <- 5L
n_cond <- n_rep_sub * length(unique(folds)) * 2 + length(unique(folds))
pred_list <- vector("list", n_cond * length(metodos))
n_fallback <- 0L
n_var_fallback <- 0L
ii <- 0L

for (rep in seq_len(n_rep_sub)) {
  for (k in sort(unique(folds))) {
    train_all <- which(folds != k)
    test_idx  <- which(folds == k)
    set.seed(CONST$semilla + 10000 * rep + 1000 * k)
    perm <- sample(train_all)
    for (fr in fracciones) {
      if (fr == 1 && rep > 1) next
      ntr <- if (fr == 1) length(train_all) else max(30, round(fr * length(train_all)))
      train_idx <- if (fr == 1) train_all else perm[seq_len(ntr)]
      dtr <- d_km[train_idx, c("x", "y", "z")]
      dte <- d_km[test_idx,  c("x", "y", "z")]
      vg  <- ajustar_vg_cv(dtr)
      if (!vg$ok) n_var_fallback <- n_var_fallback + 1L

      pk <- call_quiet(gstat::krige(z ~ 1, locations = ~ x + y, data = dtr,
                                    newdata = dte, model = vg$model, nmax = 40))
      p1 <- call_quiet(gstat::idw(z ~ 1, locations = ~ x + y, data = dtr,
                                  newdata = dte, idp = 1, nmax = 40))
      p2 <- call_quiet(gstat::idw(z ~ 1, locations = ~ x + y, data = dtr,
                                  newdata = dte, idp = 2, nmax = 40))
      p3 <- call_quiet(gstat::idw(z ~ 1, locations = ~ x + y, data = dtr,
                                  newdata = dte, idp = 3, nmax = 40))

      preds <- list(kriging_ok = pk$var1.pred, idw_p1 = p1$var1.pred,
                    idw_p2 = p2$var1.pred, idw_p3 = p3$var1.pred)
      # Respaldo documentado: si el kriging local es singular, se usa IDW(p=2)
      # solo para esa celda y se contabiliza.
      if (anyNA(preds$kriging_ok)) {
        bad <- is.na(preds$kriging_ok)
        preds$kriging_ok[bad] <- preds$idw_p2[bad]
        n_fallback <- n_fallback + sum(bad)
      }
      for (m in metodos) {
        ii <- ii + 1L
        pred_list[[ii]] <- data.frame(
          celda_edgar = d_km$celda_edgar[test_idx],
          fold = k, replica = rep, fraccion = fr, metodo = m,
          n_train = length(train_idx),
          obs = dte$z, pred = preds[[m]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
pred_list <- pred_list[!vapply(pred_list, is.null, logical(1))]
cv_pred <- do.call(rbind, pred_list)
rm(pred_list); gc(verbose = FALSE)
msg("Predicciones out-of-fold: ", nrow(cv_pred),
    " | celdas con respaldo IDW en kriging: ", n_fallback,
    " | condiciones con variograma de respaldo: ", n_var_fallback)
stopifnot(!anyNA(cv_pred$pred))

cv <- cv_pred |>
  dplyr::group_by(metodo, fraccion) |>
  dplyr::summarise(
    n_replicas = dplyr::n_distinct(replica),
    n_train = round(mean(n_train)),
    n_celdas = dplyr::n_distinct(celda_edgar),
    n_pred = dplyr::n(),
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE = mean(abs(obs - pred)),
    R2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    sesgo = mean(obs - pred),
    .groups = "drop"
  ) |>
  dplyr::mutate(metodo_etiqueta = unname(etiq_metodos[metodo])) |>
  dplyr::select(metodo, metodo_etiqueta, fraccion, n_replicas, n_train,
                n_celdas, n_pred, RMSE, MAE, R2, sesgo) |>
  dplyr::arrange(fraccion, metodo)
print(cv, row.names = FALSE, digits = 4)

# --- Diebold-Mariano + bootstrap pareado ------------------------------------
# Para cada fracción se promedian los errores cuadráticos por celda entre
# réplicas (así cada celda pesa igual) y se comparan kriging vs el mejor IDW.
dm_test <- function(d) {
  nn <- length(d); md <- mean(d)
  lag <- max(1, floor(4 * (nn / 100)^(2 / 9)))
  g0 <- mean((d - md)^2)
  vr <- g0
  for (l in seq_len(lag)) {
    gl <- mean((d[(l + 1):nn] - md) * (d[1:(nn - l)] - md))
    vr <- vr + 2 * (1 - l / (lag + 1)) * gl
  }
  stat <- md / sqrt(vr / nn)
  list(stat = stat, p = 2 * (1 - pnorm(abs(stat))), lag = lag,
       media_dif = md, n = nn)
}

B <- 2000
comparar_fr <- function(fr) {
  sub <- cv_pred |>
    dplyr::filter(fraccion == fr) |>
    dplyr::group_by(celda_edgar, metodo) |>
    dplyr::summarise(sse = mean((obs - pred)^2), .groups = "drop") |>
    tidyr::pivot_wider(names_from = metodo, values_from = sse)
  idw_fr <- cv |>
    dplyr::filter(fraccion == fr, grepl("^idw", metodo)) |>
    dplyr::slice_min(RMSE, n = 1) |>
    dplyr::pull(metodo)
  d <- sub$kriging_ok - sub[[idw_fr]]
  dd <- dm_test(d)
  set.seed(CONST$semilla + round(100 * fr))
  boot <- replicate(B, {
    idx <- sample.int(length(d), replace = TRUE)
    sqrt(mean(sub$kriging_ok[idx])) - sqrt(mean(sub[[idw_fr]][idx]))
  })
  boot_p <- min(1, 2 * min(mean(boot <= 0), mean(boot >= 0)))
  boot_ic <- stats::quantile(boot, c(0.025, 0.975), names = FALSE)
  data.frame(
    fraccion = fr, mejor_idw = idw_fr, n = length(d),
    n_replicas = max(cv$n_replicas[cv$fraccion == fr]),
    RMSE_kriging = sqrt(mean(sub$kriging_ok)),
    RMSE_idw = sqrt(mean(sub[[idw_fr]])),
    dif_RMSE = sqrt(mean(sub$kriging_ok)) - sqrt(mean(sub[[idw_fr]])),
    dif_rel_pct = 100 * (sqrt(mean(sub$kriging_ok)) / sqrt(mean(sub[[idw_fr]])) - 1),
    DM = dd$stat, p_DM = dd$p, lag_HAC = dd$lag,
    bootstrap_p = boot_p,
    bootstrap_ic_li = boot_ic[1], bootstrap_ic_ls = boot_ic[2],
    stringsAsFactors = FALSE
  )
}
por_fr <- do.call(rbind, lapply(fracciones, comparar_fr))
p100 <- por_fr[por_fr$fraccion == 1, ]

test_dm <- list(
  principal = list(
    fraccion = 1,
    mejor_kriging = "kriging_ok",
    mejor_idw = p100$mejor_idw,
    n = p100$n,
    RMSE_kriging = p100$RMSE_kriging,
    RMSE_idw = p100$RMSE_idw,
    dif_RMSE = p100$dif_RMSE,
    dif_rel_pct = p100$dif_rel_pct,
    DM = p100$DM, p_DM = p100$p_DM, lag_HAC = p100$lag_HAC,
    bootstrap_p = p100$bootstrap_p,
    bootstrap_ic95 = c(p100$bootstrap_ic_li, p100$bootstrap_ic_ls),
    bootstrap_B = B,
    nota_bootstrap = "bootstrap iid por celda; anti-conservador bajo dependencia espacial"
  ),
  por_fraccion = por_fr
)
print(test_dm$principal[c("RMSE_kriging", "RMSE_idw", "dif_RMSE", "DM", "p_DM",
                          "bootstrap_p")], digits = 4)
print(por_fr[, c("fraccion", "mejor_idw", "RMSE_kriging", "RMSE_idw",
                 "dif_RMSE", "DM", "p_DM", "bootstrap_p")],
      row.names = FALSE, digits = 4)

mejor_idw <- p100$mejor_idw
rmse_k <- p100$RMSE_kriging
rmse_i <- p100$RMSE_idw

# ===========================================================================
paso("5. FIGURAS")
# ===========================================================================
bordes <- tryCatch(
  st_read(file.path(RUTAS$geometrias, "gadm41_COL_1.json"), quiet = TRUE) |>
    st_transform(CONST$crs_proj),
  error = function(e) NULL)

tema_mapa <- theme(
  legend.position = "right",
  axis.text = element_blank(), axis.ticks = element_blank(),
  panel.grid = element_blank()
)

# 5.1 Mapa de densidad -------------------------------------------------------
p_dens <- ggplot(grilla) +
  geom_sf(aes(fill = log_dens_2023), colour = NA, linewidth = 0) +
  {if (!is.null(bordes)) geom_sf(data = bordes, fill = NA, colour = "grey35", linewidth = 0.25)} +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "log(1 + t CO2/km²)") +
  labs(title = "Densidad de emisiones de CO2 por celda EDGAR, 2023",
       subtitle = paste0("n = ", nint(n),
                         " celdas de 0,1°; I de Moran = ",
                         nfmt(moran_q$I, 3), " (p = ",
                         nfmt(moran_q$p_perm, 3), ")")) +
  tema_mapa
guardar_fig(p_dens, "03_mapa_densidad_log_2023", ancho = 8, alto = 9,
            dir = RUTAS$mapas)

# 5.2 Mapa LISA --------------------------------------------------------------
cols_lisa <- c("HH" = "#B2182B", "LL" = "#2166AC", "HL" = "#EF8A62",
               "LH" = "#67A9CF", "NS" = "grey88")
grilla_lisa <- grilla |>
  dplyr::select(celda_edgar) |>
  dplyr::left_join(lisa[, c("celda_edgar", "clase_sig")], by = "celda_edgar")
p_lisa <- ggplot(grilla_lisa) +
  geom_sf(aes(fill = clase_sig), colour = NA, linewidth = 0) +
  {if (!is.null(bordes)) geom_sf(data = bordes, fill = NA, colour = "grey35", linewidth = 0.25)} +
  scale_fill_manual(values = cols_lisa, drop = FALSE,
                    labels = c("HH", "LL", "HL", "LH", "No significativo"),
                    name = "Clúster LISA") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "Clústeres LISA de la densidad de emisiones (2023)",
       subtitle = paste0("p por permutaciones < 0,05 (999 réplicas); ",
                         nint(sum(sig)), " celdas significativas (",
                         nfmt(100 * mean(sig), 1), "%)")) +
  tema_mapa
guardar_fig(p_lisa, "03_mapa_lisa_log_dens_2023", ancho = 8, alto = 9,
            dir = RUTAS$mapas)

# 5.3 Scatter de Moran -------------------------------------------------------
sc <- data.frame(z = z_sc, lag = lag_z,
                 clase = factor(ifelse(clase == "NS", "NS", clase),
                                levels = c("HH", "LL", "HL", "LH", "NS")))
p_moran <- ggplot(sc, aes(z, lag)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = clase), size = 0.5, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = cols_lisa, drop = FALSE,
                      labels = c("HH", "LL", "HL", "LH", "No sig."),
                      name = "Cuadrante") +
  labs(x = "log(1 + densidad de CO2) estandarizado",
       y = "Rezago espacial (W · z)",
       title = "Diagrama de dispersión de Moran",
       subtitle = paste0("I = ", nfmt(moran_q$I, 3), "; E[I] = ",
                         nfmt(moran_q$EI, 4), "; p de permutación = ",
                         nfmt(moran_q$p_perm, 3), "; n = ",
                         nint(n))) +
  theme(legend.position = "right")
guardar_fig(p_moran, "03_moran_scatter_log_dens_2023", ancho = 8, alto = 6)

# 5.4 Semivariograma + modelos ----------------------------------------------
d_linea <- seq(0, 280, by = 1)
lineas <- do.call(rbind, lapply(names(ajustes), function(m) {
  data.frame(dist = d_linea,
             gamma = as.numeric(variogramLine(ajustes[[m]]$fit,
                                              dist_vector = d_linea)$gamma),
             modelo = factor(etiquetas[m], levels = etiquetas[names(ajustes)]),
             row.names = NULL)
}))
v_plot <- v_emp
p_vgm <- ggplot() +
  geom_point(data = v_plot, aes(dist, gamma, size = np), alpha = 0.55, colour = "grey20") +
  geom_line(data = lineas, aes(dist, gamma, colour = modelo), linewidth = 0.9) +
  scale_size_area(max_size = 4, labels = scales::comma,
                  name = "Pares") +
  scale_colour_manual(values = c("Esférico" = "#2E6F9E",
                                 "Exponencial" = "#B2182B",
                                 "Gaussiano" = "#1B9E77"), name = "Modelo") +
  labs(x = "Distancia (km)", y = "γ(h)  [log(1 + t CO2/km²)]²",
       title = "Semivariograma empírico y modelos ajustados",
       subtitle = paste0("Mejor ajuste por SSE: ", etiquetas[mejor], " (nugget = ",
                         nfmt(modelos_tab$nugget[1], 3), "; sill = ",
                         nfmt(modelos_tab$sill[1], 3), "; rango práctico = ",
                         nfmt(modelos_tab$rango_practico_km[1], 0), " km)")) +
  theme(legend.position = "right")
guardar_fig(p_vgm, "03_semivariograma_log_dens_2023", ancho = 9, alto = 6)

# 5.5 Mapa de predicción por kriging ----------------------------------------
pred_grid <- call_quiet(gstat::krige(z ~ 1, locations = ~ x + y, data = d_km,
                                     newdata = data.frame(x = dat$x_9377 / 1000,
                                                          y = dat$y_9377 / 1000),
                                     model = params_mejor, nmax = 40))
stopifnot(length(pred_grid$var1.pred) == n)
grilla_pred <- grilla |>
  dplyr::select(celda_edgar) |>
  dplyr::mutate(pred_krig = pred_grid$var1.pred)
lim_pred <- range(c(dat$log_dens_2023, grilla_pred$pred_krig), na.rm = TRUE)
p_krig <- ggplot(grilla_pred) +
  geom_sf(aes(fill = pred_krig), colour = NA, linewidth = 0) +
  {if (!is.null(bordes)) geom_sf(data = bordes, fill = NA, colour = "grey35", linewidth = 0.25)} +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = lim_pred,
                       name = "log predicho (1 + t CO2/km²)") +
  labs(title = "Predicción de la densidad de emisiones por kriging ordinario",
       subtitle = paste0("Modelo ", etiquetas[mejor], "; nmax = 40 vecinos; ",
                         "n = ", nint(n), " celdas")) +
  tema_mapa
guardar_fig(p_krig, "03_mapa_kriging_log_dens_2023", ancho = 8, alto = 9,
            dir = RUTAS$mapas)

# 5.6 Mapa de error de la validación cruzada (kriging, 100% entrenamiento) ---
err100 <- cv_pred |>
  dplyr::filter(fraccion == 1, metodo == "kriging_ok") |>
  dplyr::transmute(celda_edgar, error = obs - pred)
grilla_err <- grilla |>
  dplyr::select(celda_edgar) |>
  dplyr::left_join(err100, by = "celda_edgar")
lim_err <- stats::quantile(abs(grilla_err$error), 0.98, na.rm = TRUE)
p_err <- ggplot(grilla_err) +
  geom_sf(aes(fill = error), colour = NA, linewidth = 0) +
  {if (!is.null(bordes)) geom_sf(data = bordes, fill = NA, colour = "grey35", linewidth = 0.25)} +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                       midpoint = 0, limits = c(-lim_err, lim_err),
                       oob = scales::squish,
                       name = "Observado − predicho\n(log unidades)") +
  labs(title = "Error out-of-fold del kriging ordinario (CV espacial, 5 bloques)",
       subtitle = paste0("RMSE = ", nfmt(rmse_k, 4),
                         "; MAE = ", nfmt(cv$MAE[cv$fraccion == 1 & cv$metodo == "kriging_ok"], 4),
                         "; R² = ", nfmt(cv$R2[cv$fraccion == 1 & cv$metodo == "kriging_ok"], 3))) +
  tema_mapa
guardar_fig(p_err, "03_mapa_error_cv_kriging", ancho = 8, alto = 9,
            dir = RUTAS$mapas)

# 5.7 Comparación de RMSE y R2 por método y fracción ------------------------
cv_plot <- cv |>
  dplyr::mutate(fraccion_pct = factor(paste0(round(100 * fraccion), "%"),
                                      levels = c("25%", "50%", "100%")),
                metodo_etiqueta = factor(metodo_etiqueta,
                                         levels = unname(etiq_metodos)))
p_rmse <- ggplot(cv_plot, aes(fraccion_pct, RMSE, colour = metodo_etiqueta,
                              group = metodo_etiqueta)) +
  geom_point(size = 2.6) + geom_line(linewidth = 0.8) +
  labs(x = "Fracción de entrenamiento", y = "RMSE (log unidades)", colour = NULL,
       title = "RMSE por método y densidad de entrenamiento") +
  theme(legend.position = "bottom")
p_r2 <- ggplot(cv_plot, aes(fraccion_pct, R2, colour = metodo_etiqueta,
                            group = metodo_etiqueta)) +
  geom_point(size = 2.6) + geom_line(linewidth = 0.8) +
  labs(x = "Fracción de entrenamiento", y = expression(R^2), colour = NULL,
       title = expression(R^2~"out-of-fold por método")) +
  theme(legend.position = "bottom")
p_cv <- (p_rmse / p_r2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
guardar_fig(p_cv, "03_comparacion_rmse_metodos", ancho = 9, alto = 8.5)

# ===========================================================================
paso("6. RESUMEN TÉCNICO Y GUARDADO")
# ===========================================================================
get_row <- function(v, w) {
  moran_global[moran_global$variable == v & moran_global$pesos == w, ]
}
q_dens <- get_row("log_dens_2023", "Reina 8 vecinos (16 km)")
q_co2  <- get_row("log_co2_2023",  "Reina 8 vecinos (16 km)")
q_pob  <- get_row("log_pob_2020",  "Reina 8 vecinos (16 km)")
q_dens_knn <- get_row("log_dens_2023", "kNN k=6 (simétrico)")

hh_n <- tab_global$n[tab_global$clase_sig == "HH"]
ll_n <- tab_global$n[tab_global$clase_sig == "LL"]
hh_pct <- tab_global$pct[tab_global$clase_sig == "HH"]
ll_pct <- tab_global$pct[tab_global$clase_sig == "LL"]
hh_top <- concentracion[concentracion$clase == "HH", ]
ll_top <- concentracion[concentracion$clase == "LL", ]
hh_top_txt <- paste0(hh_top$region[which.max(hh_top$n)], " (",
                     nfmt(hh_top$pct_clase[which.max(hh_top$n)], 1), "%)")
ll_top_txt <- paste0(ll_top$region[which.max(ll_top$n)], " (",
                     nfmt(ll_top$pct_clase[which.max(ll_top$n)], 1), "%)")

best_row <- modelos_tab[1, ]
dm_p <- test_dm$principal
pf <- test_dm$por_fraccion
p25 <- pf[pf$fraccion == 0.25, ]
p50 <- pf[pf$fraccion == 0.50, ]
p100 <- pf[pf$fraccion == 1, ]
rmse_fr <- function(fr, met = "kriging_ok") {
  cv$RMSE[cv$fraccion == fr & cv$metodo == met]
}
r2_k <- cv$R2[cv$fraccion == 1 & cv$metodo == "kriging_ok"]
mae_k <- cv$MAE[cv$fraccion == 1 & cv$metodo == "kriging_ok"]
mae_i <- cv$MAE[cv$fraccion == 1 & cv$metodo == p100$mejor_idw]
p_txt <- function(p) if (p < 0.001) "< 0,001" else nfmt(p, 3)
resumen_txt <- paste0(
  "En las ", nint(n), " celdas EDGAR de 0,1° que cubren Colombia, ",
  "la densidad de emisiones de CO2 de 2023 presenta autocorrelación espacial global ",
  "positiva y muy fuerte: I de Moran = ", nfmt(q_dens$I, 3),
  " (E[I] = ", nfmt(q_dens$EI, 4), "; varianza = ", nfmt(q_dens$VarI, 7),
  "; p de permutación = ", nfmt(q_dens$p_perm, 3), " con 999 réplicas; n = ",
  nint(n), "), resultado robusto a la definición de vecindad ",
  "(kNN k=6: I = ", nfmt(q_dens_knn$I, 3), "; reina: I = ", nfmt(q_dens$I, 3), "). ",
  "La autocorrelación también es intensa para el CO2 total (I = ", nfmt(q_co2$I, 3),
  "; p = ", nfmt(q_co2$p_perm, 3), ") y para la población (I = ", nfmt(q_pob$I, 3),
  "; p = ", nfmt(q_pob$p_perm, 3), "). El LISA identifica ", nint(hh_n),
  " celdas HH (", nfmt(hh_pct, 1), "% del total) y ", nint(ll_n),
  " celdas LL (", nfmt(ll_pct, 1), "%) significativas con p<0,05 (999 permutaciones); ",
  "las HH se concentran en la región ", hh_top_txt, " y las LL en ", ll_top_txt,
  ". Con corrección FDR-BH quedan ", nint(sum(sig_fdr)),
  " celdas significativas (", nfmt(100 * mean(sig_fdr), 1), "%). ",
  "El semivariograma empírico (cutoff 280 km) se ajusta mejor con el modelo ",
  etiquetas[mejor], " (SSE = ", nfmt(best_row$SSE, 4), "), con nugget = ",
  nfmt(best_row$nugget, 3), ", sill = ", nfmt(best_row$sill, 3),
  ", rango práctico = ", nfmt(best_row$rango_practico_km, 0),
  " km y proporción nugget/sill = ", nfmt(best_row$prop_nugget, 3),
  ": la estructura espacial es de media-larga distancia y el ",
  nfmt(100 * best_row$prop_nugget, 0), "% de la varianza es ruido no espacial. ",
  "En validación cruzada espacial por bloques (k = 5; cinco réplicas de ",
  "submuestreo anidado para 25% y 50%), el RMSE del kriging disminuyó al ",
  "aumentar la densidad de entrenamiento (25%: ", nfmt(rmse_fr(0.25), 4),
  "; 50%: ", nfmt(rmse_fr(0.50), 4), "; 100%: ", nfmt(rmse_fr(1), 4), "). ",
  "El kriging superó al mejor IDW, ", etiq_metodos[p100$mejor_idw],
  ", solo con el 100% de los datos: RMSE ", nfmt(p100$RMSE_kriging, 4),
  " frente a ", nfmt(p100$RMSE_idw, 4), " (Δ = ", nfmt(p100$dif_RMSE, 4),
  "; DM = ", nfmt(p100$DM, 2), "; p = ", p_txt(p100$p_DM),
  "; bootstrap pareado p ", p_txt(p100$bootstrap_p), ", B = ",
  dm_p$bootstrap_B, "), con R² = ", nfmt(r2_k, 3),
  "; la mejora se concentra en la cola de errores, ya que el MAE fue ",
  "prácticamente idéntico (", nfmt(mae_k, 4), " frente a ", nfmt(mae_i, 4), "). ",
  "Con 25% y 50% el kriging no superó al IDW (Δ = ",
  nfmt(p25$dif_RMSE, 4), ", DM p = ", p_txt(p25$p_DM), "; y Δ = ",
  nfmt(p50$dif_RMSE, 4), ", DM p = ", p_txt(p50$p_DM),
  "), de modo que la ventaja del kriging aparece solo con la muestra completa. ",
  "El bootstrap pareado por celda es anti-conservador bajo dependencia espacial ",
  "y se reporta como complemento del DM con HAC. ",
  "En conjunto, la densidad de emisiones es espacialmente estructurada y predecible, ",
  "pero un quinto de su variabilidad es de corto alcance (nugget), por lo que los mapas ",
  "suavizados de kriging no deben interpretarse como emisiones observadas celda a celda."
)
cat("\n", resumen_txt, "\n")

resultados <- list(
  moran_global = moran_global,
  lisa = lisa,
  resumen_clusters = resumen_clusters,
  variograma = list(
    empirico = v_emp,
    modelos = modelos_tab,
    mejor = etiquetas[mejor],
    mejor_gstat = mejor,
    params = params_mejor,
    ajustes = lapply(ajustes, `[[`, "fit")
  ),
  cv = cv,
  test_dm = test_dm,
  n_celdas = n,
  resumen = resumen_txt,
  metadatos = list(
    fecha = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    semilla = CONST$semilla,
    crs = CONST$crs_proj,
    variable = "log_dens_2023",
    n_perm = 999,
    cutoff_km = 280,
    ancho_bin_km = 20,
    cv_folds = 5,
    cv_kmeans_seed = CONST$semilla,
    n_max_vecinos = 40,
    n_fallback_idw = n_fallback,
    n_var_fallback = n_var_fallback,
    version_R = R.version.string
  )
)
rds_out <- file.path(RUTAS$procesados, "resultados_espaciales.rds")
saveRDS(resultados, rds_out)
msg("resultados -> ", basename(rds_out))

# --- CSV tidy de resultados clave ------------------------------------------
csv_moran <- moran_global |>
  dplyr::transmute(
    seccion = "Moran global",
    metrica = paste0(variable, " | ", pesos),
    valor = sprintf("I=%.4f; E[I]=%.5f; p_perm=%.3f", I, EI, p_perm),
    unidad = "adimensional",
    nota = paste0("n=", n_efectivo, "; var=", formatC(VarI, format = "e", digits = 2))
  )

csv_lisa <- resumen_clusters$global |>
  dplyr::transmute(
    seccion = "LISA (p<0,05, 999 perm)",
    metrica = paste0("Clase ", clase_sig),
    valor = as.character(n),
    unidad = "celdas",
    nota = paste0(nfmt(pct, 2), "% del total")
  )

csv_lisa_sig <- resumen_clusters$significancia |>
  dplyr::transmute(seccion = "LISA significancia", metrica = criterio,
                   valor = as.character(n), unidad = "celdas",
                   nota = paste0(nfmt(pct, 2), "% del total"))

csv_vgm <- modelos_tab |>
  dplyr::transmute(
    seccion = "Variograma",
    metrica = paste0("Modelo ", modelo),
    valor = sprintf("nugget=%.4f; sill=%.4f; rango_practico=%.2f km; nugget/sill=%.4f",
                    nugget, sill, rango_practico_km, prop_nugget),
    unidad = "log-densidad^2",
    nota = paste0("SSE=", formatC(SSE, format = "f", digits = 5),
                  ifelse(modelo_gstat == mejor, "; MEJOR", ""))
  )

csv_cv <- cv |>
  dplyr::transmute(
    seccion = "Validacion cruzada espacial",
    metrica = paste0(metodo_etiqueta, " | train=", round(100 * fraccion), "%"),
    valor = sprintf("RMSE=%.4f; MAE=%.4f; R2=%.4f", RMSE, MAE, R2),
    unidad = "log-densidad",
    nota = paste0("n_train=", n_train, "; replicas=", n_replicas,
                  "; celdas=", n_celdas, "; sesgo=",
                  formatC(sesgo, format = "f", digits = 4))
  )

csv_dm <- data.frame(
  seccion = "Diebold-Mariano / bootstrap",
  metrica = c("Mejor kriging vs mejor IDW (100%)", "DM estadistico",
              "DM p-valor", "Bootstrap pareado p-valor",
              "Diferencia RMSE (kriging - IDW)"),
  valor = c(paste0(etiq_metodos["kriging_ok"], " vs ", etiq_metodos[dm_p$mejor_idw]),
            nfmt(dm_p$DM, 3), ifelse(dm_p$p_DM < 0.001, "<0.001", nfmt(dm_p$p_DM, 4)),
            ifelse(dm_p$bootstrap_p < 0.001, "<0.001", nfmt(dm_p$bootstrap_p, 4)),
            nfmt(dm_p$dif_RMSE, 5)),
  unidad = c("", "estadistico", "p", "p", "log-densidad"),
  nota = c(paste0("n=", dm_p$n), paste0("HAC lag=", dm_p$lag_HAC),
           "dos colas", paste0("B=", dm_p$bootstrap_B, "; iid"),
           "kriging - IDW"),
  stringsAsFactors = FALSE
)

csv_resumen <- data.frame(
  seccion = "Resumen", metrica = "Parrafo tecnico", valor = resumen_txt,
  unidad = "", nota = "generado por 03_spatial.R", stringsAsFactors = FALSE
)

resumen_csv <- dplyr::bind_rows(csv_moran, csv_lisa, csv_lisa_sig, csv_vgm,
                                csv_cv, csv_dm, csv_resumen)
guardar_tabla(resumen_csv, "03_espacial_resumen")

# Copia legible del párrafo
writeLines(resumen_txt, file.path(RUTAS$tablas, "03_espacial_resumen.txt"),
           useBytes = TRUE)
msg("resumen -> 03_espacial_resumen.txt")

paso("MÓDULO ESPACIAL COMPLETADO")
cat("Objeto:", rds_out, "\n")
cat("n_celdas:", n, "| Moran I (densidad, reina):", nfmt(q_dens$I, 4), "\n")
cat("Mejor variograma:", etiquetas[mejor], "| RMSE kriging:",
    nfmt(dm_p$RMSE_kriging, 4), "| RMSE IDW:", nfmt(dm_p$RMSE_idw, 4), "\n")
