# ============================================================================
# 06_resumen.R — Consolida los resultados de los cuatro módulos en una tabla
# única y en un bloque de texto listo para el manuscrito.
# ============================================================================
source(file.path("scripts", "00_setup.R"))

leer <- function(f) {
  p <- file.path(RUTAS$procesados, f)
  if (!file.exists(p)) return(NULL)
  readRDS(p)
}
mu <- leer("resultados_muestreo.rds")
es <- leer("resultados_espaciales.rds")
do <- leer("resultados_doe.rds")
ml <- leer("resultados_ml.rds")
eda <- leer("eda_resumen.rds")

filas <- list()
agrega <- function(modulo, metrica, valor, unidad = "") {
  filas[[length(filas) + 1]] <<- data.frame(modulo = modulo, metrica = metrica,
                                            valor = as.character(valor), unidad = unidad)
}

# --- EDA ---
if (!is.null(eda)) {
  agrega("EDA", "CO2 1990", round(eda$serie$co2_owid[eda$serie$anio == 1990], 2), "Mt")
  agrega("EDA", "CO2 2019 (máximo)", round(eda$serie$co2_owid[eda$serie$anio == 2019], 2), "Mt")
  agrega("EDA", "CO2 2024", round(eda$serie$co2_owid[eda$serie$anio == 2024], 2), "Mt")
  agrega("EDA", "Concentración top 1% de celdas", round(100 * eda$concentracion["top1"], 1), "%")
  agrega("EDA", "Celdas del marco", eda$n_celdas, "celdas")
}
# --- Muestreo ---
if (!is.null(mu) && !is.null(mu$diseno_certeza)) {
  e <- mu$diseno_certeza$estimacion
  agrega("Muestreo", "Diseño con estrato de certeza: estimación", round(e["est_mt"], 2), "Mt")
  agrega("Muestreo", "Diseño con estrato de certeza: EE", round(e["se_mt"], 2), "Mt")
  agrega("Muestreo", "Diseño con estrato de certeza: error relativo", round(e["error_rel_pct"], 2), "%")
} else if (!is.null(mu)) {
  e <- mu$estimacion_final
  agrega("Muestreo", "Total real (censo EDGAR 2023)", round(e["real_mt"], 2), "Mt")
  agrega("Muestreo", "Total estimado (n=400 GRTS estrat.)", round(e["est_mt"], 2), "Mt")
  agrega("Muestreo", "EE del estimador", round(e["ee_mt"], 3), "Mt")
  agrega("Muestreo", "Error relativo", round(e["error_rel_pct"], 2), "%")
  agrega("Muestreo", "n para CV 5% (Neyman)", mu$n_para_cv5, "celdas")
  c400 <- mu$comparacion[mu$comparacion$n == 400, ]
  for (i in seq_len(nrow(c400))) {
    agrega("Muestreo", paste0("CV n=400: ", c400$diseño[i]), round(c400$cv[i], 2), "%")
  }
}
# --- Espacial ---
if (!is.null(es)) {
  m <- es$moran_global
  i_reina <- m$I[m$variable == "log_dens_2023" & grepl("Reina", m$pesos)][1]
  p_reina <- m$p_perm[m$variable == "log_dens_2023" & grepl("Reina", m$pesos)][1]
  agrega("Espacial", "Moran I densidad (reina)", round(i_reina, 3), "-")
  agrega("Espacial", "p de Moran (permutaciones)", p_reina, "-")
  cl <- es$resumen_clusters$global
  agrega("Espacial", "Celdas HH", cl$n[cl$clase_sig == "HH"], "celdas")
  agrega("Espacial", "Celdas LL", cl$n[cl$clase_sig == "LL"], "celdas")
  vm <- es$variograma$modelos[1, ]
  agrega("Espacial", "Variograma (modelo)", vm$modelo, "-")
  agrega("Espacial", "Nugget / sill", round(vm$prop_nugget, 3), "-")
  agrega("Espacial", "Rango práctico", round(vm$rango_practico_km, 1), "km")
  cv1 <- es$cv[es$cv$fraccion == 1, ]
  agrega("Espacial", "RMSE kriging (100%)", round(cv1$RMSE[cv1$metodo == "kriging_ok"], 3), "log")
  agrega("Espacial", "RMSE mejor IDW (100%)", round(min(cv1$RMSE[grepl("idw", cv1$metodo)]), 3), "log")
  agrega("Espacial", "p Diebold-Mariano", round(es$test_dm$principal$p_DM, 3), "-")
}
# --- DOE ---
if (!is.null(do)) {
  an <- do$anova
  getp <- function(term) {
    r <- an[an$termino == term, ]
    if (!nrow(r)) return(NA)
    r$p[1]
  }
  agrega("DOE", "F tecnología", round(an$f[an$termino == "tec"][1], 1), "-")
  agrega("DOE", "p tecnología", getp("tec"), "-")
  agrega("DOE", "eta2 parcial tecnología", round(an$eta2_parc[an$termino == "tec"][1], 3), "-")
  agrega("DOE", "F penetración", round(an$f[an$termino == "pen"][1], 1), "-")
  agrega("DOE", "p penetración", getp("pen"), "-")
  agrega("DOE", "F interacción", round(an$f[an$termino == "tec:pen"][1], 1), "-")
  agrega("DOE", "p interacción", getp("tec:pen"), "-")
  agrega("DOE", "Corridas", do$resumen$n_corridas %||% nrow(do$diseno), "-")
  if (!is.null(do$rsm$base_media_mt)) agrega("DOE", "Emisiones base 2030 (clima medio)", round(do$rsm$base_media_mt, 2), "Mt")
  if (!is.null(do$rsm$optimo)) {
    agrega("DOE", "Mezcla óptima solar", round(do$rsm$optimo[["x1"]], 1), "%")
    agrega("DOE", "Mezcla óptima eólica", round(do$rsm$optimo[["x2"]], 1), "%")
    agrega("DOE", "Emisiones mínimas 2030", round(do$rsm$emisiones_minimas_mt, 2), "Mt")
    agrega("DOE", "Umbral de FNCER útil", round(do$rsm$umbral_fncer_pct, 1), "pp")
    agrega("DOE", "Reducción máxima alcanzable", round(do$rsm$reduccion_maxima_mt, 2), "Mt")
  }
}
# --- ML ---
if (!is.null(ml)) {
  pr <- ml$pronostico
  agrega("ML", "Mejor modelo de pronóstico", pr$mejor_modelo, "-")
  agrega("ML", "RMSE validación (p.p.)", round(pr$sd_cv, 3), "p.p.")
  for (i in seq_len(nrow(pr$probabilidades))) {
    agrega("ML", paste0("CO2 2030 (mediana) ", pr$probabilidades$escenario[i]),
           round(pr$probabilidades$co2_2030_mediana[i], 2), "Mt")
    agrega("ML", paste0("P(CO2 2030 < 2024) ", pr$probabilidades$escenario[i]),
           round(pr$probabilidades$p_menor_2024[i], 3), "-")
    agrega("ML", paste0("P(CO2 2030 < 2015) ", pr$probabilidades$escenario[i]),
           round(pr$probabilidades$p_menor_2015[i], 3), "-")
  }
  agrega("ML", "Moran residual RF", round(ml$celdas$moran_residual[1], 3), "-")
  agrega("ML", "RMSE espacial RF", round(ml$celdas$metricas$RMSE[ml$celdas$metricas$modelo == "RF"], 3), "log")
  agrega("ML", "RMSE aleatorio RF", round(ml$celdas$cv_espacial$RMSE[1], 3), "log")
}

resumen <- do.call(rbind, filas)
guardar_tabla(resumen, "06_resumen_global")

cat("\n", strrep("=", 72), "\nRESUMEN GLOBAL PARA EL PAPER\n", strrep("=", 72), "\n", sep = "")
print(resumen, row.names = FALSE)
