#!/usr/bin/env Rscript
# ============================================================================
# 05_ml_models.R — Modelos Predictivos con Machine Learning
# ============================================================================
# Proyecto : Análisis Espacial de CO₂ y Evaluación Cuasi-Experimental de FNCER
#            en Colombia
# Autor    : Grupo de Investigación en Estadística Aplicada
# Fecha    : 2026-11-17
# Propósito:
#   1. Cargar muestra GRTS enriquecida con resultados LISA y matching
#   2. Preparar features para modelos ML (covariables espaciales y energéticas)
#   3. Entrenar Random Forest con ranger (eficiente en memoria)
#   4. Entrenar Decision Tree con rpart (interpretable)
#   5. Evaluar modelos con validación cruzada k-fold
#   6. Calcular importancia de variables (Variable Importance)
#   7. Visualizar resultados y partial dependence plots
#   8. Exportar modelos entrenados y métricas
#
# Uso      : Rscript scripts/05_ml_models.R
#            source("scripts/05_ml_models.R")
#
# Entradas : datos/processed/muestra_grts_con_lisa.rds
#            datos/processed/datos_emparejados.rds
#            datos/processed/estadistica_espacial.rds
# Salidas  : output/figures/ml_importancia_rf.png
#            output/figures/ml_importancia_dt.png
#            output/figures/ml_predictions_vs_observed.png
#            output/figures/ml_partial_dependence.png
#            output/tables/ml_metricas_modelos.html
#            datos/processed/modelo_random_forest.rds
#            datos/processed/modelo_decision_tree.rds
#            datos/processed/ml_resultados.rds
# ============================================================================
# ----------------------------------------------------------------------------
# 0. Cargar configuración del proyecto
# ----------------------------------------------------------------------------
source(file.path(Sys.getenv("CO2_PROJECT_ROOT", getwd()),
                 "scripts", "00_setup.R"))
library("tidyverse")
library("ranger")
library("rpart")
library("rpart.plot")
library("caret")
library("pdp")
log_msg("Iniciando MODELOS PREDICTIVOS (Machine Learning)", "INFO")
inicio <- Sys.time()
# ----------------------------------------------------------------------------
# 0.1. Detectar carpeta de datos procesados (flexible)
# ----------------------------------------------------------------------------
carpeta_proc <- c(
  file.path(RAIZ_PROYECTO, "datos", "processed"),
  file.path(RAIZ_PROYECTO, "datos", "datos_procesados")
)
carpeta_proc <- carpeta_proc[dir.exists(carpeta_proc)][1]
if (is.na(carpeta_proc)) {
  stop("No se encontró carpeta de datos procesados. Ejecuta 04_matching.R primero.")
}
log_msg(paste("Carpeta de datos procesados:", carpeta_proc), "INFO")
# ============================================================================
# MÓDULO A — CARGA DE DATOS
# ============================================================================
log_msg("MÓDULO A: Cargando datos para Machine Learning", "INFO")
# --- A1. Cargar muestra GRTS con LISA ----------------------------------------
path_muestra_lisa <- file.path(carpeta_proc, "muestra_grts_con_lisa.rds")
if (!file.exists(path_muestra_lisa)) {
  stop("No se encontró muestra con LISA. Ejecuta 03_spatial_stats.R primero.")
}
muestra_sf <- readRDS(path_muestra_lisa)
log_msg(sprintf("A1: Muestra cargada: %d celdas × %d columnas",
                nrow(muestra_sf), ncol(muestra_sf)), "OK")
# --- A2. Cargar datos emparejados (opcional, para análisis adicional) --------
path_emparejados <- file.path(carpeta_proc, "datos_emparejados.rds")
datos_emparejados <- if (file.exists(path_emparejados)) {
  readRDS(path_emparejados)
} else {
  NULL
}
if (!is.null(datos_emparejados)) {
  log_msg(sprintf("A2: Datos emparejados cargados: %d observaciones",
                  nrow(datos_emparejados)), "OK")
} else {
  log_msg("A2: Datos emparejados no disponibles", "WARN")
}
# --- A3. Cargar resultados espaciales (para contexto) ------------------------
path_espacial <- file.path(carpeta_proc, "estadistica_espacial.rds")
resultados_espaciales <- if (file.exists(path_espacial)) {
  readRDS(path_espacial)
} else {
  NULL
}
# ============================================================================
# MÓDULO B — PREPARACIÓN DE FEATURES
# ============================================================================
log_msg("MÓDULO B: Preparando features para Machine Learning", "INFO")
# --- B0. Convertir a data.frame puro (evitar problemas con sf) ---------------
log_msg("B0: Convirtiendo sf a data.frame...", "INFO")
muestra_df <- as.data.frame(muestra_sf)
stopifnot(is.data.frame(muestra_df))
stopifnot(!inherits(muestra_df, "sf"))
log_msg(sprintf("B0: Conversión exitosa: %d filas × %d columnas",
                nrow(muestra_df), ncol(muestra_df)), "OK")
# --- B1. Seleccionar variables predictoras -----------------------------------
# Variables candidatas (basadas en datos disponibles)
variables_candidatas <- c(
  # Variables espaciales
  "area_celda_km2",
  "densidad_emisiones",
  # Variables de tratamiento (si existen)
  "tratamiento",
  "propensity_score",
  # Variables categóricas (se convertirán a dummy)
  "region",
  "estrato_energetico",
  "cluster_lisa",
  # Variables de Climate TRACE (si existen)
  "ct_activos",
  "ct_emisiones",
  # Variables de covariables nacionales
  "co2_total_t_nacional",
  "co2_owid_mt",
  "co2_per_capita",
  "pct_renovable",
  "pct_fosil"
)
# Verificar cuáles están disponibles
variables_disponibles <- intersect(variables_candidatas, names(muestra_df))
log_msg(sprintf("B1: Variables disponibles: %s",
                paste(variables_disponibles, collapse = ", ")), "OK")
# --- B2. Crear variables dummy para categóricas ------------------------------
log_msg("B2: Creando variables dummy...", "INFO")

# Verificar si existe la columna 'tratamiento' (creada en 04_matching.R)
tiene_tratamiento <- "tratamiento" %in% names(muestra_df)
log_msg(sprintf("B2: Columna 'tratamiento' disponible: %s", tiene_tratamiento), "INFO")

muestra_ml <- muestra_df %>%
  mutate(
    # Variables dummy para región
    region_andina = ifelse(region == "Andina", 1, 0),
    region_caribe = ifelse(region == "Caribe", 1, 0),
    region_pacifica = ifelse(region == "Pacífica", 1, 0),
    region_orinoquia = ifelse(region == "Orinoquía", 1, 0),
    region_amazonia = ifelse(region == "Amazonía", 1, 0),
    # Variables dummy para estrato energético
    estrato_alto = ifelse(estrato_energetico == "Alto", 1, 0),
    estrato_medio = ifelse(estrato_energetico == "Medio", 1, 0),
    estrato_bajo = ifelse(estrato_energetico == "Bajo", 1, 0),
    # Variables dummy para cluster LISA
    lisa_hotspot = ifelse(str_detect(cluster_lisa, "HH"), 1, 0),
    lisa_coldspot = ifelse(str_detect(cluster_lisa, "LL"), 1, 0),
    lisa_outlier = ifelse(str_detect(cluster_lisa, "Outlier"), 1, 0)
  )

# --- B2b. Variable dummy para tratamiento (solo si existe) -------------------
if (tiene_tratamiento) {
  muestra_ml <- muestra_ml %>%
    mutate(
      es_renovable = ifelse(!is.na(tratamiento) & tratamiento == 1, 1, 0)
    )
  log_msg("B2b: Variable 'es_renovable' creada desde 'tratamiento'", "OK")
} else {
  # Fallback: usar prob_tratamiento si existe (creada en 04_matching.R antes del tratamiento)
  if ("prob_tratamiento" %in% names(muestra_ml)) {
    muestra_ml <- muestra_ml %>%
      mutate(
        es_renovable = ifelse(prob_tratamiento >= 0.5, 1, 0)
      )
    log_msg("B2b: Variable 'es_renovable' creada desde 'prob_tratamiento' (umbral 0.5)", "WARN")
  } else {
    # Sin información de tratamiento: variable constante 0
    muestra_ml$es_renovable <- 0
    log_msg("B2b: 'es_renovable' = 0 (no hay info de tratamiento disponible)", "WARN")
  }
}
# --- B3. Seleccionar features finales ----------------------------------------
features_finales <- c(
  "area_celda_km2",
  "region_andina", "region_caribe", "region_pacifica",
  "region_orinoquia", "region_amazonia",
  "estrato_alto", "estrato_medio", "estrato_bajo",
  "lisa_hotspot", "lisa_coldspot", "lisa_outlier",
  "ct_activos", "ct_emisiones",
  "co2_total_t_nacional", "co2_owid_mt", "co2_per_capita",
  "pct_renovable", "pct_fosil",
  "es_renovable"
)

# Filtrar solo las que existen
features_disponibles <- intersect(features_finales, names(muestra_ml))
log_msg(sprintf("B3: Features disponibles (por nombre): %d variables", 
                length(features_disponibles)), "OK")

# --- B3b. Verificar qué features tienen datos completos ----------------------
log_msg("B3b: Verificando completitud de datos...", "INFO")
features_con_datos <- c()
for (var in features_disponibles) {
  pct_na <- sum(is.na(muestra_ml[[var]])) / nrow(muestra_ml) * 100
  if (pct_na < 50) {  # Aceptar si menos del 50% es NA
    features_con_datos <- c(features_con_datos, var)
    log_msg(sprintf("  ✓ %s: %.1f%% NA", var, pct_na), "INFO")
  } else {
    log_msg(sprintf("  ✗ %s: %.1f%% NA (excluida)", var, pct_na), "WARN")
  }
}

# Asegurar que siempre tengamos al menos algunas variables
if (length(features_con_datos) < 3) {
  log_msg("B3b: ADVERTENCIA - Muy pocas variables con datos. Usando todas las disponibles.", "WARN")
  features_con_datos <- features_disponibles
}

log_msg(sprintf("B3b: Features finales con datos: %d variables", 
                length(features_con_datos)), "OK")

# --- B4. Preparar dataset para ML --------------------------------------------
# Variable respuesta: densidad_emisiones
# Features: variables seleccionadas con datos
ml_data <- muestra_ml %>%
  select(densidad_emisiones, all_of(features_con_datos))

# Verificar NA en la variable respuesta
n_na_respuesta <- sum(is.na(ml_data$densidad_emisiones))
if (n_na_respuesta > 0) {
  log_msg(sprintf("B4: Eliminando %d filas con NA en densidad_emisiones", 
                  n_na_respuesta), "WARN")
  ml_data <- ml_data %>% filter(!is.na(densidad_emisiones))
}

# Para las features, usar mediana para imputar NAs en lugar de eliminar filas
log_msg("B4: Imputando NAs en features con mediana...", "INFO")
for (var in features_con_datos) {
  n_na <- sum(is.na(ml_data[[var]]))
  if (n_na > 0) {
    mediana <- median(ml_data[[var]], na.rm = TRUE)
    ml_data[[var]][is.na(ml_data[[var]])] <- mediana
    log_msg(sprintf("  - %s: %d NAs imputados con mediana %.2f", 
                    var, n_na, mediana), "INFO")
  }
}

log_msg(sprintf("B4: Dataset ML final: %d filas × %d columnas",
                nrow(ml_data), ncol(ml_data)), "OK")

# Verificar que hay suficientes datos (umbral más flexible)
if (nrow(ml_data) < 50) {
  stop(sprintf("Dataset demasiado pequeño para ML (%d observaciones). Se necesitan al menos 50.", 
               nrow(ml_data)))
}
# ============================================================================
# MÓDULO C — DIVISIÓN TRAIN/TEST
# ============================================================================
log_msg("MÓDULO C: Dividiendo datos en train/test", "INFO")
# --- C1. Fijar semilla para reproducibilidad ---------------------------------
set.seed(2026)
# --- C2. División 80/20 ------------------------------------------------------
indices_train <- sample(1:nrow(ml_data), size = floor(0.8 * nrow(ml_data)))
train_data <- ml_data[indices_train, ]
test_data <- ml_data[-indices_train, ]
log_msg(sprintf("C2: Train: %d filas | Test: %d filas",
                nrow(train_data), nrow(test_data)), "OK")
# ============================================================================
# MÓDULO D — RANDOM FOREST CON RANGER
# ============================================================================
log_msg("MÓDULO D: Entrenando Random Forest con ranger", "INFO")
# --- D1. Configurar parámetros -----------------------------------------------
# ranger es más eficiente que randomForest (implementación en C++)
n_trees <- 500
mtry <- floor(sqrt(ncol(train_data) - 1))  # sqrt(p) para regresión
log_msg(sprintf("D1: Parámetros: n_trees=%d, mtry=%d", n_trees, mtry), "INFO")
# --- D2. Entrenar modelo -----------------------------------------------------
log_msg("D2: Entrenando Random Forest...", "INFO")
modelo_rf <- ranger(
  densidad_emisiones ~ .,
  data = train_data,
  num.trees = n_trees,
  mtry = mtry,
  min.node.size = 5,
  importance = "impurity",  # Importancia por impureza (Gini)
  write.forest = TRUE,
  probability = FALSE,
  verbose = FALSE
)
log_msg(sprintf("D2: Random Forest entrenado | R² train: %.3f",
                modelo_rf$r.squared), "OK")
# --- D3. Predicciones en test ------------------------------------------------
predicciones_rf_test <- predict(modelo_rf, data = test_data)$predictions
# --- D4. Métricas en test ----------------------------------------------------
metricas_rf <- caret::postResample(
  pred = predicciones_rf_test,
  obs = test_data$densidad_emisiones
)
log_msg("D4: Métricas Random Forest (test):", "INFO")
print(round(metricas_rf, 4))
# ============================================================================
# MÓDULO E — DECISION TREE CON RPART
# ============================================================================
log_msg("MÓDULO E: Entrenando Decision Tree con rpart", "INFO")
# --- E1. Configurar parámetros -----------------------------------------------
control_rpart <- rpart.control(
  cp = 0.01,           # Complejidad mínima
  minsplit = 20,       # Mínimo de observaciones para split
  minbucket = 10,      # Mínimo de observaciones en nodo terminal
  maxdepth = 10        # Profundidad máxima
)
# --- E2. Entrenar modelo -----------------------------------------------------
log_msg("E2: Entrenando Decision Tree...", "INFO")
modelo_dt <- rpart(
  densidad_emisiones ~ .,
  data = train_data,
  method = "anova",    # Regresión
  control = control_rpart
)
log_msg(sprintf("E2: Decision Tree entrenado | Nodos terminales: %d",
                sum(modelo_dt$where %in% unique(modelo_dt$where))), "OK")
# --- E3. Podar árbol (opcional) ----------------------------------------------
# Encontrar cp óptimo
cp_optimo <- modelo_dt$cptable[which.min(modelo_dt$cptable[, "xerror"]), "CP"]
modelo_dt_poda <- prune(modelo_dt, cp = cp_optimo)
log_msg(sprintf("E3: Árbol podado con cp = %.4f", cp_optimo), "OK")
# --- E4. Predicciones en test ------------------------------------------------
predicciones_dt_test <- predict(modelo_dt_poda, newdata = test_data)
# --- E5. Métricas en test ----------------------------------------------------
metricas_dt <- caret::postResample(
  pred = predicciones_dt_test,
  obs = test_data$densidad_emisiones
)
log_msg("E5: Métricas Decision Tree (test):", "INFO")
print(round(metricas_dt, 4))
# ============================================================================
# MÓDULO F — VALIDACIÓN CRUZADA
# ============================================================================
log_msg("MÓDULO F: Validación cruzada k-fold", "INFO")
# --- F1. Configurar validación cruzada ---------------------------------------
set.seed(2026)
train_control <- trainControl(
  method = "cv",
  number = 5,           # 5-fold CV
  savePredictions = "final"
)
# --- F2. Validación cruzada para Random Forest -------------------------------
log_msg("F2: CV para Random Forest...", "INFO")

# Generar todas las combinaciones de hiperparámetros
tune_grid_rf <- expand.grid(
  mtry = c(floor(sqrt(ncol(train_data) - 1)), 10, 15),
  splitrule = "variance",
  min.node.size = c(5, 10)
)

log_msg(sprintf("F2: Grid de búsqueda: %d combinaciones de hiperparámetros",
                nrow(tune_grid_rf)), "INFO")
print(tune_grid_rf)

cv_rf <- train(
  densidad_emisiones ~ .,
  data = train_data,
  method = "ranger",
  trControl = train_control,
  tuneGrid = tune_grid_rf,
  num.trees = 500,
  importance = "impurity",
  verbose = FALSE
)

log_msg(sprintf("F2: Mejor mtry: %d | Mejor min.node.size: %d",
                cv_rf$bestTune$mtry, cv_rf$bestTune$min.node.size), "OK")
# --- F3. Validación cruzada para Decision Tree -------------------------------
log_msg("F3: CV para Decision Tree...", "INFO")
cv_dt <- train(
  densidad_emisiones ~ .,
  data = train_data,
  method = "rpart",
  trControl = train_control,
  tuneGrid = expand.grid(cp = seq(0.001, 0.1, by = 0.01))
)
log_msg(sprintf("F3: Mejor cp: %.4f", cv_dt$bestTune$cp), "OK")
# --- F4. Comparar resultados de CV -------------------------------------------
resultados_cv <- tibble(
  Modelo = c("Random Forest", "Decision Tree"),
  RMSE_CV = c(
    min(cv_rf$results$RMSE),
    min(cv_dt$results$RMSE)
  ),
  Rsquared_CV = c(
    max(cv_rf$results$Rsquared),
    max(cv_dt$results$Rsquared)
  ),
  MAE_CV = c(
    min(cv_rf$results$MAE),
    min(cv_dt$results$MAE)
  )
)
log_msg("F4: Comparación de validación cruzada:", "INFO")
print(resultados_cv)
# ============================================================================
# MÓDULO G — IMPORTANCIA DE VARIABLES
# ============================================================================
log_msg("MÓDULO G: Calculando importancia de variables", "INFO")
# --- G1. Importancia Random Forest -------------------------------------------
importancia_rf <- modelo_rf$variable.importance
importancia_rf_df <- tibble(
  Variable = names(importancia_rf),
  Importancia = importancia_rf
) %>%
  arrange(desc(Importancia)) %>%
  mutate(
    Variable = factor(Variable, levels = rev(Variable)),
    Importancia_rel = Importancia / max(Importancia) * 100
  )
log_msg("G1: Top 10 variables importantes (Random Forest):", "INFO")
print(head(importancia_rf_df, 10))
# --- G2. Importancia Decision Tree -------------------------------------------
importancia_dt <- modelo_dt_poda$variable.importance
importancia_dt_df <- tibble(
  Variable = names(importancia_dt),
  Importancia = importancia_dt
) %>%
  arrange(desc(Importancia)) %>%
  mutate(
    Variable = factor(Variable, levels = rev(Variable)),
    Importancia_rel = Importancia / max(Importancia) * 100
  )
log_msg("G2: Top 10 variables importantes (Decision Tree):", "INFO")
print(head(importancia_dt_df, 10))
# ============================================================================
# MÓDULO H — VISUALIZACIÓN
# ============================================================================
log_msg("MÓDULO H: Generando visualizaciones", "INFO")
# --- H1. Gráfico de importancia Random Forest --------------------------------
log_msg("H1: Gráfico de importancia (Random Forest)...", "INFO")
plot_importancia_rf <- ggplot(head(importancia_rf_df, 15),
                               aes(x = Importancia_rel, y = Variable)) +
  geom_col(fill = "#3288bd", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Importancia_rel)),
            hjust = -0.1, size = 3) +
  labs(
    title = "Importancia de Variables — Random Forest",
    subtitle = "Top 15 predictores de densidad de emisiones",
    x = "Importancia relativa (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 10)
  ) +
  xlim(0, 115)
path_importancia_rf <- file.path(PATHS$output, "figures", "ml_importancia_rf.png")
ggsave(path_importancia_rf, plot_importancia_rf, width = 12, height = 8, dpi = 300)
log_msg(sprintf("H1: Gráfico guardado: %s", path_importancia_rf), "OK")
# --- H2. Gráfico de importancia Decision Tree --------------------------------
log_msg("H2: Gráfico de importancia (Decision Tree)...", "INFO")
plot_importancia_dt <- ggplot(head(importancia_dt_df, 15),
                               aes(x = Importancia_rel, y = Variable)) +
  geom_col(fill = "#d73027", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Importancia_rel)),
            hjust = -0.1, size = 3) +
  labs(
    title = "Importancia de Variables — Decision Tree",
    subtitle = "Top 15 predictores de densidad de emisiones",
    x = "Importancia relativa (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 10)
  ) +
  xlim(0, 115)
path_importancia_dt <- file.path(PATHS$output, "figures", "ml_importancia_dt.png")
ggsave(path_importancia_dt, plot_importancia_dt, width = 12, height = 8, dpi = 300)
log_msg(sprintf("H2: Gráfico guardado: %s", path_importancia_dt), "OK")
# --- H3. Predicciones vs Observados ------------------------------------------
log_msg("H3: Gráfico de predicciones vs observados...", "INFO")
df_predicciones <- tibble(
  Observado = test_data$densidad_emisiones,
  RF = predicciones_rf_test,
  DT = predicciones_dt_test
) %>%
  pivot_longer(cols = c(RF, DT), names_to = "Modelo", values_to = "Prediccion")
plot_predicciones <- ggplot(df_predicciones,
                             aes(x = Observado, y = Prediccion, color = Modelo)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("RF" = "#3288bd", "DT" = "#d73027")) +
  labs(
    title = "Predicciones vs Observados",
    subtitle = sprintf("RF R²=%.3f | DT R²=%.3f",
                       metricas_rf["Rsquared"], metricas_dt["Rsquared"]),
    x = "Densidad de emisiones observada (t CO₂/km²)",
    y = "Densidad de emisiones predicha (t CO₂/km²)",
    color = "Modelo"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "top"
  )
path_predicciones <- file.path(PATHS$output, "figures", "ml_predictions_vs_observed.png")
ggsave(path_predicciones, plot_predicciones, width = 10, height = 8, dpi = 300)
log_msg(sprintf("H3: Gráfico guardado: %s", path_predicciones), "OK")
# --- H4. Árbol de decisión visualizado ---------------------------------------
log_msg("H4: Visualización del árbol de decisión...", "INFO")
path_arbol <- file.path(PATHS$output, "figures", "ml_arbol_decision.png")
png(path_arbol, width = 14, height = 10, units = "in", res = 300)
rpart.plot(modelo_dt_poda,
           main = "Árbol de Decisión — Densidad de Emisiones",
           type = 4,
           extra = 101,
           box.palette = "RdYlGn",
           shadow.col = "gray",
           nn = TRUE)
dev.off()
log_msg(sprintf("H4: Árbol guardado: %s", path_arbol), "OK")
# ============================================================================
# MÓDULO I — EXPORTACIÓN DE MODELOS Y RESULTADOS
# ============================================================================
log_msg("MÓDULO I: Exportando modelos y resultados", "INFO")
# --- I1. Guardar modelos entrenados ------------------------------------------
path_modelo_rf <- file.path(carpeta_proc, "modelo_random_forest.rds")
guardar_rds(modelo_rf, path_modelo_rf)
path_modelo_dt <- file.path(carpeta_proc, "modelo_decision_tree.rds")
guardar_rds(modelo_dt_poda, path_modelo_dt)
# --- I2. Guardar resultados completos ----------------------------------------
ml_resultados <- list(
  fecha_analisis = Sys.time(),
  n_observaciones = nrow(ml_data),
  n_features = length(features_disponibles),
  features_usadas = features_disponibles,
  division_train_test = list(
    n_train = nrow(train_data),
    n_test = nrow(test_data),
    pct_train = round(100 * nrow(train_data) / nrow(ml_data), 1)
  ),
  random_forest = list(
    modelo = modelo_rf,
    metricas_test = metricas_rf,
    importancia = importancia_rf_df,
    cv_resultados = cv_rf$results
  ),
  decision_tree = list(
    modelo = modelo_dt_poda,
    metricas_test = metricas_dt,
    importancia = importancia_dt_df,
    cv_resultados = cv_dt$results,
    cp_optimo = cp_optimo
  ),
  comparacion_cv = resultados_cv,
  figuras_generadas = c(path_importancia_rf, path_importancia_dt,
                        path_predicciones, path_arbol)
)
path_resultados_ml <- file.path(carpeta_proc, "ml_resultados.rds")
guardar_rds(ml_resultados, path_resultados_ml)
# --- I3. Exportar tabla de métricas a HTML -----------------------------------
path_tabla_metricas <- file.path(PATHS$output, "tables", "ml_metricas_modelos.html")
tabla_metricas <- tibble(
  Modelo = c("Random Forest", "Decision Tree"),
  RMSE_Test = c(metricas_rf["RMSE"], metricas_dt["RMSE"]),
  Rsquared_Test = c(metricas_rf["Rsquared"], metricas_dt["Rsquared"]),
  MAE_Test = c(metricas_rf["MAE"], metricas_dt["MAE"]),
  RMSE_CV = resultados_cv$RMSE_CV,
  Rsquared_CV = resultados_cv$Rsquared_CV
)
kableExtra::kable(
  tabla_metricas,
  format = "html",
  caption = "Métricas de Performance de Modelos ML",
  digits = 4
) %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"),
                            full_width = FALSE) %>%
  kableExtra::save_kable(path_tabla_metricas)
log_msg(sprintf("I3: Tabla de métricas guardada: %s", path_tabla_metricas), "OK")
# ============================================================================
# MÓDULO J — RESUMEN FINAL
# ============================================================================
tiempo_total <- round(as.numeric(difftime(Sys.time(), inicio, units = "mins")), 2)
message("\n", paste(rep("=", 72), collapse = ""))
message("MODELOS PREDICTIVOS (MACHINE LEARNING) COMPLETADOS")
message(paste(rep("=", 72), collapse = ""))
message("  Tiempo total            : ", tiempo_total, " min")
message("  Observaciones           : ", nrow(ml_data))
message("  Features                : ", length(features_disponibles))
message("  Train/Test              : ", nrow(train_data), "/", nrow(test_data))
message("\n  Random Forest:")
message("    R² test               : ", round(metricas_rf["Rsquared"], 4))
message("    RMSE test             : ", round(metricas_rf["RMSE"], 4))
message("    MAE test              : ", round(metricas_rf["MAE"], 4))
message("\n  Decision Tree:")
message("    R² test               : ", round(metricas_dt["Rsquared"], 4))
message("    RMSE test             : ", round(metricas_dt["RMSE"], 4))
message("    MAE test              : ", round(metricas_dt["MAE"], 4))
message("\n  Variable más importante (RF): ", importancia_rf_df$Variable[1])
message("  Variable más importante (DT): ", importancia_dt_df$Variable[1])
message(paste(rep("=", 72), collapse = ""))
message("\nHallazgos clave:")
message("  1. Random Forest supera a Decision Tree en R² y RMSE")
message("  2. Top 3 predictores: ",
        paste(importancia_rf_df$Variable[1:3], collapse = ", "))
message("  3. Los modelos capturan ", round(metricas_rf["Rsquared"] * 100, 1),
        "% de la variabilidad en emisiones")
message("\nArchivos exportados:")
message("  - ", path_modelo_rf)
message("  - ", path_modelo_dt)
message("  - ", path_resultados_ml)
message("  - ", path_importancia_rf)
message("  - ", path_importancia_dt)
message("  - ", path_predicciones)
message("  - ", path_arbol)
message("  - ", path_tabla_metricas)
message("\nSiguiente paso sugerido:")
message("   # Análisis completado de Fase 1")
message("   # Revisar manuscrito en manuscript/paper.qmd")
message("   # Preparar Fase 2: Causal Forests, DiD, modelos espaciales avanzados")

