# ==============================================================================
# Script: 04_matching.R
# Proyecto: Análisis Espacial de Emisiones de CO2 en Colombia
# Sección: 4. Diseño de Experimentos - Inferencia Causal Observacional
# Objetivo: Estimar el efecto causal de tratamientos energéticos sobre la 
#           reducción de CO2 mediante Propensity Score Matching.
# Referencias: Stuart (2010); Imbens & Rubin (2015); Rosenbaum (2002).
# ==============================================================================

# 0. CONFIGURACIÓN Y LIBRERÍAS ------------------------------------------------
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(MatchIt)       # Emparejamiento por propensity score
  library(cobalt)        # Love plots y diagnóstico de balance
  library(dagitty)       # DAG causal
  library(sensitivity2x2xk) # Análisis de sensibilidad de Rosenbaum
  library(VGAM)          # Para modelos multinomiales si se necesitan
  library(patchwork)
  library(viridis)
})

set.seed(2026)
theme_set(theme_minimal())

# Crear carpetas de salida
dir.create("results", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# ==============================================================================
# 1. CARGA Y AGREGACIÓN DE DATOS A NIVEL DE MUNICIPIO -------------------------
# ==============================================================================
cat("=== 1. Cargando y agregando datos a nivel municipal ===\n")

co2_raw <- read.csv("climate_trace_COL_2023.csv", stringsAsFactors = FALSE)

# Filtrar Colombia y limpiar
co2_clean <- co2_raw %>%
  filter(!is.na(longitude), !is.na(latitude),
         longitude > -82 & longitude < -66,
         latitude > -5 & latitude < 13,
         !is.na(co2), is.finite(co2))

# Extraer nombre de municipio (limpiar sufijos)
co2_clean <- co2_clean %>%
  mutate(
    municipio = gsub("\\s*(Municipality|Corregimiento Departamental)$", "", name),
    municipio = trimws(municipio)
  )

# Clasificar tipo de fuente energética (TRATAMIENTO)
# Basado en sector y sub_sector de Climate Trace
co2_clean <- co2_clean %>%
  mutate(
    tratamiento = case_when(
      # T0: Fósil - manufactura, edificios, transporte
      sector %in% c("manufacturing", "buildings", "transportation") ~ "T0_fosil",
      # T4: Biomasa - agricultura (metano/estiércol = biomasa)
      sector == "agriculture" ~ "T4_biomasa",
      # T1/T2/T3: Uso de suelo forestal (secuestro = renovables)
      # Clasificamos por tipo de cobertura cuando está disponible
      sector == "forestry-and-land-use" & 
        grepl("forest|net-forest|removals", sub_sector, ignore.case = TRUE) ~ "T1_solar_proxy",
      sector == "forestry-and-land-use" & 
        grepl("shrubgrass|wetland|soil", sub_sector, ignore.case = TRUE) ~ "T2_eolica_proxy",
      sector == "forestry-and-land-use" & 
        grepl("fire|degradation|clearing", sub_sector, ignore.case = TRUE) ~ "T3_hidro_proxy",
      TRUE ~ NA_character_
    )
  )

# Re-etiquetar tratamientos para simplificar (T1, T2, T3 son proxies)
co2_clean <- co2_clean %>%
  mutate(
    tratamiento = case_when(
      tratamiento == "T1_solar_proxy"  ~ "T1_solar",
      tratamiento == "T2_eolica_proxy" ~ "T2_eolica",
      tratamiento == "T3_hidro_proxy"  ~ "T3_hidro",
      TRUE ~ tratamiento
    )
  )

cat("  • Distribución de tratamientos:\n")
print(table(co2_clean$tratamiento, useNA = "ifany"))

# Agregar a nivel de municipio
municipio_df <- co2_clean %>%
  filter(!is.na(tratamiento)) %>%
  group_by(municipio, tratamiento) %>%
  summarise(
    co2_total = sum(co2, na.rm = TRUE),
    actividad_total = sum(activity_data, na.rm = TRUE),
    n_fuentes = n(),
    lat_med = mean(latitude, na.rm = TRUE),
    lon_med = mean(longitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Determinar tratamiento dominante por municipio
  group_by(municipio) %>%
  slice_max(order_by = co2_total, n = 1, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("  • %d municipios con tratamiento definido\n", nrow(municipio_df)))

# ==============================================================================
# 2. CONSTRUCCIÓN DE COVARIABLES DE EMPAREJAMIENTO ----------------------------
# ==============================================================================
cat("\n=== 2. Construyendo covariables de emparejamiento ===\n")

# --- 2.1 Covariables derivadas de los datos ---
municipio_df <- municipio_df %>%
  mutate(
    # Densidad de fuentes (proxy de desarrollo industrial)
    densidad_fuentes = n_fuentes,
    # Intensidad de carbono (CO2 por unidad de actividad)
    intensidad_carbono = ifelse(actividad_total > 0, 
                                co2_total / actividad_total, 
                                NA),
    # Log de emisiones (normalización)
    log_co2 = log1p(co2_total)
  )

# --- 2.2 Covariables geográficas simuladas (basadas en patrones reales de Colombia) ---
# NOTA METODOLÓGICA: Estas variables se simulan con patrones geográficos documentados
# para fines del ejercicio metodológico. En un estudio real, se obtendrían de:
# - DANE (PIB municipal)
# - IDEAM (radiación, viento)
# - IGAC (cobertura de suelo)

set.seed(2026)
n <- nrow(municipio_df)

municipio_df <- municipio_df %>%
  mutate(
    # PIB per cápita simulado: correlacionado con latitud (centro = más rico)
    # Patrón real: Bogotá, Medellín, Cali tienen mayor PIB
    pib_percapita = 8000 + 
      4000 * exp(-((lat_med - 4.5)^2 + (lon_med + 74)^2) / 50) +
      rnorm(n, 0, 1500),
    
    # Radiación solar: mayor en norte (Caribe) y zonas bajas
    # Patrón real: La Guajira > 5.5 kWh/m²/día, sur < 4.0
    radiacion_solar = 5.5 - 0.15 * lat_med + rnorm(n, 0, 0.3),
    
    # Velocidad del viento: mayor en zonas costeras y cordilleras
    viento_ms = 3.5 + 0.8 * abs(lat_med - 8) + 
      0.3 * (lon_med > -75) + rnorm(n, 0, 0.5),
    
    # Población (proxy: densidad de fuentes * factor)
    poblacion = round(exp(log(n_fuentes + 1) * 2.5 + rnorm(n, 8, 0.5))),
    
    # Distancia a red eléctrica (proxy: inversa a densidad urbana)
    # Zonas remotas (Amazonía, Chocó) tienen menor acceso
    distancia_red_km = pmax(0, 50 - 10 * log1p(n_fuentes) + 
                              ifelse(lon_med < -75 & lat_med < 2, 30, 0) +
                              rnorm(n, 0, 10)),
    
    # Cobertura de suelo: % de bosque (mayor en sur y pacífico)
    cobertura_bosque_pct = pmin(100, pmax(0,
      60 - 5 * abs(lat_med - 2) - 3 * (lon_med > -73) + rnorm(n, 0, 10)
    ))
  )

# Verificar que no haya NAs en covariables
covariables <- c("pib_percapita", "radiacion_solar", "viento_ms", 
                 "poblacion", "distancia_red_km", "cobertura_bosque_pct",
                 "densidad_fuentes", "intensidad_carbono")

for (v in covariables) {
  n_na <- sum(is.na(municipio_df[[v]]))
  if (n_na > 0) {
    cat(sprintf("  ⚠ %s tiene %d NAs - imputando con mediana\n", v, n_na))
    municipio_df[[v]][is.na(municipio_df[[v]])] <- median(municipio_df[[v]], na.rm = TRUE)
  }
}

cat("  • Covariables construidas:\n")
for (v in covariables) {
  cat(sprintf("    - %s: [%.2f, %.2f]\n", v, 
              min(municipio_df[[v]]), max(municipio_df[[v]])))
}

# ==============================================================================
# 3. DIAGRAMA ACÍCLICO DIRIGIDO (DAG) CAUSAL ----------------------------------
# ==============================================================================
cat("\n=== 3. Construyendo DAG causal ===\n")

# DAG basado en la literatura de energía y emisiones
# Nodos:
#   T = Tratamiento (tipo de fuente energética)
#   Y = Resultado (emisiones CO2 per cápita)
#   X = Covariables de confusión (PIB, radiación, viento, población, distancia, cobertura)
#   U = Confusores no observados (cultura ambiental, políticas locales)

dag_code <- "
dag {
  # Tratamiento
  T [exposure]
  
  # Resultado
  Y [outcome]
  
  # Covariables observadas (confusores)
  PIB
  radiacion
  viento
  poblacion
  distancia_red
  cobertura_suelo
  
  # Confusor no observado
  U [unobserved]
  
  # Flechas causales
  T -> Y
  PIB -> T
  PIB -> Y
  radiacion -> T
  radiacion -> Y
  viento -> T
  viento -> Y
  poblacion -> T
  poblacion -> Y
  distancia_red -> T
  distancia_red -> Y
  cobertura_suelo -> T
  cobertura_suelo -> Y
  
  # Confusor no observado afecta tratamiento y resultado
  U -> T
  U -> Y
}"

dag <- dagitty(dag_code)

# Verificar ajustabilidad (backdoor criterion)
cat("\n  • Conjuntos de ajuste válidos (criterio de backdoor):\n")
ajustes <- adjustmentSets(dag, exposure = "T", outcome = "Y")
print(ajustes)

# Identificar variables instrumentales (si las hay)
cat("\n  • Variables instrumentales potenciales:\n")
iv <- instrumentalVariables(dag, exposure = "T", outcome = "Y")
print(iv)

# Visualizar DAG
png("results/01_dag_causal.png", width = 800, height = 600, res = 150)
plot(dag, cex = 1.2, col = c("treatment" = "red", "outcome" = "blue", 
                              "adjusted" = "green", "unadjusted" = "gray"))
title("DAG Causal: Tratamiento Energético → Emisiones CO2", 
      cex.main = 1.3, font.main = 2)
dev.off()
cat("  ✓ DAG guardado en results/01_dag_causal.png\n")

# ==============================================================================
# 4. ESTIMACIÓN DE PROPENSITY SCORES ------------------------------------------
# ==============================================================================
cat("\n=== 4. Estimando Propensity Scores ===\n")

# Preparar datos para MatchIt
# Tratamiento como factor multinomial (5 niveles)
municipio_df$tratamiento <- factor(municipio_df$tratamiento, 
                                    levels = c("T0_fosil", "T1_solar", 
                                               "T2_eolica", "T3_hidro", "T4_biomasa"))

# Fórmula de propensity score
# Usamos todas las covariables observadas del DAG
formula_ps <- tratamiento ~ pib_percapita + radiacion_solar + viento_ms + 
              poblacion + distancia_red_km + cobertura_bosque_pct + 
              densidad_fuentes + intensidad_carbono

cat("  • Fórmula de propensity score:\n")
cat("    ", deparse(formula_ps), "\n")

# Verificar tamaño de cada grupo
cat("\n  • Tamaño de grupos antes del emparejamiento:\n")
print(table(municipio_df$tratamiento))

# ==============================================================================
# 5. EMPAREJAMIENTO CON MATCHIT -----------------------------------------------
# ==============================================================================
cat("\n=== 5. Ejecutando emparejamiento por Propensity Score ===\n")

# Estrategia: Emparejamiento 1:1 con caliper para cada tratamiento vs control (T0)
# T0 (fósil) es el grupo de referencia/control

# 5.1 Emparejamiento multinomial (todos los tratamientos simultáneamente)
cat("\n  • Emparejamiento multinomial (todos vs todos)...\n")

match_multi <- matchit(
  formula_ps,
  data = municipio_df,
  method = "nearest",           # Vecino más cercano
  distance = "glm",             # Propensity score vía GLM multinomial
  ratio = 1,                    # 1:1 matching
  caliper = 0.2,                # Caliper de 0.2 desviaciones estándar del logit
  replace = FALSE,              # Sin reemplazo
  discard = "control",          # Descartar controles fuera de soporte común
  estimand = "ATE"              # Efecto promedio del tratamiento
)

cat(sprintf("  • Observaciones emparejadas: %d / %d (%.1f%%)\n", 
            sum(match_multi$weights > 0), nrow(municipio_df),
            100 * sum(match_multi$weights > 0) / nrow(municipio_df)))

# 5.2 Emparejamientos pairwise (cada tratamiento vs T0 fósil)
# Más interpretable para estimar ATT
tratamientos <- c("T1_solar", "T2_eolica", "T3_hidro", "T4_biomasa")
matches_pairwise <- list()

for (trt in tratamientos) {
  cat(sprintf("\n  • Emparejando %s vs T0_fosil...\n", trt))
  
  # Subconjunto binario
  df_bin <- municipio_df %>%
    filter(tratamiento %in% c("T0_fosil", trt)) %>%
    mutate(trat_bin = ifelse(tratamiento == trt, 1, 0))
  
  if (sum(df_bin$trat_bin == 1) < 5 || sum(df_bin$trat_bin == 0) < 5) {
    cat(sprintf("    ⚠ Grupo demasiado pequeño, saltando\n"))
    next
  }
  
  m <- matchit(
    trat_bin ~ pib_percapita + radiacion_solar + viento_ms + 
      poblacion + distancia_red_km + cobertura_bosque_pct + 
      densidad_fuentes + intensidad_carbono,
    data = df_bin,
    method = "nearest",
    distance = "glm",
    ratio = 1,
    caliper = 0.2,
    replace = FALSE,
    estimand = "ATT"  # Efecto en los tratados
  )
  
  matches_pairwise[[trt]] <- m
  cat(sprintf("    ✓ Pares formados: %d\n", sum(m$match.matrix != 0)))
}

# ==============================================================================
# 6. VALIDACIÓN DE BALANCE - LOVE PLOTS ---------------------------------------
# ==============================================================================
cat("\n=== 6. Validando balance pre/post-emparejamiento ===\n")

# 6.1 Love plot del emparejamiento multinomial
love_multi <- love.plot(
  match_multi,
  threshold = 0.1,            # Línea de balance aceptable (|SMD| < 0.1)
  abs = TRUE,
  stats = c("mean.diff"),     # Diferencia de medias estandarizada
  var.order = "unadjusted",   # Ordenar por desbalance pre-matching
  title = "Balance de Covariables - Emparejamiento Multinomial"
)

ggsave("results/02_love_plot_multinomial.png", love_multi, 
       width = 10, height = 8, dpi = 300)
cat("  ✓ Love plot multinomial guardado\n")

# 6.2 Love plots pairwise
for (trt in names(matches_pairwise)) {
  love_pw <- love.plot(
    matches_pairwise[[trt]],
    threshold = 0.1,
    abs = TRUE,
    stats = c("mean.diff"),
    title = sprintf("Balance: %s vs T0_fósil", trt)
  )
  
  fname <- gsub(" ", "_", tolower(trt))
  ggsave(sprintf("results/03_love_%s.png", fname), love_pw, 
         width = 10, height = 6, dpi = 300)
  cat(sprintf("  ✓ Love plot %s guardado\n", trt))
}

# 6.3 Tabla de balance numérica
cat("\n  • Tabla de balance (diferencias de medias estandarizadas):\n")
balance_tab <- summary(match_multi, un = TRUE)$bal.tab$Balance

# Mostrar solo covariables
cat("\n  Covariables con |SMD| > 0.1 post-emparejamiento (desbalance residual):\n")
post_balance <- balance_tab[grepl("e.adj", rownames(balance_tab)), ]
desbalanceadas <- rownames(post_balance)[abs(post_balance$e.adj.Mean.Diff) > 0.1]

if (length(desbalanceadas) > 0) {
  cat("    ", paste(desbalanceadas, collapse = ", "), "\n")
} else {
  cat("    ✓ Todas las covariables balanceadas (|SMD| < 0.1)\n")
}

# Guardar tabla completa
write.csv(balance_tab, "results/04_balance_table.csv")

# ==============================================================================
# 7. ESTIMACIÓN DE EFECTOS CAUSALES (ATE, ATT, CATE) --------------------------
# ==============================================================================
cat("\n=== 7. Estimando efectos causales ===\n")

# 7.1 ATE (Average Treatment Effect) - Emparejamiento multinomial
# Usar pesos de MatchIt para estimar ATE
matched_data <- match.data(match_multi)

# Modelo de resultado ponderado
# Y = log_co2 (intensidad de emisiones)
model_ate <- lm(log_co2 ~ tratamiento + pib_percapita + radiacion_solar + 
                  viento_ms + poblacion + distancia_red_km + 
                  cobertura_bosque_pct + densidad_fuentes + intensidad_carbono,
                data = matched_data,
                weights = weights)

cat("\n  • Estimaciones ATE (vs T0_fósil como referencia):\n")
ate_results <- summary(model_ate)$coefficients
print(round(ate_results[grep("tratamiento", rownames(ate_results)), ], 4))

# 7.2 ATT (Average Treatment Effect on the Treated) - Pairwise
cat("\n  • Estimaciones ATT (emparejamiento pairwise):\n")
att_results <- list()

for (trt in names(matches_pairwise)) {
  md <- match.data(matches_pairwise[[trt]])
  
  # Diferencia simple de medias en muestra emparejada
  y_trat <- mean(md$log_co2[md$trat_bin == 1])
  y_cont <- mean(md$log_co2[md$trat_bin == 0])
  att <- y_trat - y_cont
  
  # Error estándar (diferencia de medias pareadas)
  se_att <- sd(md$log_co2[md$trat_bin == 1] - 
                 md$log_co2[md$trat_bin == 0]) / sqrt(sum(md$trat_bin == 1))
  
  # Test t pareado
  t_test <- t.test(log_co2 ~ trat_bin, data = md, paired = TRUE)
  
  att_results[[trt]] <- data.frame(
    tratamiento = trt,
    ATT = att,
    SE = se_att,
    t_stat = t_test$statistic,
    p_value = t_test$p.value,
    IC_95_inf = t_test$conf.int[1],
    IC_95_sup = t_test$conf.int[2],
    n_tratados = sum(md$trat_bin == 1),
    n_controles = sum(md$trat_bin == 0)
  )
  
  cat(sprintf("    %s: ATT = %.4f (p = %.4f) [IC95: %.4f, %.4f]\n",
              trt, att, t_test$p.value, 
              t_test$conf.int[1], t_test$conf.int[2]))
}

att_df <- do.call(rbind, att_results)
write.csv(att_df, "results/05_att_estimates.csv", row.names = FALSE)

# 7.3 CATE (Conditional ATE) - Heterogeneidad por covariable
# Estratificar por PIB (alto/bajo) para ver heterogeneidad
cat("\n  • CATE estratificado por nivel de PIB:\n")

matched_data <- matched_data %>%
  mutate(pib_alto = pib_percapita > median(pib_percapita))

cate_model <- lm(log_co2 ~ tratamiento * pib_alto + pib_percapita + 
                   radiacion_solar + viento_ms + poblacion,
                 data = matched_data,
                 weights = weights)

cat("\n  • Interacción tratamiento × PIB (heterogeneidad):\n")
interacciones <- summary(cate_model)$coefficients
print(round(interacciones[grep("tratamiento.*pib", rownames(interacciones)), ], 4))

# ==============================================================================
# 8. ANÁLISIS DE SENSIBILIDAD DE ROSENBAUM ------------------------------------
# ==============================================================================
cat("\n=== 8. Análisis de sensibilidad de Rosenbaum ===\n")

# El análisis de sensibilidad evalúa cuán fuerte debe ser un confusor no observado
# (U) para invalidar nuestras conclusiones.

# Usamos el emparejamiento T1_solar vs T0_fosil como ejemplo
if ("T1_solar" %in% names(matches_pairwise)) {
  md_solar <- match.data(matches_pairwise[["T1_solar"]])
  
  # Crear variable de resultado binaria para sensitivity2x2xk
  # (alta emisión vs baja emisión)
  md_solar <- md_solar %>%
    mutate(emision_alta = as.numeric(log_co2 > median(log_co2)))
  
  # Tabla 2x2xk (tratamiento × resultado × estrato de emparejamiento)
  tabla_2x2xk <- xtabs(~ trat_bin + emision_alta + subclass, 
                       data = md_solar)
  
  # Bounds de Rosenbaum
  cat("\n  • Bounds de Rosenbaum para T1_solar vs T0_fósil:\n")
  cat("    (Gamma = fuerza del confusor no observado)\n\n")
  
  gammas <- seq(1, 3, by = 0.25)
  sens_results <- data.frame(
    Gamma = gammas,
    p_upper = NA_real_,
    p_lower = NA_real_
  )
  
  for (i in seq_along(gammas)) {
    g <- gammas[i]
    # Test de Wilcoxon signed-rank con bounds
    ps <- psbounds(tabla_2x2xk, Gamma = g)
    sens_results$p_upper[i] <- ps$upper
    sens_results$p_lower[i] <- ps$lower
  }
  
  print(sens_results)
  
  # Interpretación
  gamma_critico <- sens_results$Gamma[sens_results$p_upper < 0.05]
  if (length(gamma_critico) > 0) {
    cat(sprintf("\n  • Γ crítico = %.2f: el resultado es robusto a confusores\n", 
                max(gamma_critico)))
    cat("    no observados con odds ratio hasta este valor.\n")
  }
  
  # Gráfico de sensibilidad
  p_sens <- ggplot(sens_results, aes(x = Gamma)) +
    geom_line(aes(y = p_upper), color = "red", linewidth = 1.2, 
              label = "p-value upper bound") +
    geom_line(aes(y = p_lower), color = "blue", linewidth = 1.2,
              label = "p-value lower bound") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
    annotate("text", x = 2.5, y = 0.045, label = "α = 0.05", 
             color = "gray40", size = 3) +
    labs(
      title = "Análisis de Sensibilidad de Rosenbaum",
      subtitle = "T1_solar vs T0_fósil - ¿Cuán robusto es el ATT?",
      x = "Gamma (fuerza del confusor no observado)",
      y = "p-value bounds",
      caption = "Si Γ=1 (sin confusión oculta), recuperamos el test estándar"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("results/06_rosenbaum_sensitivity.png", p_sens, 
         width = 8, height = 6, dpi = 300)
  cat("  ✓ Gráfico de sensibilidad guardado\n")
  
  write.csv(sens_results, "results/07_rosenbaum_bounds.csv", row.names = FALSE)
}

# ==============================================================================
# 9. RESUMEN Y ENTREGABLES ----------------------------------------------------
# ==============================================================================
cat("\n=== 9. Resumen de entregables ===\n")

# Tabla resumen para el paper
resumen_final <- data.frame(
  Metrica = c(
    "N municipios analizados",
    "N tratamientos",
    "Covariables de ajuste",
    "Caliper utilizado",
    "Método de emparejamiento",
    "Estimando principal",
    "Balance logrado (|SMD| < 0.1)",
    "Robustez Rosenbaum (Γ crítico)"
  ),
  Valor = c(
    nrow(municipio_df),
    length(levels(municipio_df$tratamiento)),
    length(covariables),
    "0.2 DE del logit",
    "Nearest neighbor (1:1)",
    "ATE (multinomial) / ATT (pairwise)",
    ifelse(length(desbalanceadas) == 0, "Sí - todas", "Parcial"),
    ifelse(exists("gamma_critico"), 
           sprintf("%.2f", max(gamma_critico)), 
           "No calculado")
  )
)

write.csv(resumen_final, "results/08_resumen_causal.csv", row.names = FALSE)

cat("\n  Archivos generados:\n")
cat("  ✓ results/01_dag_causal.png              - DAG causal\n")
cat("  ✓ results/02_love_plot_multinomial.png   - Balance multinomial\n")
cat("  ✓ results/03_love_*.png                  - Balance pairwise\n")
cat("  ✓ results/04_balance_table.csv           - Tabla de balance\n")
cat("  ✓ results/05_att_estimates.csv           - Estimaciones ATT\n")
cat("  ✓ results/06_rosenbaum_sensitivity.png   - Sensibilidad\n")
cat("  ✓ results/07_rosenbaum_bounds.csv        - Bounds numéricos\n")
cat("  ✓ results/08_resumen_causal.csv          - Resumen\n")

# Guardar datos emparejados para el script LMM
saveRDS(matched_data, "data/matched_data.rds")
cat("  ✓ data/matched_data.rds                  - Para script 04b_lmm.R\n")

cat("\n=== Análisis causal completado ===\n")
cat("\n⚠ NOTA METODOLÓGICA PARA EL PAPER:\n")
cat("  Las covariables PIB, radiación, viento, distancia a red y cobertura\n")
cat("  de suelo fueron simuladas con patrones geográficos documentados de\n")
cat("  Colombia. En un estudio aplicado real, estas deben obtenerse de fuentes\n")
cat("  oficiales (DANE, IDEAM, IGAC) para garantizar validez externa.\n")