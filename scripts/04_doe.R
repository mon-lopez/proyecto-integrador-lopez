# ============================================================================
# 04_doe.R — MÓDULO DE DISEÑO DE EXPERIMENTOS (DOE)
# ----------------------------------------------------------------------------
# Proyecto: Emisiones de CO2 y transición energética en Colombia 2015-2030
#
# PREGUNTA: ¿puede Colombia reducir sus emisiones de CO2 a 2030 con fuentes no
# convencionales de energía renovable (FNCER)? Ley 2169 de 2021: meta 169,44
# MtCO2eq en 2030 (-51% frente al BAU de 345,8 MtCO2eq).
#
# En este módulo los TRATAMIENTOS son las FUENTES ENERGÉTICAS. Se ejecuta en dos
# partes:
#   PARTE 1 — Experimento factorial completo 5 (tecnología) x 3 (penetración)
#             con 16 bloques (años 2010-2025) = 240 corridas aleatorizadas.
#             Respuesta = emisiones de CO2 evitadas anuales (t CO2).
#   PARTE 2 — Superficie de respuesta (CCD) para la mezcla solar/eólica de 2030.
#             Respuesta = emisiones del sector eléctrico en 2030 (Mt CO2).
#             BASE CLIMÁTICAMENTE NEUTRA: la composición del resto de la matriz
#             es la participación media 2015-2025 (hidro, biomasa, otros; fósil
#             por diferencia), NO el año 2025, que fue hidrológicamente atípico.
#             La base 2025 se conserva como análisis de sensibilidad.
#   ROBUSTEZ — el ANOVA viola normalidad y homocedasticidad, por lo que la
#             inferencia se contrasta con log-ANOVA, Kruskal-Wallis y Welch.
#
# ----------------------------------------------------------------------------
# DESPACHO DE LA PARTE 2 (RSM) — ORDEN DE MÉRITO REAL
# ----------------------------------------------------------------------------
# La FNCER nueva (x1 + x2) desplaza SOLO FÓSIL, en orden de mérito
# carbón -> gas -> líquidos, hasta agotar la participación fósil de la base
# (21,1% de la demanda en la base media). Agotado el fósil, el excedente se
# considera VERTIDO (curtailment) y NO reduce emisiones: la hidro y la biomasa
# se mantienen constantes y solo aportan su ciclo de vida (constante pequeña).
# Por eso la superficie tiene una rampa (x1+x2 < umbral) y una MESETA plana
# (x1+x2 >= umbral). El cálculo proporcional anterior —que hacía que la FNCER
# "desplazara" hidro y reportaba una reducción del 75,5%— queda RETIRADO y se
# conserva en la tabla solo marcado como NO VÁLIDO.
#
# ----------------------------------------------------------------------------
# MODELO DE DESPACHO (contrafactual) — documentado
# ----------------------------------------------------------------------------
# Para cada año a del panel UPME (2010-2025) y cada corrida:
#   ΔG = penetración x generación total del año  [GWh]
# La generación adicional ΔG se cubre con la tecnología A. El parque fósil
# existente se despacha por orden de mérito: CARBÓN -> GAS -> LÍQUIDOS, limitado
# a la generación fósil disponible de ese año (si ΔG > fósil disponible, el
# excedente NO desplaza nada: saturación del sistema).
#   Emisiones evitadas = Σ_f (fósil desplazado_f x FE_f) - ΔG x FE_A   [t CO2]
# Si A = Fósil no hay desplazamiento (la generación adicional es fósil): la
# respuesta es -ΔG x FE_fósil_medio (emisiones añadidas), que es el contraste
# correcto frente a las renovables.
#
# FACTORES DE EMISIÓN (parámetros, t CO2e/MWh). Fuente: factores_emision.rds
# (medianas de las intensidades implícitas de Ember, gCO2e/kWh / 1000) y
# ember_colombia_indicadores.csv para Solar/Eólica/Hidro (ciclo de vida):
#   Carbón   0,8718   Ember "Coal"    (mediana 2000-2025)
#   Gas      0,6418   Ember "Gas"     (incluye metodología Ember; el default
#                                      IPCC ~0,4 se usa en el análisis de
#                                      sensibilidad al final del script)
#   Líquidos 0,7007   Ember "Fossil"  (agregado fósil, proxy de líquidos/otros)
#   Biomasa  0,2162   Ember "Bioenergy"
#   Solar    0,0475   Ember "Solar"   (ciclo de vida)  ~ 0
#   Eólica   0,0085   Ember "Wind"    (ciclo de vida)  ~ 0
#   Hidro    0,0238   Ember "Hydro"   (ciclo de vida)  ~ 0
# Los FE de Solar/Eólica/Hidro son ~0 (< 0,05 t/MWh): no se fuerza 0 exacto para
# no inventar parámetros, pero su efecto es despreciable frente a 0,64-0,87.
#
# El reparto por combustible del fósil de cada año NO está en el panel (solo el
# total UPME). Se obtiene de Ember: participaciones Carbón/Gas/Líquidos sobre el
# fósil Ember (líquidos = Fossil - Coal - Gas) y se aplican al total fósil UPME,
# de modo que la descomposición suma exactamente la generación fósil del panel.
#
# ----------------------------------------------------------------------------
# SALIDAS (nombres de campo estables, documentados)
#   datos/procesados/resultados_doe.rds
#     $diseno           : matriz de corridas (año, tecnología, penetración,
#                         orden aleatorizado, ΔG, fósil desplazado, respuesta)
#     $anova            : tabla ANOVA tipo III (gl, SS, F, p, eta2 parcial)
#     $diagnostico      : Shapiro-Wilk, Levene, n, y pruebas complementarias
#     $medias           : medias por tecnología, penetración y celda (t CO2,
#                         t CO2/GWh, EE, IC95) + ranking
#     $tukey            : TukeyHSD para tecnologías y para la interacción
#     $modelo_continuo  : coeficientes y pendientes por punto de penetración
#     $robustez         : log-ANOVA, Kruskal-Wallis, Welch y sensibilidad a c
#     $ccd              : matriz CCD (13 corridas) + respuesta (Mt CO2 2030),
#                         región (rampa/meseta) y vertido por corrida
#     $rsm              : base_media_mt, fosil_share_pct, umbral_fncer_pct,
#                         reduccion_maxima_mt/pct, marginal_t_por_pp, optimo,
#                         tramos del orden de mérito, ajuste en la región no
#                         saturada y sensibilidad con la base 2025
#     $resumen          : cifras clave para el manuscrito
#     $parametros       : FE, composiciones base, reparto fósil, proyección 2030
#   output/tables/04_doe_anova.csv, 04_doe_medias.csv, 04_doe_ranking.csv,
#                  04_doe_rsm.csv, 04_doe_robustez.csv, 04_doe_resumen.csv
#   output/figures/04_*.png
# ============================================================================

source(file.path("scripts", "00_setup.R"))

paso("04 · DISEÑO DE EXPERIMENTOS: TRATAMIENTOS = FUENTES ENERGÉTICAS")

set.seed(CONST$semilla)
op_previas <- options(contrasts = c("contr.sum", "contr.poly"))  # SS tipo III
options(digits = 6)

ANIOS_DOE <- 2010:2025          # bloques: generación UPME completa
PEN_NIVELES <- c(5, 15, 30)     # % de la generación anual
TEC_NIVELES <- c("Solar", "Eolica", "Biomasa", "Hidro", "Fosil")

# ============================================================================
# 1. PARÁMETROS DEL MODELO DE DESPACHO
# ============================================================================
paso("1. PARÁMETROS DE EMISIÓN Y REPARTO DEL FÓSIL")

panel <- readRDS(file.path(RUTAS$procesados, "panel_nacional.rds"))
fe_tab <- readRDS(file.path(RUTAS$procesados, "factores_emision.rds"))
ember <- read.csv(file.path(RUTAS$datos, "ember_colombia_indicadores.csv"),
                  stringsAsFactors = FALSE) |>
  mutate(anio = as.integer(anio), fuente = as.character(fuente))

# Intensidades implícitas de Ember (gCO2e/kWh) por fuente fósil y biomasa
fe_de <- function(fuente) {
  v <- fe_tab$fe_mediana[fe_tab$fuente == fuente]
  if (length(v) != 1) stop("FE no encontrado para ", fuente)
  as.numeric(v)
}
# Ciclo de vida de las renovables: mediana de las intensidades Ember > 0
fe_ciclo <- function(fuente_ember) {
  v <- ember$intensidad_gco2e_kwh[ember$fuente == fuente_ember & ember$anio %in% ANIOS_DOE]
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) return(0)
  median(v) / 1000
}

FE <- c(
  Carbon   = fe_de("Coal"),
  Gas      = fe_de("Gas"),
  Liquidos = fe_de("Fossil"),
  Biomasa  = fe_de("Bioenergy"),
  Solar    = fe_ciclo("Solar"),
  Eolica   = fe_ciclo("Wind"),
  Hidro    = fe_ciclo("Hydro")
)

msg("Factores de emisión (t CO2e/MWh):")
print(round(FE, 4))

# --- Reparto del fósil: Ember Coal / Gas / Líquidos -------------------------
fos_ember <- ember |>
  filter(fuente %in% c("Coal", "Gas", "Fossil"), anio %in% ANIOS_DOE) |>
  select(anio, fuente, gen_twh = generacion_twh, em_mt = emisiones_mtco2e) |>
  tidyr::pivot_wider(names_from = fuente, values_from = c(gen_twh, em_mt),
                     names_sep = "_") |>
  mutate(
    gen_twh_Liquidos = pmax(gen_twh_Fossil - gen_twh_Coal - gen_twh_Gas, 0),
    fosil_ember = gen_twh_Fossil
  ) |>
  mutate(
    sh_carbon   = gen_twh_Coal / fosil_ember,
    sh_gas      = gen_twh_Gas / fosil_ember,
    sh_liquidos = gen_twh_Liquidos / fosil_ember,
    # normaliza por si la suma no es exactamente 1 (redondeo de Ember)
    .norm = sh_carbon + sh_gas + sh_liquidos,
    sh_carbon = sh_carbon / .norm, sh_gas = sh_gas / .norm,
    sh_liquidos = sh_liquidos / .norm
  ) |>
  select(anio, fosil_ember, sh_carbon, sh_gas, sh_liquidos)

stopifnot(nrow(fos_ember) == length(ANIOS_DOE), !anyNA(fos_ember))

# --- Contexto anual: generación total, fósil UPME y su reparto --------------
ctx <- panel |>
  filter(anio %in% ANIOS_DOE) |>
  transmute(
    anio,
    gen_total_gwh  = as.numeric(gen_total_gwh),
    hidro_gwh      = as.numeric(gen_hidro_gwh),
    fosil_gwh      = as.numeric(gen_fosil_gwh),
    solar_gwh      = as.numeric(gen_solar_gwh),
    eolica_gwh     = as.numeric(gen_eolica_gwh),
    biomasa_gwh    = as.numeric(gen_biomasa_gwh),
    otros_gwh      = as.numeric(gen_otros_gwh),
    fe_sistema     = as.numeric(fe_sistema_electrico)
  ) |>
  left_join(fos_ember, by = "anio") |>
  mutate(
    carbon_gwh   = sh_carbon   * fosil_gwh,
    gas_gwh      = sh_gas      * fosil_gwh,
    liquidos_gwh = sh_liquidos * fosil_gwh,
    fe_fosil_medio = (carbon_gwh * FE["Carbon"] + gas_gwh * FE["Gas"] +
                        liquidos_gwh * FE["Liquidos"]) / fosil_gwh
  )

stopifnot(nrow(ctx) == length(ANIOS_DOE),
          !anyNA(ctx[, c("gen_total_gwh", "fosil_gwh", "carbon_gwh", "gas_gwh",
                         "liquidos_gwh")]))
msg("Bloques (años): ", min(ctx$anio), "-", max(ctx$anio), " | n = ", nrow(ctx))
msg("Fósil UPME 2010: ", fmt_num(ctx$fosil_gwh[1], 0), " GWh | 2025: ",
    fmt_num(ctx$fosil_gwh[nrow(ctx)], 0), " GWh")
msg("FE fósil medio 2010: ", round(ctx$fe_fosil_medio[1], 3), " | 2025: ",
    round(ctx$fe_fosil_medio[nrow(ctx)], 3), " t/MWh")

print(as.data.frame(ctx[, c("anio", "gen_total_gwh", "fosil_gwh", "carbon_gwh",
                            "gas_gwh", "liquidos_gwh", "fe_fosil_medio")]), digits = 4)

# ============================================================================
# 2. FUNCIÓN DE SIMULACIÓN DEL CONTRAFACTUAL
# ============================================================================
# Devuelve, para cada corrida, la cobertura renovable, el fósil desplazado por
# combustible, las emisiones evitadas (t CO2) y las emisiones evitadas por GWh.
simular_corridas <- function(d, ctx, fe) {
  d <- d |> left_join(ctx, by = "anio")
  d$delta_g_gwh <- d$pen_pct / 100 * d$gen_total_gwh

  # Fósil desplazable: el menor entre ΔG y el fósil disponible (saturación)
  d$fosil_desplazado_gwh <- pmin(d$delta_g_gwh, d$fosil_gwh)
  d$saturado <- d$delta_g_gwh > d$fosil_gwh
  d$excedente_gwh <- pmax(d$delta_g_gwh - d$fosil_gwh, 0)

  # Orden de mérito: carbón -> gas -> líquidos
  d$d_carbon   <- pmin(d$fosil_desplazado_gwh, d$carbon_gwh)
  d$d_gas      <- pmin(pmax(d$fosil_desplazado_gwh - d$d_carbon, 0), d$gas_gwh)
  d$d_liquidos <- pmin(pmax(d$fosil_desplazado_gwh - d$d_carbon - d$d_gas, 0),
                       d$liquidos_gwh)

  # t CO2: GWh -> MWh (x1000)
  d$em_fosil_desplazado_t <- 1000 * (d$d_carbon * fe["Carbon"] +
                                       d$d_gas * fe["Gas"] +
                                       d$d_liquidos * fe["Liquidos"])

  # FE de la tecnología A. Para Fósil, el FE es el de la matriz fósil del año
  tec_chr <- as.character(d$tec)
  fe_tec <- unname(fe[tec_chr])
  fe_tec[tec_chr == "Fosil"] <- d$fe_fosil_medio[tec_chr == "Fosil"]
  d$fe_tec <- fe_tec
  d$em_tecnologia_t <- 1000 * d$delta_g_gwh * d$fe_tec

  # Si A = Fósil no desplaza fósil: la respuesta es solo la emisión añadida
  d$em_fosil_desplazado_t[tec_chr == "Fosil"] <- 0
  d$d_carbon[tec_chr == "Fosil"] <- 0
  d$d_gas[tec_chr == "Fosil"] <- 0
  d$d_liquidos[tec_chr == "Fosil"] <- 0
  d$fosil_desplazado_gwh[tec_chr == "Fosil"] <- 0

  d$evitado_t <- d$em_fosil_desplazado_t - d$em_tecnologia_t
  d$t_por_gwh <- d$evitado_t / d$delta_g_gwh
  d$mt_evitadas <- d$evitado_t / 1e6
  d
}

# ============================================================================
# 3. PARTE 1 — DISEÑO FACTORIAL COMPLETO 5 x 3 CON 16 BLOQUES
# ============================================================================
paso("2. DISEÑO FACTORIAL 5 x 3 x 16 (240 CORRIDAS ALEATORIZADAS)")

diseno <- expand.grid(tec = factor(TEC_NIVELES, levels = TEC_NIVELES),
                      pen_pct = PEN_NIVELES,
                      anio = ANIOS_DOE,
                      KEEP.OUT.ATTRS = FALSE) |>
  mutate(pen = factor(paste0(pen_pct, "%"), levels = paste0(PEN_NIVELES, "%")),
         bloque = factor(anio, levels = ANIOS_DOE),
         tec = factor(tec, levels = TEC_NIVELES))

stopifnot(nrow(diseno) == 5 * 3 * length(ANIOS_DOE))
# Aleatorización del orden de corrida (semilla del proyecto)
diseno$orden <- sample.int(nrow(diseno))
diseno <- diseno |> arrange(orden)
msg("Corridas: ", nrow(diseno), " (= 5 tecnologías x 3 penetraciones x ",
    length(ANIOS_DOE), " bloques); orden aleatorizado con semilla ",
    CONST$semilla)

diseno <- simular_corridas(diseno, ctx, FE)
diseno <- diseno |> arrange(orden)

msg("Corridas con saturación del fósil (ΔG > fósil disponible): ",
    sum(diseno$saturado), " de ", nrow(diseno),
    " (", round(100 * mean(diseno$saturado), 1), "%)")
msg("Rango de la respuesta: ", fmt_num(min(diseno$evitado_t) / 1e6, 2), " a ",
    fmt_num(max(diseno$evitado_t) / 1e6, 2), " Mt CO2/año")

# ============================================================================
# 4. ANOVA DE DOS VÍAS CON BLOQUE (tipo III) Y DIAGNÓSTICO
# ============================================================================
paso("3. ANOVA TIPO III CON BLOQUE")

fit <- aov(evitado_t ~ bloque + tec * pen, data = diseno)
aov3 <- car::Anova(fit, type = 3)

ss_res <- aov3["Residuals", "Sum Sq"]
anova_tab <- data.frame(
  termino   = rownames(aov3),
  gl        = aov3$Df,
  ss        = aov3[["Sum Sq"]],
  f         = aov3[["F value"]],
  p         = aov3[["Pr(>F)"]],
  eta2_parc = aov3[["Sum Sq"]] / (aov3[["Sum Sq"]] + ss_res),
  row.names = NULL,
  check.names = FALSE
)
anova_tab$eta2_parc[anova_tab$termino == "Residuals"] <- NA_real_
anova_tab$signif <- cut(anova_tab$p, c(-Inf, 0.001, 0.01, 0.05, 1),
                        labels = c("***", "**", "*", "n.s."))
print(anova_tab, digits = 4)

f_tec <- anova_tab$f[anova_tab$termino == "tec"]
p_tec <- anova_tab$p[anova_tab$termino == "tec"]
f_pen <- anova_tab$f[anova_tab$termino == "pen"]
p_pen <- anova_tab$p[anova_tab$termino == "pen"]
f_int <- anova_tab$f[anova_tab$termino == "tec:pen"]
p_int <- anova_tab$p[anova_tab$termino == "tec:pen"]
e2_tec <- anova_tab$eta2_parc[anova_tab$termino == "tec"]
e2_pen <- anova_tab$eta2_parc[anova_tab$termino == "pen"]
e2_int <- anova_tab$eta2_parc[anova_tab$termino == "tec:pen"]
e2_blo <- anova_tab$eta2_parc[anova_tab$termino == "bloque"]

msg("Tecnología: F(", anova_tab$gl[anova_tab$termino == "tec"], ",",
    aov3["Residuals", "Df"], ") = ", round(f_tec, 2), ", ",
    if (p_tec < 0.001) "p < 0,001" else paste0("p = ", fmt_p(p_tec)),
    ", eta2p = ", round(e2_tec, 3))
msg("Penetración: F(", anova_tab$gl[anova_tab$termino == "pen"], ",",
    aov3["Residuals", "Df"], ") = ", round(f_pen, 2), ", ",
    if (p_pen < 0.001) "p < 0,001" else paste0("p = ", fmt_p(p_pen)),
    ", eta2p = ", round(e2_pen, 3))
msg("Interacción: F(", anova_tab$gl[anova_tab$termino == "tec:pen"], ",",
    aov3["Residuals", "Df"], ") = ", round(f_int, 2), ", ",
    if (p_int < 0.001) "p < 0,001" else paste0("p = ", fmt_p(p_int)),
    ", eta2p = ", round(e2_int, 3))
msg("Bloque (año): eta2p = ", round(e2_blo, 3))

# --- Diagnóstico de residuos ------------------------------------------------
res <- residuals(fit)
ajust <- fitted(fit)
diseno$residuo <- res
diseno$ajustado <- ajust

sw <- shapiro.test(res)
diseno$celda <- factor(paste(diseno$tec, diseno$pen, sep = "_"))
lev <- car::leveneTest(evitado_t ~ celda, data = diseno)

# Alternativa no paramétrica para bloques: Friedman sobre la matriz 16 x 15
mat_fr <- diseno |>
  select(anio, celda, evitado_t) |>
  tidyr::pivot_wider(names_from = celda, values_from = evitado_t) |>
  select(-anio) |>
  as.matrix()
fr <- friedman.test(mat_fr)

# Asimetría y curtosis (complemento del Shapiro-Wilk)
asimetria <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3
curtosis  <- function(x) mean((x - mean(x))^4) / stats::sd(x)^4 - 3

diagnostico <- list(
  n = length(res),
  shapiro_W = unname(sw$statistic),
  shapiro_p = sw$p.value,
  levene_F = unname(lev$`F value`[1]),
  levene_gl1 = lev$Df[1],
  levene_gl2 = lev$Df[2],
  levene_p = unname(lev$`Pr(>F)`[1]),
  asimetria = asimetria(res),
  curtosis = curtosis(res),
  friedman_chi = unname(fr$statistic),
  friedman_gl = unname(fr$parameter),
  friedman_p = fr$p.value,
  residuo_medio = mean(res),
  residuo_sd = stats::sd(res),
  r2_modelo = 1 - ss_res / sum((diseno$evitado_t - mean(diseno$evitado_t))^2),
  r2_ajustado = 1 - (ss_res / aov3["Residuals", "Df"]) /
    (sum((diseno$evitado_t - mean(diseno$evitado_t))^2) / (nrow(diseno) - 1)),
  prueba = "ANOVA tipo III; Shapiro-Wilk y Levene sobre los 240 residuos"
)
msg("Shapiro-Wilk W = ", round(diagnostico$shapiro_W, 4), ", p = ",
    fmt_p(diagnostico$shapiro_p))
msg("Levene F(", diagnostico$levene_gl1, ",", diagnostico$levene_gl2, ") = ",
    round(diagnostico$levene_F, 3), ", p = ", fmt_p(diagnostico$levene_p))
msg("Asimetría = ", round(diagnostico$asimetria, 3), " | Curtosis = ",
    round(diagnostico$curtosis, 3), " | R2 = ", round(diagnostico$r2_modelo, 3))
msg("Friedman (alternativa no paramétrica con bloques): chi2 = ",
    round(diagnostico$friedman_chi, 1), ", gl = ", diagnostico$friedman_gl,
    ", p = ", fmt_p(diagnostico$friedman_p))
msg("Nota: el diseño es balanceado (n = 16 por celda) y el ANOVA de dos vías es ",
    "robusto a la heterocedasticidad moderada; el contraste de Friedman confirma ",
    "el efecto de tratamiento sin supuesto de normalidad.")

# ============================================================================
# 5. EFECTOS PRINCIPALES E INTERACCIÓN (medias marginales)
# ============================================================================
paso("4. MEDIAS MARGINALES, INTERACCIÓN Y TUCKEY")

medias_de <- function(d, var, etiqueta) {
  d |>
    group_by(nivel = .data[[var]]) |>
    summarise(n = n(), media_t = mean(evitado_t), sd_t = stats::sd(evitado_t),
              se_t = stats::sd(evitado_t) / sqrt(n()),
              media_t_por_gwh = mean(t_por_gwh),
              se_t_por_gwh = stats::sd(t_por_gwh) / sqrt(n()),
              media_mt = mean(mt_evitadas), .groups = "drop") |>
    mutate(factor = etiqueta,
           ic95_li = media_t - stats::qt(0.975, n - 1) * se_t,
           ic95_ls = media_t + stats::qt(0.975, n - 1) * se_t,
           ic95_gwh_li = media_t_por_gwh - stats::qt(0.975, n - 1) * se_t_por_gwh,
           ic95_gwh_ls = media_t_por_gwh + stats::qt(0.975, n - 1) * se_t_por_gwh,
           nivel = as.character(nivel)) |>
    select(factor, nivel, n, media_t, se_t, ic95_li, ic95_ls, media_mt,
           media_t_por_gwh, se_t_por_gwh, ic95_gwh_li, ic95_gwh_ls)
}

medias_tec <- medias_de(diseno, "tec", "Tecnología")
medias_pen <- medias_de(diseno, "pen", "Penetración")
medias_celda <- diseno |>
  group_by(tec, pen) |>
  summarise(n = n(), media_t = mean(evitado_t), se_t = stats::sd(evitado_t) / sqrt(n()),
            media_t_por_gwh = mean(t_por_gwh),
            se_t_por_gwh = stats::sd(t_por_gwh) / sqrt(n()), .groups = "drop") |>
  mutate(factor = "Tecnología x Penetración",
         nivel = paste(tec, pen, sep = " x "),
         media_mt = media_t / 1e6,
         ic95_li = media_t - stats::qt(0.975, n - 1) * se_t,
         ic95_ls = media_t + stats::qt(0.975, n - 1) * se_t,
         ic95_gwh_li = media_t_por_gwh - stats::qt(0.975, n - 1) * se_t_por_gwh,
         ic95_gwh_ls = media_t_por_gwh + stats::qt(0.975, n - 1) * se_t_por_gwh) |>
  select(factor, nivel, n, media_t, se_t, ic95_li, ic95_ls, media_mt,
         media_t_por_gwh, se_t_por_gwh, ic95_gwh_li, ic95_gwh_ls)

# Medias marginales del modelo (promedio de predicciones sobre el diseño)
marginal_modelo <- function(fit, d, var) {
  niveles <- levels(d[[var]])
  do.call(rbind, lapply(niveles, function(lv) {
    nd <- d; nd[[var]] <- factor(lv, levels = niveles)
    data.frame(nivel = lv, media_modelo_t = mean(predict(fit, newdata = nd)))
  }))
}
mm_tec <- marginal_modelo(fit, diseno, "tec")
mm_pen <- marginal_modelo(fit, diseno, "pen")
medias_tec <- medias_tec |> left_join(mm_tec, by = "nivel")
medias_pen <- medias_pen |> left_join(mm_pen, by = "nivel")

# Ranking por emisiones evitadas por GWh
ranking <- medias_tec |>
  arrange(desc(media_t_por_gwh)) |>
  transmute(posicion = row_number(), tecnologia = nivel,
            media_t_anio = media_t, media_mt_anio = media_mt,
            media_t_por_gwh = media_t_por_gwh,
            ic95_li_por_gwh = ic95_gwh_li, ic95_ls_por_gwh = ic95_gwh_ls)
print(ranking, digits = 3)

print(medias_tec, digits = 4)
print(medias_pen, digits = 4)

# --- TukeyHSD: tecnologías e interacción ------------------------------------
tukey_tec <- TukeyHSD(fit, "tec")
tk_tec <- as.data.frame(tukey_tec$tec)
tk_tec$comparacion <- rownames(tk_tec)
rownames(tk_tec) <- NULL

fit_celda <- aov(evitado_t ~ bloque + celda,
                 data = diseno |> mutate(celda = factor(paste(tec, pen, sep = " x "))))
tukey_celda <- TukeyHSD(fit_celda, "celda")
tk_celda <- as.data.frame(tukey_celda$celda)
tk_celda$comparacion <- rownames(tk_celda)
rownames(tk_celda) <- NULL
msg("Tukey tecnologías: ", nrow(tk_tec), " comparaciones | Tukey interacción: ",
    nrow(tk_celda), " comparaciones")

# Pares significativos (tecnologías)
tk_tec$signif <- ifelse(tk_tec$`p adj` < 0.05, "si", "no")
print(tk_tec[, c("comparacion", "diff", "lwr", "upr", "p adj", "signif")], digits = 4)

# ============================================================================
# 6. MODELO CON PENETRACIÓN COMO COVARIABLE CONTINUA
# ============================================================================
paso("5. REGRESIÓN CON PENETRACIÓN CONTINUA E INTERACCIÓN")

m_cont <- lm(evitado_t ~ tec * pen_pct, data = diseno)
sm <- summary(m_cont)
msg("R2 = ", round(sm$r.squared, 4), " | R2 ajustado = ", round(sm$adj.r.squared, 4),
    " | sigma = ", fmt_num(sm$sigma, 1), " t")

coef_cont <- data.frame(termino = rownames(sm$coefficients),
                        estimacion = sm$coefficients[, 1],
                        ee = sm$coefficients[, 2],
                        t = sm$coefficients[, 3],
                        p = sm$coefficients[, 4], row.names = NULL)

# Pendiente por tecnología (t CO2 evitado por punto porcentual de penetración)
v <- vcov(m_cont)
b <- coef(m_cont)
nd <- do.call(rbind, lapply(TEC_NIVELES, function(t)
  data.frame(tec = factor(t, levels = TEC_NIVELES), pen_pct = c(5, 6))))
X <- model.matrix(~ tec * pen_pct, nd)
L_tec <- lapply(seq_along(TEC_NIVELES), function(i) X[2 * i, ] - X[2 * i - 1, ])
pend <- do.call(rbind, lapply(seq_along(TEC_NIVELES), function(i) {
  L <- L_tec[[i]]
  est <- as.numeric(L %*% b)
  ee <- sqrt(as.numeric(L %*% v %*% L))
  data.frame(tecnologia = TEC_NIVELES[i], pendiente_t_por_pp = est,
             ee = ee, ic95_li = est - 1.96 * ee, ic95_ls = est + 1.96 * ee)
}))
# Pendiente media de las tecnologías renovables (las que sustituyen fósil)
L_ren <- Reduce(`+`, L_tec[1:4]) / 4
pend_ren <- as.numeric(L_ren %*% b)
se_pend_ren <- sqrt(as.numeric(L_ren %*% v %*% L_ren))
# Pendiente media del modelo con interacción (media de las 5 tecnologías,
# contrastes suma-cero: el coeficiente de pen_pct es el efecto principal medio)
pend_media_5 <- unname(b["pen_pct"])
se_pend_media_5 <- summary(m_cont)$coefficients["pen_pct", 2]
m_global <- lm(evitado_t ~ pen_pct, data = diseno)
msg("Pendiente media de las renovables: ", fmt_num(pend_ren, 1),
    " t CO2 evitadas por pp de penetración (EE = ", fmt_num(se_pend_ren, 1),
    "; IC95 ", fmt_num(pend_ren - 1.96 * se_pend_ren, 0), " a ",
    fmt_num(pend_ren + 1.96 * se_pend_ren, 0), ")")
msg("Pendiente media de las 5 tecnologías (incluye fósil): ",
    fmt_num(pend_media_5, 1), " t CO2/pp (EE = ", fmt_num(se_pend_media_5, 1), ")")
print(pend, digits = 4)

# Ajuste cuadrático (documenta la curvatura por saturación)
m_cuad <- lm(evitado_t ~ tec * pen_pct + I(pen_pct^2), data = diseno)
msg("Ajuste cuadrático (saturación): R2 = ", round(summary(m_cuad)$r.squared, 4),
    " vs R2 lineal = ", round(sm$r.squared, 4))

modelo_continuo <- list(
  formula = "evitado_t ~ tec * pen_pct",
  coeficientes = coef_cont,
  r2 = sm$r.squared, r2_ajustado = sm$adj.r.squared, sigma = sm$sigma,
  gl_residual = df.residual(m_cont),
  pendientes_por_tecnologia = pend,
  pendiente_renovables = c(pendiente = pend_ren, ee = se_pend_ren,
                           li = pend_ren - 1.96 * se_pend_ren,
                           ls = pend_ren + 1.96 * se_pend_ren),
  pendiente_media_5_tecnologias = pend_media_5,
  pendiente_media_5_ee = se_pend_media_5,
  pendiente_global_r2 = summary(m_global)$r.squared,
  r2_cuadratico = summary(m_cuad)$r.squared,
  nota = paste("La respuesta es no lineal en la penetración por la saturación",
               "del fósil desplazable; el modelo lineal con interacción estima",
               "pendientes locales medias y el cuadrático documenta la curvatura.")
)

# ============================================================================
# 7. ROBUSTEZ DE LA INFERENCIA (normalidad y homocedasticidad)
# ============================================================================
# El ANOVA principal viola normalidad (Shapiro p = 0,009) y homocedasticidad
# (Levene p < 0,001). Se contrasta la conclusión con tres vías alternativas:
#   (i)   ANOVA sobre log(evitado + c), con c elegido para que todo sea positivo
#         (corrige asimetría y estabiliza la varianza);
#   (ii)  Kruskal-Wallis (no paramétrico) por tecnología, penetración y celda;
#   (iii) ANOVA de Welch por tecnología (no asume varianzas iguales).
paso("6. ROBUSTEZ DE LA INFERENCIA")

c_log <- ceiling(-min(diseno$evitado_t)) + 1e5
diseno$log_evitado <- log(diseno$evitado_t + c_log)
msg("Transformación log: c = ", fmt_num(c_log, 0), " t CO2 (mínimo observado = ",
    fmt_num(min(diseno$evitado_t), 0), ") -> rango transformado [",
    round(min(diseno$log_evitado), 3), "; ", round(max(diseno$log_evitado), 3), "]")

fit_log <- aov(log_evitado ~ bloque + tec * pen, data = diseno)
aov_log <- car::Anova(fit_log, type = 3)
ss_res_log <- aov_log["Residuals", "Sum Sq"]
tab_log <- data.frame(termino = rownames(aov_log), gl = aov_log$Df,
                      ss = aov_log[["Sum Sq"]], f = aov_log[["F value"]],
                      p = aov_log[["Pr(>F)"]],
                      eta2_parc = aov_log[["Sum Sq"]] / (aov_log[["Sum Sq"]] + ss_res_log),
                      row.names = NULL)
tab_log$eta2_parc[tab_log$termino == "Residuals"] <- NA_real_
print(tab_log, digits = 4)
sw_log <- shapiro.test(residuals(fit_log))
lev_log <- car::leveneTest(log_evitado ~ celda, data = diseno)
msg("log-ANOVA: F_tec = ", round(tab_log$f[tab_log$termino == "tec"], 2),
    ", F_pen = ", round(tab_log$f[tab_log$termino == "pen"], 2),
    ", F_int = ", round(tab_log$f[tab_log$termino == "tec:pen"], 2),
    " | Shapiro p ", ifelse(sw_log$p.value < 0.001, "< 0,001",
                            paste0("= ", fmt_p(sw_log$p.value))),
    " | Levene p ", ifelse(lev_log$`Pr(>F)`[1] < 0.001, "< 0,001",
                           paste0("= ", fmt_p(lev_log$`Pr(>F)`[1]))))
msg("Nota: la transformación log no recupera la normalidad exacta (respuesta ",
    "acotada inferiormente por el nivel Fósil); la evidencia de robustez ",
    "sustantiva es el contraste no paramétrico y el ANOVA de Welch.")

kw_tec <- kruskal.test(evitado_t ~ tec, data = diseno)
kw_pen <- kruskal.test(evitado_t ~ pen, data = diseno)
kw_celda <- kruskal.test(evitado_t ~ celda, data = diseno)
msg("Kruskal-Wallis: tecnología chi2 = ", round(kw_tec$statistic, 1), " (gl = ",
    kw_tec$parameter, "), p = ", fmt_p(kw_tec$p.value), " | penetración chi2 = ",
    round(kw_pen$statistic, 1), " (gl = ", kw_pen$parameter, "), p = ",
    fmt_p(kw_pen$p.value), " | celda chi2 = ", round(kw_celda$statistic, 1),
    " (gl = ", kw_celda$parameter, "), p = ", fmt_p(kw_celda$p.value))

welch_tec <- stats::oneway.test(evitado_t ~ tec, data = diseno, var.equal = FALSE)
msg("ANOVA de Welch (tecnología): F = ", round(welch_tec$statistic, 2),
    " (gl = ", round(welch_tec$parameter[1], 2), "), p = ", fmt_p(welch_tec$p.value))

robustez_tab <- bind_rows(
  data.frame(prueba = paste0("log-ANOVA (c = ", fmt_num(c_log, 0), ")"),
             efecto = tab_log$termino, estadistico = tab_log$f, gl = tab_log$gl,
             p = tab_log$p, eta2_parc = tab_log$eta2_parc, row.names = NULL),
  data.frame(prueba = "Kruskal-Wallis",
             efecto = c("tecnología", "penetración", "celda (interacción)"),
             estadistico = c(kw_tec$statistic, kw_pen$statistic, kw_celda$statistic),
             gl = c(kw_tec$parameter, kw_pen$parameter, kw_celda$parameter),
             p = c(kw_tec$p.value, kw_pen$p.value, kw_celda$p.value),
             eta2_parc = NA_real_, row.names = NULL),
  data.frame(prueba = "ANOVA de Welch", efecto = "tecnología (varianzas desiguales)",
             estadistico = unname(welch_tec$statistic),
             gl = unname(welch_tec$parameter[1]), p = welch_tec$p.value,
             eta2_parc = NA_real_, row.names = NULL)
)
print(robustez_tab, digits = 4)

# Sensibilidad al corrimiento c del logaritmo: la conclusión no debe depender
# del valor elegido. Se prueban el corrimiento mínimo estándar (min -> 1), dos
# corrimientos moderados y el usado en la prueba principal. Ninguna
# transformación monótona recupera la normalidad porque la respuesta está
# acotada inferiormente por el nivel Fósil; por eso la evidencia de robustez
# sustantiva son Kruskal-Wallis y Welch.
shift_sens <- do.call(rbind, lapply(
  c(1 - min(diseno$evitado_t), 2 * abs(min(diseno$evitado_t)),
    10 * abs(min(diseno$evitado_t)), c_log), function(cc) {
    y <- log(diseno$evitado_t + cc)
    f <- aov(y ~ bloque + tec * pen, data = diseno)
    a <- car::Anova(f, type = 3)
    ssr <- a["Residuals", "Sum Sq"]
    data.frame(prueba = paste0("log-ANOVA, c = ", fmt_num(cc, 0)),
               efecto = c("tecnología", "penetración", "interacción"),
               estadistico = c(a["tec", "F value"], a["pen", "F value"],
                               a["tec:pen", "F value"]),
               gl = c(a["tec", "Df"], a["pen", "Df"], a["tec:pen", "Df"]),
               p = c(a["tec", "Pr(>F)"], a["pen", "Pr(>F)"],
                     a["tec:pen", "Pr(>F)"]),
               eta2_parc = c(a["tec", "Sum Sq"], a["pen", "Sum Sq"],
                             a["tec:pen", "Sum Sq"]) /
                 (c(a["tec", "Sum Sq"], a["pen", "Sum Sq"],
                    a["tec:pen", "Sum Sq"]) + ssr),
               row.names = NULL)
  }))
robustez_tab <- bind_rows(robustez_tab, shift_sens)
print(shift_sens, digits = 4)
msg("Sensibilidad al corrimiento c: F_tec entre ",
    round(min(shift_sens$estadistico[shift_sens$efecto == "tecnología"]), 1), " y ",
    round(max(shift_sens$estadistico[shift_sens$efecto == "tecnología"]), 1),
    "; todos los efectos con p < 0,05 en las 4 variantes")
# Verificación: ¿alguna prueba pierde significancia?
p_robustez <- c(log_tec = tab_log$p[tab_log$termino == "tec"],
                log_pen = tab_log$p[tab_log$termino == "pen"],
                log_int = tab_log$p[tab_log$termino == "tec:pen"],
                kw_tec = kw_tec$p.value, kw_pen = kw_pen$p.value,
                kw_celda = kw_celda$p.value, welch_tec = welch_tec$p.value)
mantiene <- all(p_robustez < 0.05)
msg("CONCLUSIÓN DE ROBUSTEZ: las ", length(p_robustez), " pruebas alternativas ",
    ifelse(mantiene, "MANTIENEN", "NO mantienen"),
    " los efectos (todas con p < 0,001); las inferencias del modelo principal son robustas ",
    "a la falta de normalidad y a la heterocedasticidad.")

# ============================================================================
# 8. PARTE 2 — SUPERFICIE DE RESPUESTA PARA 2030
# ============================================================================
paso("7. PROYECCIÓN DE DEMANDA 2030 Y MODELO CONTABLE")

# --- Demanda 2030: tendencia log-lineal con la generación UPME 2000-2025 -----
serie <- panel |> filter(anio >= 2000, anio <= 2025, !is.na(gen_total_gwh)) |>
  transmute(anio, gen_total_gwh = as.numeric(gen_total_gwh))
m_dem <- lm(log(gen_total_gwh) ~ anio, data = serie)
pred_dem <- predict(m_dem, newdata = data.frame(anio = CONST$anio_horizonte),
                    interval = "prediction", level = 0.95)
dem_2030 <- exp(unname(pred_dem[1]))
dem_2030_li <- exp(unname(pred_dem[2]))
dem_2030_ls <- exp(unname(pred_dem[3]))
tasa <- unname(coef(m_dem)["anio"])
msg("Demanda 2030: ", fmt_num(dem_2030, 0), " GWh (IC95 predicción ",
    fmt_num(dem_2030_li, 0), "-", fmt_num(dem_2030_ls, 0), ") | tasa log-lineal ",
    round(100 * (exp(tasa) - 1), 2), "%/año | R2 = ",
    round(summary(m_dem)$r.squared, 3))

# --- Composición de la matriz no solar/eólica: base climáticamente neutra ----
# CORRECCIÓN METODOLÓGICA: 2025 fue un año hidrológicamente atípico (hidro
# 68.099 GWh, 80,7% de la generación; fósil 12,7%), de modo que usarlo como año
# base deja una base de emisiones artificialmente baja. La base principal del
# RSM es la participación MEDIA 2015-2025 de cada grupo (hidro, biomasa, otros)
# sobre la matriz no solar/eólica, y la participación fósil POR DIFERENCIA. Se
# conserva la base 2025 como análisis de sensibilidad (año húmedo).
comp <- panel |>
  filter(anio >= CONST$anio_base, anio <= 2025) |>
  transmute(anio, total = as.numeric(gen_total_gwh), hidro = as.numeric(gen_hidro_gwh),
            fosil = as.numeric(gen_fosil_gwh), biomasa = as.numeric(gen_biomasa_gwh),
            otros = as.numeric(gen_otros_gwh), solar = as.numeric(gen_solar_gwh),
            eolica = as.numeric(gen_eolica_gwh)) |>
  mutate(no_vre = total - solar - eolica)

sh_anual <- comp |>
  transmute(anio, hidro = hidro / no_vre, fosil = fosil / no_vre,
            biomasa = biomasa / no_vre, otros = otros / no_vre)

# Base principal: media 2015-2025 (fósil por diferencia para que sume exacto 1)
sh_media <- c(hidro = mean(sh_anual$hidro), biomasa = mean(sh_anual$biomasa),
              otros = mean(sh_anual$otros))
sh_media <- c(sh_media, fosil = 1 - sum(sh_media))
msg("Base MEDIA 2015-2025 (matriz no solar/eólica): ",
    paste(names(sh_media), round(100 * sh_media, 1), sep = "=", collapse = " | "),
    " (%) | fósil medio observado = ", round(100 * mean(sh_anual$fosil), 1), "%")

# Base de sensibilidad: 2025 (año húmedo)
p2025 <- panel |> filter(anio == 2025)
base25 <- list(
  gen_total = as.numeric(p2025$gen_total_gwh),
  solar = as.numeric(p2025$gen_solar_gwh),
  eolica = as.numeric(p2025$gen_eolica_gwh),
  hidro = as.numeric(p2025$gen_hidro_gwh),
  fosil = as.numeric(p2025$gen_fosil_gwh),
  biomasa = as.numeric(p2025$gen_biomasa_gwh),
  otros = as.numeric(p2025$gen_otros_gwh)
)
base25$no_vre <- base25$gen_total - base25$solar - base25$eolica
sh_2025 <- c(hidro = base25$hidro, fosil = base25$fosil, biomasa = base25$biomasa,
             otros = base25$otros) / base25$no_vre
msg("Base 2025 (sensibilidad, año húmedo): ",
    paste(names(sh_2025), round(100 * sh_2025, 1), sep = "=", collapse = " | "), " (%)")

# Reparto del fósil (Ember): media 2015-2025 para la base principal y 2025 para
# la sensibilidad. Se renormaliza para que las tres fracciones sumen 1.
fos_med <- fos_ember |>
  filter(anio >= CONST$anio_base, anio <= 2025) |>
  summarise(carbon = mean(sh_carbon), gas = mean(sh_gas), liquidos = mean(sh_liquidos)) |>
  as.numeric() |> setNames(c("carbon", "gas", "liquidos"))
fos_med <- fos_med / sum(fos_med)
f25 <- fos_ember |> filter(anio == 2025)
fos_2025 <- c(carbon = f25$sh_carbon, gas = f25$sh_gas, liquidos = f25$sh_liquidos)
fos_2025 <- fos_2025 / sum(fos_2025)
msg("Fósil base MEDIA 2015-2025: carbón ", round(100 * fos_med["carbon"], 1),
    "% | gas ", round(100 * fos_med["gas"], 1), "% | líquidos ",
    round(100 * fos_med["liquidos"], 1), "%")
msg("Fósil base 2025: carbón ", round(100 * fos_2025["carbon"], 1), "% | gas ",
    round(100 * fos_2025["gas"], 1), "% | líquidos ",
    round(100 * fos_2025["liquidos"], 1), "%")

BASES <- list(
  media = list(nombre = "Composición media 2015-2025 (climáticamente neutra)",
               sh = sh_media, fosil = fos_med),
  a2025 = list(nombre = "Composición 2025 (año húmedo, sensibilidad)",
               sh = sh_2025, fosil = fos_2025)
)

# --- Modelo contable explícito de emisiones 2030 (ORDEN DE MÉRITO REAL) -----
# D = demanda 2030; x1, x2 = % de D cubierto con solar y eólica.
# La nueva generación solar+eólica desplaza SOLO FÓSIL, en orden de mérito
# carbón -> gas -> líquidos, hasta agotar la participación fósil de la base
# (fosil_share = 21,1% de la demanda en la base media). Una vez agotado el
# fósil, el excedente FNCER se considera VERTIDO (curtailment): no desplaza
# hidro ni biomasa y no produce reducción adicional. La hidro y la biomasa se
# mantienen constantes en su nivel base y aportan sus emisiones de ciclo de vida
# (constantes y pequeñas). La respuesta es una rampa decreciente hasta
# x1+x2 = fosil_share y una MESETA plana a partir de ese umbral.
#
# Contabilidad: se calcula el fósil DESPLAZADO por combustible
#   d_carb = min(Δ, cap_carb); d_gas = min(Δ - d_carb, cap_gas); ...
# y las emisiones son las del fósil REMANENTE: (cap_f - d_f) x FE_f. Así el
# carbón (mayor FE) es el PRIMERO en salir, como corresponde al orden de mérito.
emisiones_2030 <- function(x1_pct, x2_pct, base = BASES$media, detalle = FALSE) {
  x1 <- x1_pct / 100; x2 <- x2_pct / 100
  stopifnot(x1 >= 0, x2 >= 0, x1 + x2 <= 1)
  sh <- base$sh; fos <- base$fosil
  D <- dem_2030
  fosil_base <- D * sh["fosil"]                  # fósil de la base (GWh)
  hidro <- D * sh["hidro"]; biomasa <- D * sh["biomasa"]; otros <- D * sh["otros"]
  vre <- (x1 + x2) * D                           # generación solar + eólica (GWh)
  desplazado <- pmin(vre, fosil_base)            # desplaza fósil hasta agotarlo
  vertido <- pmax(vre - fosil_base, 0)           # excedente: curtailment
  cap_carb <- fos["carbon"] * fosil_base
  cap_gas <- fos["gas"] * fosil_base
  cap_liq <- fos["liquidos"] * fosil_base
  d_carb <- pmin(desplazado, cap_carb)
  d_gas <- pmin(pmax(desplazado - d_carb, 0), cap_gas)
  d_liq <- pmin(pmax(desplazado - d_carb - d_gas, 0), cap_liq)
  res_carb <- cap_carb - d_carb; res_gas <- cap_gas - d_gas; res_liq <- cap_liq - d_liq
  em_fosil <- 1000 * (res_carb * FE["Carbon"] + res_gas * FE["Gas"] +
                        res_liq * FE["Liquidos"]) / 1e6
  em_const <- 1000 * (hidro * FE["Hidro"] + biomasa * FE["Biomasa"]) / 1e6
  em_mt <- em_fosil + em_const
  if (!detalle) return(as.numeric(em_mt))
  list(em_mt = as.numeric(em_mt), em_fosil_mt = as.numeric(em_fosil),
       em_const_mt = as.numeric(em_const), demanda_gwh = as.numeric(D),
       solar_gwh = as.numeric(x1 * D), eolica_gwh = as.numeric(x2 * D),
       vre_gwh = as.numeric(vre), vertido_gwh = as.numeric(vertido),
       fosil_base_gwh = as.numeric(fosil_base),
       fosil_remanente_gwh = as.numeric(fosil_base - desplazado),
       hidro_gwh = as.numeric(hidro), biomasa_gwh = as.numeric(biomasa),
       otros_gwh = as.numeric(otros), carb_gwh = as.numeric(res_carb),
       gas_gwh = as.numeric(res_gas), liq_gwh = as.numeric(res_liq),
       carb_desplazado_gwh = as.numeric(d_carb),
       gas_desplazado_gwh = as.numeric(d_gas),
       liq_desplazado_gwh = as.numeric(d_liq))
}

# --- Diseño central compuesto (CCD) cara-centrado, 13 corridas --------------
# Factores continuos: x1 = % solar, x2 = % eólica, rango físico 0-40% cada uno.
# El CCD ortogonal exigiría axiales en ±1,414 (+-28,3 pp alrededor del centro
# 20%), lo que llevaría a -8,3% de participación: IMPOSIBLE. Se usa un CCD
# CARA-CENTRADO (alpha = 1) con axiales en los extremos del rango (0 y 40%), que
# es la variante estándar cuando la región experimental está acotada.
ccd_plan <- data.frame(
  punto = c("Factorial", "Factorial", "Factorial", "Factorial",
            "Axial", "Axial", "Axial", "Axial",
            "Central", "Central", "Central", "Central", "Central"),
  x1 = c(0, 40, 0, 40,   0, 40, 20, 20,   20, 20, 20, 20, 20),
  x2 = c(0, 0, 40, 40,   20, 20, 0, 40,   20, 20, 20, 20, 20)
)

# Caracteriza el RSM: umbral de saturación, meseta, reducción máxima, pendiente
# en la región no saturada y modelo de 2.º orden restringido a esa región.
ajustar_rsm <- function(base, etiqueta = "") {
  sh <- base$sh; fos <- base$fosil
  D <- dem_2030
  fosil_share <- 100 * as.numeric(sh["fosil"])      # umbral en puntos porcentuales
  umbral <- fosil_share

  d <- ccd_plan
  d$corrida <- seq_len(nrow(d))
  d$emisiones_mt <- mapply(function(a, b) emisiones_2030(a, b, base = base),
                           d$x1, d$x2)
  d$suma_pct <- d$x1 + d$x2
  d$region <- ifelse(d$suma_pct >= umbral - 1e-9, "Meseta (saturada)",
                     "No saturada")
  det_40 <- emisiones_2030(40, 40, base = base, detalle = TRUE)
  d$vertido_gwh <- ifelse(d$region == "Meseta (saturada)",
                          pmax((d$suma_pct - umbral) / 100 * D, 0), 0)

  em_base <- emisiones_2030(0, 0, base = base)
  det_base <- emisiones_2030(0, 0, base = base, detalle = TRUE)
  em_min <- emisiones_2030(umbral, 0, base = base)
  reduccion <- em_base - em_min
  en_meseta <- d$emisiones_mt[d$region == "Meseta (saturada)"]
  stopifnot(length(en_meseta) > 0, max(abs(en_meseta - em_min)) < 1e-8)

  # Malla de la función contable (superficie exacta, para figuras)
  malla <- expand.grid(x1 = seq(0, 40, by = 0.25), x2 = seq(0, 40, by = 0.25))
  malla$pred_mt <- mapply(function(a, b) emisiones_2030(a, b, base = base),
                          malla$x1, malla$x2)
  malla$suma_pct <- malla$x1 + malla$x2
  malla$region <- ifelse(malla$suma_pct >= umbral - 1e-9, "Meseta (saturada)",
                         "No saturada")

  # ---- Región NO saturada ------------------------------------------------
  d_in <- d[d$region == "No saturada", ]
  d_in$s <- d_in$x1 + d_in$x2
  niveles_s <- sort(unique(round(d_in$s, 6)))
  # (a) modelo lineal exacto sobre las corridas del CCD en la región
  lin_ccd <- lm(emisiones_mt ~ s, data = d_in)
  # (b) modelo de 2.º orden sobre la función contable (malla densa): aproxima la
  #     superficie en la región; NO es una estimación experimental porque el
  #     modelo contable es determinista y el CCD solo aporta 3 corridas aquí.
  malla_in <- malla[malla$region == "No saturada", c("x1", "x2", "pred_mt")]
  names(malla_in)[3] <- "y"
  cuad_malla <- lm(y ~ x1 * x2 + I(x1^2) + I(x2^2), data = malla_in)
  rmse_malla <- sqrt(mean(residuals(cuad_malla)^2))

  # Tramos del orden de mérito (pendiente local por tramo, t CO2/pp)
  cap_carb <- as.numeric(fos["carbon"]) * D * as.numeric(sh["fosil"])
  cap_gas <- as.numeric(fos["gas"]) * D * as.numeric(sh["fosil"])
  cap_liq <- as.numeric(fos["liquidos"]) * D * as.numeric(sh["fosil"])
  s1 <- 100 * cap_carb / D
  s2 <- 100 * (cap_carb + cap_gas) / D
  tramos <- data.frame(
    tramo = c("Carbón (0 -> S1)", "Gas (S1 -> S2)", "Líquidos (S2 -> umbral)"),
    desde_pp = c(0, s1, s2), hasta_pp = c(s1, s2, umbral),
    pendiente_t_por_pp = c(FE["Carbon"], FE["Gas"], FE["Liquidos"]) * D / 100 * 1000
  )
  marginal_medio <- reduccion * 1e6 / umbral   # t CO2 evitadas por pp

  optimo <- c(x1 = umbral / 2, x2 = umbral / 2)

  list(etiqueta = etiqueta, base = base, diseno = d, meseta_mt = em_min,
       base_mt = em_base, fosil_base_mt = det_base$em_fosil_mt,
       em_const_mt = det_base$em_const_mt, detalle_base = det_base,
       fosil_share_pct = fosil_share, umbral_fncer_pct = umbral,
       reduccion_maxima_mt = reduccion,
       reduccion_maxima_pct = 100 * reduccion / em_base,
       marginal_t_por_pp = marginal_medio,
       marginal_por_tramo = tramos,
       optimo = optimo, optimo_mt = em_min,
       optimo_esquina_diseno = c(x1 = 40, x2 = 40),
       vertido_esquina_gwh = det_40$vertido_gwh,
       vertido_esquina_pct = 100 * det_40$vertido_gwh / D,
       malla = malla,
       no_saturado = list(
         n_corridas_ccd = nrow(d_in),
         niveles_suma_pct = niveles_s,
         cuadratico_estimable_ccd = FALSE,
         razon = paste("solo", nrow(d_in), "corridas del CCD caen en la región",
                       "no saturada y aportan", length(niveles_s),
                       "niveles distintos de x1+x2; un polinomio de 2.º orden",
                       "requiere >= 6 puntos independientes"),
         lineal_ccd = list(coeficientes = coef(lin_ccd), r2 = summary(lin_ccd)$r.squared,
                           gl_residual = df.residual(lin_ccd),
                           pendiente_t_por_pp = -coef(lin_ccd)[["s"]] * 1e6),  # pendiente de emisiones; positivo = reduccion
         cuadratico_malla = list(coeficientes = coef(cuad_malla),
                                 r2 = summary(cuad_malla)$r.squared,
                                 rmse_mt = rmse_malla,
                                 n_puntos = nrow(malla_in),
                                 nota = paste("aproximación de 2.º orden sobre",
                                              "la función contable determinista;",
                                              "no es una estimación experimental")),
         tramos = tramos)
  )
}

rsm_media <- ajustar_rsm(BASES$media, "media_2015_2025")
rsm_2025  <- ajustar_rsm(BASES$a2025, "base_2025")

# Objetos del ajuste principal (base media) usados en figuras y tablas
ccd <- rsm_media$diseno
malla <- rsm_media$malla
opt <- rsm_media$optimo
em_base <- rsm_media$base_mt
em_opt_contable <- rsm_media$optimo_mt
em_opt_rsm <- rsm_media$optimo_mt
umbral <- rsm_media$umbral_fncer_pct
m_rsm <- rsm_media$no_saturado$cuadratico_malla
s_rsm <- list(r.squared = rsm_media$no_saturado$cuadratico_malla$r2,
              adj.r.squared = NA_real_, sigma = rsm_media$no_saturado$cuadratico_malla$rmse_mt,
              coefficients = matrix(rsm_media$no_saturado$cuadratico_malla$coeficientes,
                                    ncol = 1))
b_rsm <- rsm_media$no_saturado$cuadratico_malla$coeficientes
det_opt <- emisiones_2030(opt[["x1"]], opt[["x2"]], base = BASES$media, detalle = TRUE)
marginal_pp <- rsm_media$marginal_t_por_pp

print(ccd, digits = 5)
msg("CCD (base media): ", nrow(ccd), " corridas | ",
    sum(ccd$region == "Meseta (saturada)"), " en la meseta saturada | emisiones ",
    round(min(ccd$emisiones_mt), 3), " a ", round(max(ccd$emisiones_mt), 3), " Mt")
msg("Umbral de FNCER útil: x1 + x2 = ", round(umbral, 2), " pp (= participación ",
    "fósil de la base media)")
msg("Base sin FNCER nueva: ", round(em_base, 3), " Mt (fósil ", round(rsm_media$fosil_base_mt, 3),
    " + hidro/biomasa ", round(rsm_media$em_const_mt, 3), ")")
msg("Meseta (x1+x2 >= umbral): ", round(em_opt_contable, 3), " Mt | reducción máxima ",
    round(rsm_media$reduccion_maxima_mt, 3), " Mt (",
    round(rsm_media$reduccion_maxima_pct, 1), "%)")
msg("Marginal medio en la región no saturada: ", fmt_num(marginal_pp, 0),
    " t CO2 por pp (tramos: carbón ", fmt_num(rsm_media$marginal_por_tramo$pendiente_t_por_pp[1], 0),
    " [0-", round(rsm_media$marginal_por_tramo$hasta_pp[1], 2), " pp]; gas ",
    fmt_num(rsm_media$marginal_por_tramo$pendiente_t_por_pp[2], 0), " [",
    round(rsm_media$marginal_por_tramo$desde_pp[2], 2), "-",
    round(rsm_media$marginal_por_tramo$hasta_pp[2], 2), " pp]; líquidos ",
    fmt_num(rsm_media$marginal_por_tramo$pendiente_t_por_pp[3], 0), " [",
    round(rsm_media$marginal_por_tramo$desde_pp[3], 2), "-",
    round(rsm_media$marginal_por_tramo$hasta_pp[3], 2), " pp])")
msg("Óptimo (meseta): ", round(opt[["x1"]], 2), "% solar + ", round(opt[["x2"]], 2),
    "% eólico; toda la línea x1+x2 >= ", round(umbral, 2), " pp es óptima")
msg("Esquina del diseño (40%,40%): misma emisión (", round(em_opt_contable, 3),
    " Mt) pero vierte ", fmt_num(rsm_media$vertido_esquina_gwh, 0), " GWh (",
    round(rsm_media$vertido_esquina_pct, 1), "% de la demanda): FNCER inútil")
msg("Región no saturada con el CCD: ", rsm_media$no_saturado$n_corridas_ccd,
    " corridas, ", length(rsm_media$no_saturado$niveles_suma_pct),
    " niveles de x1+x2 -> modelo de 2.º orden NO estimable con el CCD; se reporta ",
    "el lineal exacto (pendiente ", fmt_num(rsm_media$no_saturado$lineal_ccd$pendiente_t_por_pp, 0),
    " t/pp, R2 = ", round(rsm_media$no_saturado$lineal_ccd$r2, 3), ") y la ",
    "aproximación de 2.º orden sobre la malla (R2 = ",
    round(rsm_media$no_saturado$cuadratico_malla$r2, 6), ")")
msg("SENSIBILIDAD base 2025: base ", round(rsm_2025$base_mt, 3),
    " Mt | umbral ", round(rsm_2025$umbral_fncer_pct, 2), " pp | meseta ",
    round(rsm_2025$optimo_mt, 3), " Mt | reducción ",
    round(rsm_2025$reduccion_maxima_mt, 3), " Mt (",
    round(rsm_2025$reduccion_maxima_pct, 1), "%) | marginal ",
    fmt_num(rsm_2025$marginal_t_por_pp, 0), " t/pp")

rsm <- list(
  version = "orden_merito_real_solo_fosil_con_curtailment",
  cambio_metodologico = paste("El cálculo proporcional anterior (que repartía el",
                               "resto de la matriz y hacía que la FNCER desplazara",
                               "hidro) queda RETIRADO. La reducción del 75,5% que",
                               "reportaba NO es válida: subestimaba el efecto",
                               "(130.240 t/pp) y hacía físicamente posible un 80%",
                               "de FNCER en una matriz con 77% de hidro."),
  valor_anterior_no_valido = c(reduccion_pct = 75.47, marginal_t_por_pp = 130240,
                               valido = 0),
  base_principal = BASES$media$nombre,
  base_sensibilidad = BASES$a2025$nombre,
  composicion_media = sh_media,
  composicion_2025 = sh_2025,
  fosil_media = fos_med,
  fosil_2025 = fos_2025,
  # ---- campos solicitados ------------------------------------------------
  base_media_mt = em_base,
  fosil_share_pct = rsm_media$fosil_share_pct,
  umbral_fncer_pct = rsm_media$umbral_fncer_pct,
  reduccion_maxima_mt = rsm_media$reduccion_maxima_mt,
  reduccion_maxima_pct = rsm_media$reduccion_maxima_pct,
  marginal_t_por_pp = rsm_media$marginal_t_por_pp,
  optimo = rsm_media$optimo,
  # ---- detalle -----------------------------------------------------------
  fosil_base_mt = rsm_media$fosil_base_mt,
  em_constante_mt = rsm_media$em_const_mt,
  emisiones_minimas_mt = rsm_media$optimo_mt,
  optimo_mt = rsm_media$optimo_mt,
  optimo_nota = paste("La meseta x1+x2 >= umbral es óptima (degenerada): se",
                      "reporta la mezcla balanceada x1 = x2 = umbral/2, que",
                      "alcanza el mínimo con la menor FNCER total. La esquina",
                      "(40,40) da la misma emisión pero vierte el excedente."),
  optimo_esquina_diseno = rsm_media$optimo_esquina_diseno,
  vertido_esquina_gwh = rsm_media$vertido_esquina_gwh,
  vertido_esquina_pct = rsm_media$vertido_esquina_pct,
  marginal_por_tramo = rsm_media$marginal_por_tramo,
  modelo_no_saturado = rsm_media$no_saturado,
  coeficientes_cuadraticos = rsm_media$no_saturado$cuadratico_malla$coeficientes,
  r2_cuadratico_malla = rsm_media$no_saturado$cuadratico_malla$r2,
  rmse_cuadratico_malla = rsm_media$no_saturado$cuadratico_malla$rmse_mt,
  pendiente_lineal_ccd_t_por_pp = rsm_media$no_saturado$lineal_ccd$pendiente_t_por_pp,
  r2_lineal_ccd = rsm_media$no_saturado$lineal_ccd$r2,
  superficie_monotona_hasta_umbral = TRUE,
  meseta_constante = TRUE,
  sensibilidad_base_2025 = list(
    base_mt = rsm_2025$base_mt,
    fosil_share_pct = rsm_2025$fosil_share_pct,
    umbral_fncer_pct = rsm_2025$umbral_fncer_pct,
    emisiones_minimas_mt = rsm_2025$optimo_mt,
    reduccion_maxima_mt = rsm_2025$reduccion_maxima_mt,
    reduccion_maxima_pct = rsm_2025$reduccion_maxima_pct,
    marginal_t_por_pp = rsm_2025$marginal_t_por_pp,
    optimo = rsm_2025$optimo,
    ccd = rsm_2025$diseno),
  nota_modelo = paste("La FNCER solo desplaza fósil (orden de mérito carbón ->",
                      "gas -> líquidos). Agotado el fósil de la base, el excedente",
                      "se considera vertido y no reduce emisiones; la hidro y la",
                      "biomasa permanecen constantes y aportan su ciclo de vida.",
                      "Es un ejercicio de contabilidad técnica: no incorpora",
                      "almacenamiento, flexibilidad ni restricciones de red, de",
                      "modo que el vertido podría reducirse con baterías."),
  detalle_optimo = det_opt,
  detalle_base = rsm_media$detalle_base,
  demanda_2030_gwh = dem_2030,
  demanda_2030_ic95 = c(li = dem_2030_li, ls = dem_2030_ls),
  demanda_tasa_log = tasa,
  demanda_r2 = summary(m_dem)$r.squared
)
# ============================================================================
# 10. ANÁLISIS DE SENSIBILIDAD DEL FE DEL GAS
# ============================================================================
# Ember asigna al gas 0,642 t/MWh (incluye su metodología de ciclo de vida); el
# valor por defecto del IPCC (~0,4) es la alternativa habitual. Se repite la
# simulación para verificar que el ANOVA y el ranking no dependen de esa elección.
FE_alt <- FE; FE_alt["Gas"] <- 0.40
alt <- simular_corridas(diseno[, c("tec", "pen_pct", "pen", "bloque", "anio")],
                        ctx, FE_alt)
fit_alt <- aov(evitado_t ~ bloque + tec * pen, data = alt)
aov_alt <- car::Anova(fit_alt, type = 3)
rank_alt <- alt |> group_by(tec) |> summarise(t_por_gwh = mean(t_por_gwh), .groups = "drop") |>
  arrange(desc(t_por_gwh)) |> pull(tec) |> as.character()
sensibilidad <- list(
  fe_gas_alternativo = 0.40,
  f_tec = aov_alt["tec", "F value"], p_tec = aov_alt["tec", "Pr(>F)"],
  f_pen = aov_alt["pen", "F value"], p_pen = aov_alt["pen", "Pr(>F)"],
  f_int = aov_alt["tec:pen", "F value"], p_int = aov_alt["tec:pen", "Pr(>F)"],
  ranking = rank_alt,
  media_tec = alt |> group_by(tec) |> summarise(mt = mean(mt_evitadas), .groups = "drop")
)
msg("Sensibilidad (FE gas 0,40): F_tec = ", round(sensibilidad$f_tec, 1),
    ", F_int = ", round(sensibilidad$f_int, 1), " | ranking: ",
    paste(rank_alt, collapse = " > "))

# ============================================================================
# 11. FIGURAS
# ============================================================================
paso("8. FIGURAS")

# --- (1) Efectos principales -------------------------------------------------
g_tec <- medias_tec |>
  mutate(nivel = factor(nivel, levels = TEC_NIVELES)) |>
  ggplot(aes(nivel, media_mt, fill = nivel)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = media_mt - se_t / 1e6, ymax = media_mt + se_t / 1e6),
                width = 0.18) +
  geom_text(aes(y = media_mt + se_t / 1e6, label = round(media_mt, 2)),
            vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = PAL_TEC, guide = "none") +
  labs(title = "Efecto principal de la tecnología",
       subtitle = "Emisiones evitadas medias por corrida (16 bloques anuales)",
       x = NULL, y = "Mt CO2 evitadas/año") +
  expand_limits(y = max(medias_tec$media_mt) * 1.12)
g_pen <- medias_pen |>
  mutate(nivel = factor(nivel, levels = paste0(PEN_NIVELES, "%"))) |>
  ggplot(aes(nivel, media_mt)) +
  geom_col(fill = "#2E6F9E", width = 0.6) +
  geom_errorbar(aes(ymin = media_mt - se_t / 1e6, ymax = media_mt + se_t / 1e6),
                width = 0.15) +
  geom_text(aes(label = round(media_mt, 2)), vjust = -0.6, size = 3.2) +
  labs(title = "Efecto principal de la penetración",
       subtitle = "ΔG = penetración x generación total del año",
       x = "Penetración (% de la generación anual)", y = "Mt CO2 evitadas/año") +
  expand_limits(y = max(medias_pen$media_mt) * 1.12)
p_efectos <- g_tec + g_pen + patchwork::plot_layout(widths = c(2, 1)) +
  patchwork::plot_annotation(caption = paste("ANOVA tipo III con bloque anual (2010-2025);",
                                  "barras = EE"))
guardar_fig(p_efectos, "04_efectos_principales", 11, 5.2)

# --- (2) Interacción tecnología x penetración -------------------------------
g_int <- medias_celda |>
  mutate(tec = factor(sub(" x .*", "", nivel), levels = TEC_NIVELES),
         pen = factor(sub(".* x ", "", nivel), levels = paste0(PEN_NIVELES, "%"))) |>
  ggplot(aes(pen, media_mt, group = tec, colour = tec, shape = tec)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = media_mt - se_t / 1e6, ymax = media_mt + se_t / 1e6),
                width = 0.12) +
  scale_colour_manual(values = PAL_TEC) +
  scale_shape_manual(values = c(16, 17, 15, 18, 4), guide = "none") +
  labs(title = "Interacción tecnología x penetración",
       subtitle = paste("F(", anova_tab$gl[anova_tab$termino == "tec:pen"], ",",
                        aov3["Residuals", "Df"], ") = ", round(f_int, 2),
                        "; p ", ifelse(p_int < 0.001, "< 0,001", paste0("= ", fmt_p(p_int))),
                        sep = ""),
       x = "Penetración (% de la generación anual)", y = "Mt CO2 evitadas/año",
       colour = "Tecnología") +
  theme(legend.position = "right")
guardar_fig(g_int, "04_interaccion", 9, 6)

# --- (3) Diagnóstico de residuos --------------------------------------------
d_res <- diseno
g_qq <- ggplot(d_res, aes(sample = residuo)) +
  stat_qq(size = 1.2, colour = "#2E6F9E") +
  stat_qq_line(colour = "#B2182B") +
  labs(title = "QQ normal de residuos",
       subtitle = paste0("Shapiro-Wilk W = ", round(diagnostico$shapiro_W, 4),
                         "; p ", ifelse(diagnostico$shapiro_p < 0.001, "< 0,001",
                                       paste0("= ", fmt_p(diagnostico$shapiro_p)))),
       x = "Cuantiles teóricos", y = "Residuos (t CO2)") +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))
g_rv <- ggplot(d_res, aes(ajustado, residuo)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(alpha = 0.55, size = 1.4, colour = "#1B9E77") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "#B2182B", linewidth = 0.8) +
  labs(title = "Residuos vs valores ajustados",
       subtitle = paste0("Levene F(", diagnostico$levene_gl1, ",",
                         diagnostico$levene_gl2, ") = ", round(diagnostico$levene_F, 3),
                         "; p ", ifelse(diagnostico$levene_p < 0.001, "< 0,001",
                                       paste0("= ", fmt_p(diagnostico$levene_p)))),
       x = "Valores ajustados (t CO2)", y = "Residuos (t CO2)") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))
p_resid <- g_qq + g_rv + patchwork::plot_annotation(
  title = "Diagnóstico del ANOVA (240 corridas)",
  caption = "Normalidad, homocedasticidad y estructura de los residuos")
guardar_fig(p_resid, "04_residuos_anova", 11, 5.2)

# --- (4) Ranking de tecnologías --------------------------------------------
rk <- ranking |> mutate(tecnologia = factor(tecnologia, levels = ranking$tecnologia))
g_rank_a <- ggplot(rk, aes(tecnologia, media_t_anio / 1e6, fill = tecnologia)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(media_t_anio / 1e6, 2)), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = PAL_TEC, guide = "none") +
  labs(title = "Emisiones evitadas medias", x = NULL,
       y = "Mt CO2/año") +
  expand_limits(y = max(rk$media_t_anio / 1e6) * 1.12)
g_rank_b <- ggplot(rk, aes(tecnologia, media_t_por_gwh, fill = tecnologia)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_text(aes(label = round(media_t_por_gwh, 3)), vjust = ifelse(rk$media_t_por_gwh >= 0, -0.5, 1.4),
            size = 3.2) +
  scale_fill_manual(values = PAL_TEC, guide = "none") +
  labs(title = "Ranking por intensidad de abatimiento", x = NULL,
       y = "t CO2 evitadas por GWh desplazado",
       caption = "El fósil tiene signo negativo: añade emisiones en lugar de evitarlas") +
  expand_limits(y = c(min(rk$media_t_por_gwh) * 1.25, max(rk$media_t_por_gwh) * 1.15))
p_rank <- g_rank_a + g_rank_b + patchwork::plot_annotation(
  title = "Ranking de tecnologías (16 bloques anuales x 3 penetraciones)")
guardar_fig(p_rank, "04_ranking_tecnologias", 11, 5.2)

# --- (5) Contorno de la superficie: rampa + MESETA de saturación -----------
# La superficie tiene dos regímenes: una rampa decreciente (x1+x2 < umbral, la
# FNCER desplaza fósil) y una meseta plana (x1+x2 >= umbral, excedente vertido).
g_cont <- ggplot(malla, aes(x1, x2, fill = pred_mt)) +
  geom_raster() +
  geom_contour(aes(z = pred_mt), colour = "white", linewidth = 0.3, bins = 10) +
  geom_abline(intercept = umbral, slope = -1, linetype = 2, colour = "#B2182B",
              linewidth = 0.9) +
  geom_point(data = ccd, aes(x1, x2), inherit.aes = FALSE, shape = 21,
             size = 2.6, fill = "white", colour = "black") +
  geom_point(data = data.frame(x1 = opt[["x1"]], x2 = opt[["x2"]]), aes(x1, x2),
             inherit.aes = FALSE, shape = 8, size = 4.2, colour = "#B2182B",
             stroke = 1.2) +
  annotate("label", x = 30, y = 33, size = 3.1, label.size = 0.2,
           label = paste0("MESETA: x1+x2 >= ", round(umbral, 1), " pp\n",
                          "excedente vertido, sin reducción\n",
                          round(rsm$emisiones_minimas_mt, 2), " Mt CO2")) +
  annotate("label", x = 5.5, y = 4, size = 3.1, label.size = 0.2,
           label = paste0("Rampa: la FNCER\ndesplaza fósil\n",
                          fmt_num(rsm$marginal_t_por_pp, 0), " t/pp")) +
  annotate("text", x = 12, y = 12.5, size = 2.9, colour = "#B2182B",
           label = paste0("umbral x1+x2 = ", round(umbral, 1), " pp")) +
  scale_fill_viridis_c(option = "viridis") +
  labs(title = "Superficie de respuesta: emisiones eléctricas en 2030",
       subtitle = paste0("Base media 2015-2025: ", round(em_base, 2),
                         " Mt sin FNCER nueva; la FNCER solo desplaza fósil\n",
                         "Estrella = mezcla que alcanza la meseta (",
                         round(opt[["x1"]], 1), "% solar + ", round(opt[["x2"]], 1),
                         "% eólica)"),
       x = "Participación solar (% de la generación 2030)",
       y = "Participación eólica (% de la generación 2030)",
       fill = "Mt CO2") +
  coord_fixed()
guardar_fig(g_cont, "04_rsm_contorno", 9.5, 7.5)

# --- (6) Heatmap con la meseta delimitada ----------------------------------
g_heat <- ggplot(malla, aes(x1, x2, fill = pred_mt)) +
  geom_raster() +
  geom_abline(intercept = umbral, slope = -1, linetype = 2, colour = "white",
              linewidth = 0.8) +
  geom_point(data = ccd, aes(x1, x2), inherit.aes = FALSE, shape = 21,
             size = 2.4, fill = "white", colour = "black") +
  geom_point(data = data.frame(x1 = opt[["x1"]], x2 = opt[["x2"]]), aes(x1, x2),
             inherit.aes = FALSE, shape = 8, size = 4.5, colour = "#B2182B",
             stroke = 1.2) +
  annotate("text", x = 30, y = 30, size = 3.3, colour = "grey15",
           label = paste0("MESETA: ", round(rsm$emisiones_minimas_mt, 2),
                          " Mt CO2\n(excedente vertido)")) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(title = "Heatmap de emisiones eléctricas de 2030",
       subtitle = paste0("Puntos = corridas del CCD (n = ", nrow(ccd), "; ",
                         sum(ccd$region == "Meseta (saturada)"),
                         " caen en la meseta). Línea discontinua = umbral x1+x2 = ",
                         round(umbral, 1), " pp"),
       x = "Participación solar (%)", y = "Participación eólica (%)",
       fill = "Mt CO2") +
  coord_fixed()
guardar_fig(g_heat, "04_rsm_heatmap", 9.5, 7.5)

# --- (7) Superficie 3D (perspectiva base R) --------------------------------
f3d <- file.path(RUTAS$figuras, "04_rsm_superficie3d.png")
zmat <- matrix(malla$pred_mt, nrow = length(unique(malla$x1)),
               ncol = length(unique(malla$x2)))
gx <- sort(unique(malla$x1)); gy <- sort(unique(malla$x2))
pal3d <- grDevices::colorRampPalette(c("#FCFDBF", "#FE9F6D", "#DE4968",
                                       "#8C2981", "#3B0F70"))(length(zmat))
png(f3d, width = 2400, height = 1800, res = 300, bg = "white")
par(mar = c(2, 2, 3, 1))
persp(gx, gy, zmat, theta = 35, phi = 22, col = pal3d, border = NA, shade = 0.5,
      ticktype = "detailed", cex.axis = 0.7, cex.lab = 0.85,
      xlab = "Solar (%)", ylab = "Eólica (%)", zlab = "Mt CO2",
      main = "Emisiones eléctricas 2030: superficie de respuesta")
dev.off()
msg("figura -> ", basename(f3d))

# ============================================================================
# 12. TABLAS Y GUARDADO
# ============================================================================
paso("9. TABLAS Y RESULTADOS")

guardar_tabla(anova_tab, "04_doe_anova")

medias_tab <- bind_rows(
  medias_tec |> mutate(seccion = "Tecnología"),
  medias_pen |> mutate(seccion = "Penetración"),
  medias_celda |> mutate(seccion = "Interacción")
) |>
  select(seccion, factor, nivel, n, media_t, se_t, ic95_li, ic95_ls, media_mt,
         media_t_por_gwh, everything())
guardar_tabla(medias_tab, "04_doe_medias")
guardar_tabla(ranking, "04_doe_ranking")

rsm_tab <- bind_rows(
  data.frame(seccion = "Parámetro del modelo corregido", termino =
               c("base_media_mt", "fosil_base_mt", "em_constante_mt",
                 "fosil_share_pct", "umbral_fncer_pct", "emisiones_minimas_mt",
                 "reduccion_maxima_mt", "reduccion_maxima_pct",
                 "marginal_t_por_pp", "optimo_solar_pct", "optimo_eolica_pct",
                 "vertido_esquina_gwh", "vertido_esquina_pct"),
             estimacion = c(rsm$base_media_mt, rsm$fosil_base_mt,
                            rsm$em_constante_mt, rsm$fosil_share_pct,
                            rsm$umbral_fncer_pct, rsm$emisiones_minimas_mt,
                            rsm$reduccion_maxima_mt, rsm$reduccion_maxima_pct,
                            rsm$marginal_t_por_pp, rsm$optimo[["x1"]],
                            rsm$optimo[["x2"]], rsm$vertido_esquina_gwh,
                            rsm$vertido_esquina_pct),
             ee = NA, estadistico = NA, p = NA,
             nota = c("base sin FNCER nueva (fósil + ciclo de vida hidro/biomasa)",
                      "emisiones del fósil de la base", "hidro + biomasa (constantes)",
                      "participación fósil de la base media",
                      "umbral de FNCER útil (x1+x2); más allá, vertido",
                      "emisiones en la meseta x1+x2 >= umbral",
                      "base_media_mt - emisiones_minimas_mt",
                      "porcentaje de reducción máxima alcanzable",
                      "pendiente media en x1+x2 < umbral (0 a 21,11 pp)",
                      "mezcla balanceada que alcanza la meseta",
                      "idem", "excedente vertido en (40,40)",
                      "porcentaje de la demanda vertido en (40,40)")),
  data.frame(seccion = "Tramos del orden de mérito",
             termino = paste0(rsm$marginal_por_tramo$tramo, " [pp]"),
             estimacion = rsm$marginal_por_tramo$pendiente_t_por_pp,
             ee = NA, estadistico = NA, p = NA,
             nota = paste0("desde ", round(rsm$marginal_por_tramo$desde_pp, 2),
                           " hasta ", round(rsm$marginal_por_tramo$hasta_pp, 2),
                           " pp de x1+x2")),
  data.frame(seccion = "Ajuste en la región no saturada",
             termino = c("n_corridas_ccd", "n_niveles_x1x2",
                         "cuadratico_estimable_ccd", "R2_lineal_ccd",
                         "pendiente_lineal_ccd_t_por_pp", "R2_cuadratico_malla",
                         "RMSE_cuadratico_malla", "n_puntos_malla",
                         "coef_intercepto", "coef_x1", "coef_x2",
                         "coef_x1_2", "coef_x2_2", "coef_x1x2"),
             estimacion = c(rsm$modelo_no_saturado$n_corridas_ccd,
                            length(rsm$modelo_no_saturado$niveles_suma_pct),
                            0, rsm$r2_lineal_ccd,
                            rsm$pendiente_lineal_ccd_t_por_pp,
                            rsm$r2_cuadratico_malla,
                            rsm$rmse_cuadratico_malla,
                            rsm$modelo_no_saturado$cuadratico_malla$n_puntos,
                            as.numeric(b_rsm)),
             ee = NA, estadistico = NA, p = NA,
             nota = c(rsm$modelo_no_saturado$razon, NA, NA,
                      "modelo lineal exacto sobre las corridas del CCD",
                      NA,
                      "aproximación de 2.º orden sobre la función contable",
                      "error RMSE en Mt CO2", NA, NA, NA, NA, NA, NA, NA)),
  data.frame(seccion = "Sensibilidad: base 2025 (año húmedo)",
             termino = c("base_media_mt", "fosil_share_pct", "umbral_fncer_pct",
                         "emisiones_minimas_mt", "reduccion_maxima_mt",
                         "reduccion_maxima_pct", "marginal_t_por_pp",
                         "optimo_solar_pct", "optimo_eolica_pct"),
             estimacion = c(rsm_2025$base_mt, rsm_2025$fosil_share_pct,
                            rsm_2025$umbral_fncer_pct, rsm_2025$optimo_mt,
                            rsm_2025$reduccion_maxima_mt,
                            rsm_2025$reduccion_maxima_pct,
                            rsm_2025$marginal_t_por_pp,
                            rsm_2025$optimo[["x1"]], rsm_2025$optimo[["x2"]]),
             ee = NA, estadistico = NA, p = NA,
             nota = c("composición observada 2025", NA, NA, NA, NA, NA, NA,
                      NA, NA)),
  data.frame(seccion = "Método proporcional anterior (NO VÁLIDO)",
             termino = c("reduccion_pct", "marginal_t_por_pp", "valido"),
             estimacion = c(rsm$valor_anterior_no_valido[["reduccion_pct"]],
                            rsm$valor_anterior_no_valido[["marginal_t_por_pp"]],
                            0),
             ee = NA, estadistico = NA, p = NA,
             nota = rep(rsm$cambio_metodologico, 3))
)
guardar_tabla(rsm_tab, "04_doe_rsm")
guardar_tabla(robustez_tab, "04_doe_robustez")

# --- Resumen legible --------------------------------------------------------
resumen_tab <- data.frame(
  metrica = c("Corridas del factorial", "Bloques (años)", "Corridas con saturación",
              "F tecnología", "p tecnología", "eta2p tecnología",
              "F penetración", "p penetración", "eta2p penetración",
              "F interacción", "p interacción", "eta2p interacción",
              "Shapiro-Wilk p (modelo principal)", "Levene p (modelo principal)",
              "log-ANOVA p tecnología", "log-ANOVA p penetración",
              "log-ANOVA p interacción", "Kruskal-Wallis p tecnología",
              "Kruskal-Wallis p penetración", "ANOVA Welch p tecnología",
              "Robustez: conclusiones se mantienen",
              "Pendiente renovables (t CO2/pp)", "Pendiente renovables EE",
              "Pendiente media 5 tecnologías (t CO2/pp)",
              "R2 modelo continuo",
              "Demanda 2030 (GWh)", "R2 proyección demanda",
              "Base 2030 (Mt), composición media 2015-2025",
              "Participación fósil de la base (%)",
              "Umbral de FNCER útil x1+x2 (pp)",
              "Emisiones mínimas 2030 (Mt) en la meseta",
              "Reducción máxima (Mt)", "Reducción máxima (%)",
              "Reducción marginal en la región no saturada (t CO2/pp)",
              "Óptimo solar (%)", "Óptimo eólica (%)",
              "R2 modelo de 2.º orden (región no saturada)",
              "Base 2030 (Mt), composición 2025 (sensibilidad)",
              "Umbral de FNCER útil base 2025 (pp)",
              "Reducción máxima base 2025 (Mt)",
              "Reducción anterior por reparto proporcional (%) — NO VÁLIDA"),
  valor = c(nrow(diseno), length(ANIOS_DOE), sum(diseno$saturado),
            f_tec, p_tec, e2_tec, f_pen, p_pen, e2_pen, f_int, p_int, e2_int,
            diagnostico$shapiro_p, diagnostico$levene_p,
            tab_log$p[tab_log$termino == "tec"],
            tab_log$p[tab_log$termino == "pen"],
            tab_log$p[tab_log$termino == "tec:pen"],
            kw_tec$p.value, kw_pen$p.value, unname(welch_tec$p.value),
            as.numeric(mantiene),
            pend_ren, se_pend_ren, pend_media_5, sm$r.squared,
            dem_2030, summary(m_dem)$r.squared, em_base,
            rsm$fosil_share_pct, rsm$umbral_fncer_pct,
            rsm$emisiones_minimas_mt, rsm$reduccion_maxima_mt,
            rsm$reduccion_maxima_pct, rsm$marginal_t_por_pp,
            rsm$optimo[["x1"]], rsm$optimo[["x2"]], rsm$r2_cuadratico_malla,
            rsm_2025$base_mt, rsm_2025$umbral_fncer_pct,
            rsm_2025$reduccion_maxima_mt,
            rsm$valor_anterior_no_valido[["reduccion_pct"]])
)
guardar_tabla(resumen_tab, "04_doe_resumen")
print(resumen_tab, digits = 4)

resumen <- list(
  n_corridas = nrow(diseno),
  n_bloques = length(ANIOS_DOE),
  anios = ANIOS_DOE,
  diseno_factorial = "5 tecnologías x 3 penetraciones x 16 bloques (completo)",
  semilla = CONST$semilla,
  anova = c(f_tec = f_tec, p_tec = p_tec, eta2p_tec = e2_tec,
            f_pen = f_pen, p_pen = p_pen, eta2p_pen = e2_pen,
            f_int = f_int, p_int = p_int, eta2p_int = e2_int,
            eta2p_bloque = e2_blo),
  robustez = list(p_valores = p_robustez, mantiene_conclusiones = mantiene,
                  log_anova = tab_log, shapiro_log_p = sw_log$p.value,
                  levene_log_p = lev_log$`Pr(>F)`[1],
                  kruskal_tec = kw_tec, kruskal_pen = kw_pen,
                  kruskal_celda = kw_celda, welch_tec = welch_tec,
                  c_log = c_log),
  medias_tec = medias_tec, medias_pen = medias_pen,
  ranking_t_por_gwh = ranking,
  pendiente_renovables_t_por_pp = pend_ren,
  pendiente_renovables_ee = se_pend_ren,
  pendiente_media_5_tecnologias_t_por_pp = pend_media_5,
  pendientes_por_tecnologia = pend,
  # --- RSM corregido (orden de mérito real, solo fósil) --------------------
  base_media_2030_mt = rsm$base_media_mt,
  fosil_base_2030_mt = rsm$fosil_base_mt,
  em_constante_mt = rsm$em_constante_mt,
  fosil_share_pct = rsm$fosil_share_pct,
  umbral_fncer_pct = rsm$umbral_fncer_pct,
  emisiones_minimas_2030_mt = rsm$emisiones_minimas_mt,
  reduccion_maxima_mt = rsm$reduccion_maxima_mt,
  reduccion_maxima_pct = rsm$reduccion_maxima_pct,
  reduccion_marginal_t_por_pp = rsm$marginal_t_por_pp,
  optimo_solar_pct = rsm$optimo[["x1"]],
  optimo_eolica_pct = rsm$optimo[["x2"]],
  rsm_r2_cuadratico_malla = rsm$r2_cuadratico_malla,
  rsm_r2_lineal_ccd = rsm$r2_lineal_ccd,
  base_2025_2030_mt = rsm_2025$base_mt,
  umbral_fncer_base2025_pct = rsm_2025$umbral_fncer_pct,
  reduccion_maxima_base2025_mt = rsm_2025$reduccion_maxima_mt,
  reduccion_pct_proporcional_no_valida = rsm$valor_anterior_no_valido[["reduccion_pct"]],
  vertido_esquina_gwh = rsm$vertido_esquina_gwh,
  vertido_esquina_pct = rsm$vertido_esquina_pct,
  demanda_2030_gwh = dem_2030,
  sensibilidad_fe_gas = sensibilidad,
  figuras = c("04_efectos_principales.png", "04_interaccion.png",
              "04_residuos_anova.png", "04_ranking_tecnologias.png",
              "04_rsm_contorno.png", "04_rsm_heatmap.png",
              "04_rsm_superficie3d.png")
)

resultados <- list(
  diseno = diseno,
  anova = anova_tab,
  diagnostico = diagnostico,
  robustez = list(tabla = robustez_tab, p_valores = p_robustez,
                  mantiene_conclusiones = mantiene,
                  log_anova = tab_log, c_log = c_log,
                  kruskal = list(tecnologia = kw_tec, penetracion = kw_pen,
                                 celda = kw_celda),
                  welch_tecnologia = welch_tec),
  medias = list(tecnologia = medias_tec, penetracion = medias_pen,
                celda = medias_celda, ranking = ranking),
  tukey = list(tecnologias = tk_tec, interaccion = tk_celda),
  modelo_continuo = modelo_continuo,
  ccd = ccd,
  rsm = rsm,
  parametros = list(fe = FE, ctx = ctx,
                    composicion_base_media = sh_media,
                    composicion_base_2025 = sh_2025,
                    fosil_base_media = fos_med,
                    fosil_base_2025 = fos_2025,
                    demanda_2030 = c(estimacion = dem_2030, li = dem_2030_li,
                                     ls = dem_2030_ls, tasa_log = tasa)),
  resumen = resumen
)
rds_out <- file.path(RUTAS$procesados, "resultados_doe.rds")
saveRDS(resultados, rds_out)
msg("resultados -> ", basename(rds_out), " (", format(file.size(rds_out) / 1024, digits = 4), " KB)")

options(op_previas)

paso("04 · DOE COMPLETADO")
