# ============================================================================
# 01b_eda.R — Análisis exploratorio de datos (EDA) del proyecto
# ----------------------------------------------------------------------------
# Produce las figuras y tablas descriptivas que sostienen el paper:
#   - Serie nacional de CO2 y su composición por combustible (1990-2024)
#   - Matriz eléctrica y penetración FNCER (2000-2025)
#   - Emisiones evitadas por FNCER
#   - Distribución espacial y concentración de emisiones (grilla 0,1°)
#   - Perfil regional y proyectos FNCER
# Salidas: output/figures, output/maps, output/tables, datos/procesados/eda_resumen.rds
# ============================================================================
source(file.path("scripts", "00_setup.R"))

panel <- readRDS(file.path(RUTAS$procesados, "panel_nacional.rds"))
grilla <- readRDS(file.path(RUTAS$procesados, "grilla_co2.rds"))
municipios <- readRDS(file.path(RUTAS$procesados, "municipios.rds"))
proy <- readRDS(file.path(RUTAS$procesados, "proyectos_fncer.rds"))
adm1 <- st_read(file.path(RUTAS$geometrias, "gadm41_COL_1.json"), quiet = TRUE)

paso("1. SERIE NACIONAL DE CO2")
serie <- panel |>
  filter(anio >= 1990, anio <= 2024) |>
  select(anio, co2_owid, co2_carbon, co2_petroleo, co2_gas, co2_cemento, co2_quema,
         ghg_total_owid, poblacion, pib)

largo <- serie |>
  select(anio, Carbon = co2_carbon, Petroleo = co2_petroleo, Gas = co2_gas,
         Cemento = co2_cemento, Quema = co2_quema) |>
  pivot_longer(-anio, names_to = "combustible", values_to = "mt") |>
  mutate(combustible = factor(combustible,
                              levels = c("Carbon", "Petroleo", "Gas", "Cemento", "Quema")))

p1 <- ggplot(largo, aes(anio, mt, fill = combustible)) +
  geom_area(alpha = .9, colour = "white", linewidth = .2) +
  geom_line(data = serie, aes(anio, co2_owid), inherit.aes = FALSE,
            colour = "black", linewidth = .9) +
  scale_fill_brewer(palette = "YlOrBr", direction = -1, name = "Combustible") +
  labs(title = "Emisiones de CO2 de Colombia, 1990-2024",
       subtitle = "Área: composición por fuente (Global Carbon Budget). Línea negra: total nacional",
       x = NULL, y = "Mt CO2 / año",
       caption = "Fuente: OWID / Global Carbon Budget (2026)") +
  scale_x_continuous(breaks = seq(1990, 2024, 5)) +
  tema()
guardar_fig(p1, "eda_01_serie_co2", 9, 5.5)

p1b <- serie |>
  mutate(indice = 100 * co2_owid / co2_owid[anio == 1990]) |>
  filter(anio >= 1990) |>
  ggplot(aes(anio, indice)) +
  geom_hline(yintercept = 100, linetype = 2, colour = "grey50") +
  geom_line(colour = "#B2182B", linewidth = 1) +
  geom_point(size = 1.5, colour = "#B2182B") +
  annotate("text", x = 1991.5, y = 165, hjust = 0, size = 3.3,
           label = sprintf("2019: %.0f (máximo)", 100 * serie$co2_owid[serie$anio == 2019] / serie$co2_owid[serie$anio == 1990])) +
  labs(title = "Índice de emisiones de CO2 (1990 = 100)",
       subtitle = "El crecimiento se desacopla parcialmente del PIB desde 2015",
       x = NULL, y = "Índice (1990 = 100)") +
  tema()
guardar_fig(p1b, "eda_02_indice_co2", 9, 4.5)

paso("2. MATRIZ ELÉCTRICA Y FNCER")
mix <- panel |>
  filter(anio >= 2000, anio <= 2025, meses_reportados == 12) |>
  select(anio, Hidro = gen_hidro_gwh, Biomasa = gen_biomasa_gwh, Solar = gen_solar_gwh,
         Eolica = gen_eolica_gwh, Fosil = gen_fosil_gwh, Otros = gen_otros_gwh) |>
  pivot_longer(-anio, names_to = "fuente", values_to = "gwh") |>
  mutate(fuente = factor(fuente, levels = c("Hidro", "Biomasa", "Solar", "Eolica", "Fosil", "Otros")))

p2 <- ggplot(mix, aes(anio, gwh, fill = fuente)) +
  geom_col(width = .85) +
  scale_fill_manual(values = c(Hidro = "#1F78B4", Biomasa = "#7A9E3B", Solar = "#F4C430",
                               Eolica = "#4A90D9", Fosil = "#5B5B5B", Otros = "#C7C7C7"),
                    name = NULL) +
  labs(title = "Generación eléctrica por fuente en Colombia, 2000-2025",
       subtitle = "La hidroeléctrica domina la matriz; solar y eólica crecen desde 2017",
       x = NULL, y = "GWh / año", caption = "Fuente: UPME (generación por fuente)") +
  tema()
guardar_fig(p2, "eda_03_mix_electrico", 9, 5.5)

p3 <- panel |>
  filter(anio >= 2000, anio <= 2025, meses_reportados == 12) |>
  select(anio, `FNCER no convencional` = pct_fncer_nc, `Renovable total (incl. hidro)` = pct_renovable) |>
  pivot_longer(-anio) |>
  ggplot(aes(anio, value, colour = name)) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_colour_manual(values = c("FNCER no convencional" = "#F4C430",
                                 "Renovable total (incl. hidro)" = "#1F78B4"), name = NULL) +
  labs(title = "Participación de las renovables en la generación eléctrica",
       subtitle = "Porcentaje de la generación total (UPME)",
       x = NULL, y = "% de la generación") +
  tema()
guardar_fig(p3, "eda_04_pct_renovable", 9, 4.5)

p4 <- panel |>
  filter(anio >= 2010, anio <= 2025) |>
  ggplot(aes(anio, emisiones_evitadas_mt)) +
  geom_col(fill = "#1B9E77", width = .8) +
  geom_text(aes(label = ifelse(anio %% 3 == 0 | anio >= 2022, round(emisiones_evitadas_mt, 2), "")),
            vjust = -.4, size = 2.9) +
  labs(title = "Emisiones de CO2 evitadas por FNCER no convencionales",
       subtitle = "GWh FNCER × intensidad de emisión del sistema eléctrico",
       x = NULL, y = "Mt CO2 evitadas / año",
       caption = "Cálculo propio con generación UPME e intensidad Ember") +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  tema()
guardar_fig(p4, "eda_05_evitadas", 9, 4.8)

paso("3. DISTRIBUCIÓN ESPACIAL Y CONCENTRACIÓN")

p5 <- ggplot(grilla) +
  geom_sf(aes(fill = dens_2023), colour = NA, linewidth = 0) +
  geom_sf(data = adm1, fill = NA, colour = "white", linewidth = .25, alpha = .6) +
  scale_fill_viridis_c(trans = "log1p", option = "inferno", name = "t CO2/km²",
                       breaks = c(0, 1, 10, 100, 1000, 10000),
                       labels = c("0", "1", "10", "100", "1.000", "10.000")) +
  labs(title = "Densidad de emisiones de CO2 por celda de 0,1°, 2023",
       subtitle = "Escala logarítmica; fronteras departamentales en blanco",
       caption = "Fuente: EDGAR v2026 (JRC), celdas de ~11 km") +
  tema() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
guardar_fig(p5, "eda_06_mapa_densidad", 8.5, 9)

p6 <- ggplot(grilla) +
  geom_sf(aes(fill = d_co2), colour = NA, linewidth = 0) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                       trans = "pseudo_log", midpoint = 0,
                       name = "Δ t CO2",
                       labels = label_number(big.mark = " ")) +
  labs(title = "Cambio de emisiones entre 2015 y 2023",
       subtitle = "Rojo: aumento; azul: reducción (escala pseudo-logarítmica)",
       caption = "Fuente: EDGAR v2026") +
  tema() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())
guardar_fig(p6, "eda_07_mapa_cambio", 8.5, 9)

p7 <- ggplot(grilla, aes(region, dens_2023 + 1, fill = region)) +
  geom_boxplot(outlier.alpha = .15, show.legend = FALSE) +
  scale_y_log10(labels = label_number(big.mark = " ")) +
  scale_fill_manual(values = PAL_REGION) +
  labs(title = "Densidad de emisiones por región natural (2023)",
       subtitle = "t CO2/km² + 1, escala logarítmica",
       x = NULL, y = "t CO2/km² (log)") +
  tema()
guardar_fig(p7, "eda_08_boxplot_region", 8.5, 4.8)

# Curva de concentración (Lorenz) de emisiones entre celdas
ord <- sort(grilla$co2_2023, decreasing = TRUE)
lorenz <- tibble(
  prop_celdas = seq_along(ord) / length(ord),
  prop_emis = cumsum(ord) / sum(ord)
)
# Gini sobre la curva ascendente (la curva del gráfico está en orden descendente)
x_asc <- sort(grilla$co2_2023); p_asc <- seq_along(x_asc) / length(x_asc)
L_asc <- cumsum(x_asc) / sum(x_asc)
gini <- 1 - 2 * sum((L_asc[-1] + L_asc[-length(L_asc)]) / 2 * diff(p_asc))
p8 <- ggplot(lorenz, aes(prop_celdas, prop_emis)) +
  geom_abline(linetype = 2, colour = "grey60") +
  geom_line(colour = "#B2182B", linewidth = 1.1) +
  scale_x_continuous(labels = percent_format()) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Concentración espacial de las emisiones de CO2 (2023)",
       subtitle = sprintf("El 1%% de las celdas concentra el %.1f%% de las emisiones; Gini = %.3f",
                          100 * lorenz$prop_emis[round(0.01 * nrow(lorenz))], gini),
       x = "% de celdas (ordenadas de mayor a menor emisión)",
       y = "% acumulado de emisiones") +
  tema()
guardar_fig(p8, "eda_09_lorenz", 7.5, 5.5)

paso("4. PERFIL REGIONAL Y PROYECTOS FNCER")
tabla_region <- grilla |>
  st_drop_geometry() |>
  group_by(region) |>
  summarise(
    celdas = n(),
    area_km2 = sum(area_km2),
    co2_2015_mt = sum(co2_2015) / 1e6,
    co2_2023_mt = sum(co2_2023) / 1e6,
    share_2023 = 100 * sum(co2_2023) / sum(grilla$co2_2023),
    dens_media = mean(dens_2023),
    dens_mediana = median(dens_2023),
    poblacion_m = sum(poblacion_2020) / 1e6,
    co2_pc = sum(co2_2023) / sum(poblacion_2020),
    var_pct = 100 * (sum(co2_2023) / sum(co2_2015) - 1)
  ) |>
  arrange(desc(co2_2023_mt)) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))
guardar_tabla(tabla_region, "01_eda_region")

p9 <- tabla_region |>
  select(region, `2015` = co2_2015_mt, `2023` = co2_2023_mt) |>
  pivot_longer(-region) |>
  ggplot(aes(reorder(region, -value), value, fill = name)) +
  geom_col(position = "dodge", width = .7) +
  scale_fill_manual(values = c("2015" = "#9ECAE1", "2023" = "#08519C"), name = NULL) +
  labs(title = "Emisiones por región natural: 2015 vs 2023",
       x = NULL, y = "Mt CO2 / año") +
  tema()
guardar_fig(p9, "eda_10_region_barras", 8.5, 4.8)

proy_tab <- proy |>
  filter(!is.na(mun_id), year <= 2025) |>
  count(year, tipo, wt = capacidad_mw, name = "mw") |>
  tidyr::complete(year, tipo, fill = list(mw = 0))
p10 <- ggplot(proy_tab, aes(year, mw, fill = tipo)) +
  geom_col(width = .8) +
  scale_fill_manual(values = c(Solar = "#F4C430", `Eólico` = "#4A90D9"), name = NULL) +
  labs(title = "Capacidad FNCER declarada por año de entrada (MinEnergía)",
       subtitle = "Suma de la capacidad de los proyectos solares y eólicos",
       x = NULL, y = "MW declarados") +
  tema()
guardar_fig(p10, "eda_11_proyectos_fncer", 8.5, 4.5)

# Scatter CO2 vs PIB (desacoplamiento)
p11 <- serie |>
  filter(anio >= 1990) |>
  ggplot(aes(pib / 1e12, co2_owid)) +
  geom_path(colour = "grey60", arrow = arrow(length = unit(2, "mm"), type = "closed")) +
  geom_point(aes(colour = anio), size = 2.2) +
  scale_colour_viridis_c(name = "Año") +
  labs(title = "CO2 vs PIB: trayectoria de desacoplamiento",
       subtitle = "Cada punto es un año (1990-2024); la flecha indica el sentido temporal",
       x = "PIB (billones de USD PPP 2011)", y = "Mt CO2") +
  tema()
guardar_fig(p11, "eda_12_co2_pib", 8, 5.5)

paso("5. TABLA DESCRIPTIVA Y GUARDADO")
descriptivos <- tibble::tribble(
  ~indicador, ~anio, ~valor, ~unidad,
  "CO2 nacional", 1990, round(serie$co2_owid[serie$anio == 1990], 1), "Mt CO2",
  "CO2 nacional", 2015, round(serie$co2_owid[serie$anio == 2015], 1), "Mt CO2",
  "CO2 nacional (máximo)", 2019, round(serie$co2_owid[serie$anio == 2019], 1), "Mt CO2",
  "CO2 nacional", 2023, round(serie$co2_owid[serie$anio == 2023], 1), "Mt CO2",
  "CO2 nacional", 2024, round(serie$co2_owid[serie$anio == 2024], 1), "Mt CO2",
  "GEI totales", 2024, round(serie$ghg_total_owid[serie$anio == 2024], 1), "Mt CO2e",
  "CO2 per cápita", 2024, round(panel$co2_pc[panel$anio == 2024], 2), "t CO2/hab",
  "Generación eléctrica", 2024, round(panel$gen_total_gwh[panel$anio == 2024] / 1000, 1), "TWh",
  "Participación FNCER no convencional", 2024, round(panel$pct_fncer_nc[panel$anio == 2024], 2), "% generación",
  "Emisiones evitadas por FNCER", 2024, round(panel$emisiones_evitadas_mt[panel$anio == 2024], 2), "Mt CO2",
  "Celdas del marco muestral", NA, nrow(grilla), "celdas 0,1°",
  "Concentración (top 1% celdas)", 2023, round(100 * sum(head(ord, 97)) / sum(ord), 1), "% emisiones",
  "Municipios con proyecto FNCER", 2023, sum(municipios$tratado), "municipios"
)
guardar_tabla(descriptivos, "01_eda_descriptivos")

eda_resumen <- list(
  serie = serie, tabla_region = tabla_region, descriptivos = descriptivos,
  gini = gini, concentracion = c(top1 = sum(head(ord, 97)) / sum(ord),
                                 top5 = sum(head(ord, 483)) / sum(ord),
                                 top10 = sum(head(ord, 966)) / sum(ord)),
  mix_2024 = panel |> filter(anio == 2024) |>
    select(gen_total_gwh, gen_hidro_gwh, gen_fosil_gwh, gen_solar_gwh,
           gen_eolica_gwh, gen_biomasa_gwh, pct_fncer_nc, pct_renovable),
  n_celdas = nrow(grilla), n_municipios = nrow(municipios)
)
saveRDS(eda_resumen, file.path(RUTAS$procesados, "eda_resumen.rds"))

paso("EDA COMPLETADO")
cat("CO2 1990:", round(serie$co2_owid[serie$anio == 1990], 1),
    "| 2019:", round(serie$co2_owid[serie$anio == 2019], 1),
    "| 2024:", round(serie$co2_owid[serie$anio == 2024], 1), "Mt\n")
cat("Concentración top1%:", round(100 * eda_resumen$concentracion["top1"], 1), "%\n")
print(tabla_region)
