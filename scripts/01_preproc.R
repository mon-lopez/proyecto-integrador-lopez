# ============================================================================
# 01_preproc.R — Construcción de la capa de datos del proyecto
# ----------------------------------------------------------------------------
# Genera, a partir de fuentes primarias verificables:
#   (1) grilla_co2.rds       : celdas EDGAR 0,1° sobre Colombia (marco muestral)
#   (2) panel_nacional.rds   : panel anual 1990-2026 (CO2, energía, FNCER, GEI)
#   (3) municipios.rds       : panel municipal 2015->2023 (unidad cuasi-experimental)
#   (4) proyectos_fncer.rds  : proyectos FNCER geolocalizados por municipio
#   (5) resumen_cobertura.csv: auditoría de la integración de fuentes
# ============================================================================
source(file.path("scripts", "00_setup.R"))

paso("1. GEOMETRÍAS Y MARCO MUESTRAL (grilla EDGAR 0,1°)")

adm0 <- st_read(file.path(RUTAS$datos, "gadm41_COL_0.json"), quiet = TRUE) |> st_geometry()
adm1 <- st_read(file.path(RUTAS$geometrias, "gadm41_COL_1.json"), quiet = TRUE)
adm2 <- st_read(file.path(RUTAS$geometrias, "gadm41_COL_2.json"), quiet = TRUE)
tierra <- st_union(st_geometry(adm1))
msg("ADM1: ", nrow(adm1), " departamentos | ADM2: ", nrow(adm2), " municipios")
msg("Área terrestre: ", fmt_num(as.numeric(st_area(st_transform(tierra, CONST$crs_proj))) / 1e6, 0), " km²")

# Rásteres EDGAR (0,1°; t CO2/celda/año)
anios_edgar <- c(2015, 2023, 2025)
r_edgar <- lapply(anios_edgar, function(a) {
  rast(file.path(RAIZ, "emisiones_co2_espacial", "datos",
                 sprintf("co2_grid_edgar_%d.tif", a)))
}) |> setNames(paste0("y", anios_edgar))
plantilla <- crop(r_edgar$y2023, ext(vect(tierra)) + 0.2)
msg("Plantilla: ", nrow(plantilla), " x ", ncol(plantilla), " celdas de ", res(plantilla)[1], "°")

# Fracción de territorio colombiano en cada celda (evita sumar celdas de países vecinos)
cob <- rasterize(vect(tierra), plantilla, cover = TRUE, background = 0)
cobv <- values(cob)[, 1]
sel <- which(cobv > 0.01)
xy_cel <- xyFromCell(plantilla, sel)
msg("Celdas con cobertura colombiana > 1%: ", length(sel))

mascara <- ifel(cob > 0.01, 1, NA)
pol <- as.polygons(mascara, na.rm = TRUE, dissolve = FALSE)
grilla <- st_as_sf(pol)[, "layer"]
grilla$celda_edgar <- sel
grilla$layer <- NULL
stopifnot(nrow(grilla) == length(sel))
chk_pol <- suppressWarnings(st_coordinates(st_centroid(st_geometry(grilla[1:5, ]))))
stopifnot(max(abs(chk_pol - xy_cel[1:5, ])) < 1e-6)

# Emisiones por año (valor EDGAR íntegro de cada celda del marco)
for (a in anios_edgar) {
  idx <- cellFromXY(r_edgar[[paste0("y", a)]], xy_cel)
  v <- r_edgar[[paste0("y", a)]][idx][, 1]
  v[is.na(v)] <- 0
  grilla[[paste0("co2_", a)]] <- v
}
grilla$cobertura <- cobv[sel]
msg("Emisiones EDGAR atribuidas: 2015 = ", fmt_num(sum(grilla$co2_2015) / 1e6, 2),
    " Mt | 2023 = ", fmt_num(sum(grilla$co2_2023) / 1e6, 2),
    " Mt | 2025 = ", fmt_num(sum(grilla$co2_2025) / 1e6, 2), " Mt")

# Geometría proyectada (EPSG:9377) y áreas
grilla <- st_transform(grilla, CONST$crs_proj)
grilla$area_km2 <- as.numeric(st_area(grilla)) / 1e6
cent <- st_centroid(st_geometry(grilla))
grilla$x_9377 <- st_coordinates(cent)[, 1]
grilla$y_9377 <- st_coordinates(cent)[, 2]

# Centroides geográficos para uniones
cent_geo <- st_transform(cent, CONST$crs_geo)
grilla$lon <- st_coordinates(cent_geo)[, 1]
grilla$lat <- st_coordinates(cent_geo)[, 2]

# Departamento y municipio (centroide; respaldo: vecino más cercano)
adm2 <- st_transform(adm2, CONST$crs_geo)
adm1 <- st_transform(adm1, CONST$crs_geo)
idx_dep <- st_intersects(cent_geo, adm1, sparse = TRUE)
idx_mun <- st_intersects(cent_geo, adm2, sparse = TRUE)
falta <- lengths(idx_mun) == 0
if (any(falta)) {
  idx_mun[falta] <- st_nearest_feature(cent_geo[falta], adm2)
  idx_dep[falta] <- st_nearest_feature(cent_geo[falta], adm1)
  msg("Celdas sin municipio por centroide (asignadas al más cercano): ", sum(falta))
}
grilla$depto <- adm1$NAME_1[vapply(idx_dep, `[`, integer(1), 1)]
grilla$municipio <- adm2$NAME_2[vapply(idx_mun, `[`, integer(1), 1)]
grilla$mun_id <- adm2$GID_2[vapply(idx_mun, `[`, integer(1), 1)]

# Región natural
normaliza <- function(x) toupper(gsub("[^A-Za-z]", "", iconv(x, to = "ASCII//TRANSLIT")))
region_depto <- c(
  "LAGUAJIRA"="Caribe", "ATLANTICO"="Caribe", "MAGDALENA"="Caribe", "CESAR"="Caribe",
  "CORDOBA"="Caribe", "SUCRE"="Caribe", "BOLIVAR"="Caribe",
  "SANANDRESYPROVIDENCIA"="Caribe", "SANANDRES"="Caribe",
  "ANTIOQUIA"="Andina", "CALDAS"="Andina", "RISARALDA"="Andina", "QUINDIO"="Andina",
  "TOLIMA"="Andina", "HUILA"="Andina", "CUNDINAMARCA"="Andina", "BOYACA"="Andina",
  "SANTANDER"="Andina", "NORTEDESANTANDER"="Andina", "BOGOTADC"="Andina",
  "CHOCO"="Pacifica", "VALLEDELCAUCA"="Pacifica", "CAUCA"="Pacifica", "NARINO"="Pacifica",
  "META"="Orinoquia", "CASANARE"="Orinoquia", "ARAUCA"="Orinoquia", "VICHADA"="Orinoquia",
  "AMAZONAS"="Amazonia", "CAQUETA"="Amazonia", "PUTUMAYO"="Amazonia", "GUAVIARE"="Amazonia",
  "GUAINIA"="Amazonia", "VAUPES"="Amazonia"
)
grilla$region <- unname(region_depto[normaliza(grilla$depto)])
if (anyNA(grilla$region)) {
  warning("Sin región: ", paste(unique(grilla$depto[is.na(grilla$region)]), collapse = ", "))
}
msg("Regiones: ", paste(names(table(grilla$region)), table(grilla$region), sep = "=", collapse = " | "))

# ---- Población: WorldPop 1 km calibrada a los totales de OWID ---------------
# WorldPop se usa por su patrón espacial (urbano/rural); el nivel se calibra
# al total nacional oficial para que la suma coincida con la población del país.
owid_pob <- read.csv(file.path(RUTAS$datos, "owid-co2-data.csv")) |>
  subset(country == "Colombia", select = c(year, population))
for (a in c(2015, 2020)) {
  wp <- rast(file.path(RUTAS$brutos, sprintf("worldpop_col_%d_1km.tif", a)))
  rp <- resample(wp, plantilla, method = "sum")       # agregación 1 km -> 0,1°
  p <- values(rp)[sel, 1]
  p[is.na(p)] <- 0
  total_oficial <- owid_pob$population[owid_pob$year == a]
  factor_cal <- total_oficial / sum(p)
  grilla[[paste0("poblacion_", a)]] <- p * factor_cal
  msg("Población ", a, ": cruda = ", fmt_num(sum(p) / 1e6, 2),
      " M | calibrada a OWID = ", fmt_num(total_oficial / 1e6, 2),
      " M (factor ", round(factor_cal, 3), ")")
  rm(wp, rp); gc(verbose = FALSE)
}

# ---- Climate TRACE 2023: emisiones por celda y macro-sector -----------------
ct <- read.csv(file.path(RUTAS$datos, "climate_trace_COL_2023.csv"), stringsAsFactors = FALSE)
ct <- ct[!is.na(ct$emissionsQuantity) & !is.na(ct$lon) & !is.na(ct$lat), ]
ct$celda_edgar <- cellFromXY(plantilla, as.matrix(ct[, c("lon", "lat")]))
ct <- ct[!is.na(ct$celda_edgar), ]
ct$grupo <- dplyr::case_when(
  ct$sector == "power" ~ "power",
  ct$sector == "transportation" ~ "transporte",
  ct$sector == "buildings" ~ "edificaciones",
  ct$sector == "manufacturing" ~ "industria",
  ct$sector == "fossil-fuel-operations" ~ "ffo",
  ct$sector == "forestry-and-land-use" ~ "lulucf",
  ct$sector == "waste" ~ "residuos",
  TRUE ~ "otros"
)
ct_agg <- ct |>
  group_by(celda_edgar, grupo) |>
  summarise(co2 = sum(emissionsQuantity, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = grupo, values_from = co2, values_fill = 0)
names(ct_agg)[-1] <- paste0("ct_", names(ct_agg)[-1])
grilla <- grilla |> left_join(st_drop_geometry(ct_agg), by = "celda_edgar")
cols_ct <- grep("^ct_", names(grilla), value = TRUE)
for (cc in cols_ct) grilla[[cc]][is.na(grilla[[cc]])] <- 0
msg("Climate TRACE 2023 total: ",
    fmt_num(sum(as.matrix(st_drop_geometry(grilla[, cols_ct]))) / 1e6, 2), " Mt")

# ---- Proyectos FNCER (MinEnergía) y tratamiento municipal -------------------
paso("2. PROYECTOS FNCER Y DEFINICIÓN DEL TRATAMIENTO")
proy <- read.csv(file.path(RUTAS$datos, "minenergia_meta_fncer.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
proy <- proy |>
  mutate(
    tipo = factor(tipo, levels = c("Solar", "Eólico")),
    year = as.integer(substr(fecha_estimada_fpo, 1, 4)),
    capacidad_mw = as.numeric(capacidad),
    emisiones_co2_ton = as.numeric(emisiones_co2_ton_a_o),
    dep_n = normaliza(departamento),
    mun_n = normaliza(municipio)
  ) |>
  filter(!is.na(year), !is.na(capacidad_mw))

adm2$dep_n <- normaliza(adm2$NAME_1)
adm2$mun_n <- normaliza(adm2$NAME_2)
adm2$key <- paste(adm2$dep_n, adm2$mun_n)
proy$key <- paste(proy$dep_n, proy$mun_n)

proy$mun_idx <- NA_integer_
for (i in seq_len(nrow(proy))) {
  hit <- which(adm2$key == proy$key[i])
  if (length(hit)) {
    proy$mun_idx[i] <- hit[1]
  } else {
    cand <- which(adm2$dep_n == proy$dep_n[i])
    if (length(cand)) {
      d <- utils::adist(adm2$mun_n[cand], proy$mun_n[i])
      if (min(d) <= 3) proy$mun_idx[i] <- cand[which.min(d)]
    }
  }
}
msg("Proyectos FNCER: ", nrow(proy), " | geolocalizados: ", sum(!is.na(proy$mun_idx)),
    " (", round(100 * mean(!is.na(proy$mun_idx))), "%)")
proy$mun_id <- adm2$GID_2[proy$mun_idx]
proy$depto_gadm <- adm2$NAME_1[proy$mun_idx]

mun_trat <- proy |>
  filter(year <= 2023, !is.na(mun_id)) |>
  group_by(mun_id) |>
  summarise(n_proy = n(), cap_mw = sum(capacidad_mw),
            n_solar = sum(tipo == "Solar"), n_eolico = sum(tipo == "Eólico"),
            .groups = "drop") |>
  mutate(tratado = 1L)
msg("Municipios tratados (proyecto FNCER ≤ 2023): ", nrow(mun_trat),
    " | capacidad declarada: ", fmt_num(sum(mun_trat$cap_mw), 1), " MW")

grilla <- grilla |>
  left_join(mun_trat, by = "mun_id") |>
  mutate(n_proy = tidyr::replace_na(n_proy, 0), cap_mw = tidyr::replace_na(cap_mw, 0),
         n_solar = tidyr::replace_na(n_solar, 0), n_eolico = tidyr::replace_na(n_eolico, 0),
         tratado = as.integer(tidyr::replace_na(tratado, 0L) > 0))

mun_pol <- st_transform(adm2[adm2$GID_2 %in% mun_trat$mun_id, ], CONST$crs_proj)
d <- st_distance(cent, st_geometry(mun_pol))
grilla$dist_fncer_km <- as.numeric(apply(d, 1, min)) / 1000
rm(d); gc(verbose = FALSE)

# ---- Variables derivadas ---------------------------------------------------
grilla <- grilla |>
  mutate(
    dens_2015 = co2_2015 / area_km2,
    dens_2023 = co2_2023 / area_km2,
    dens_2025 = co2_2025 / area_km2,
    d_co2 = co2_2023 - co2_2015,
    crec_co2 = ifelse(co2_2015 > 0, (co2_2023 - co2_2015) / co2_2015, NA_real_),
    log_co2_2023 = log1p(co2_2023),
    log_dens_2023 = log1p(dens_2023),
    pop_dens_2020 = poblacion_2020 / area_km2,
    urbano = as.integer(pop_dens_2020 > 500),
    co2_pc_2023 = ifelse(poblacion_2020 > 0, co2_2023 / poblacion_2020, NA_real_),
    grupo_emision = factor(dplyr::ntile(dens_2023, 3), levels = 1:3,
                           labels = c("Baja", "Media", "Alta")),
    estrato = paste(region, grupo_emision, sep = "_"),
    region = factor(region, levels = names(PAL_REGION))
  ) |>
  select(celda_edgar, area_km2, cobertura, lon, lat, x_9377, y_9377, depto, municipio, mun_id,
         region, co2_2015, co2_2023, co2_2025, dens_2015, dens_2023, dens_2025,
         d_co2, crec_co2, log_co2_2023, log_dens_2023, poblacion_2015, poblacion_2020,
         pop_dens_2020, urbano, co2_pc_2023, everything())

saveRDS(grilla, file.path(RUTAS$procesados, "grilla_co2.rds"))
msg("grilla_co2.rds: ", nrow(grilla), " celdas x ", ncol(grilla), " variables")

# ---- Panel nacional --------------------------------------------------------
paso("3. PANEL NACIONAL 1990-2026")
owid <- read.csv(file.path(RUTAS$datos, "owid-co2-data.csv"), stringsAsFactors = FALSE)
owid <- owid |>
  filter(country == "Colombia") |>
  transmute(
    anio = year,
    co2_owid = co2, co2_luc_owid = co2_including_luc, ghg_total_owid = total_ghg,
    ch4_owid = methane, n2o_owid = nitrous_oxide,
    co2_carbon = coal_co2, co2_petroleo = oil_co2, co2_gas = gas_co2,
    co2_cemento = cement_co2, co2_quema = flaring_co2,
    poblacion = population, pib = gdp,
    energia_primaria = primary_energy_consumption,
    co2_pc = co2_per_capita, co2_pib = co2_per_gdp,
    energia_pc = energy_per_capita, energia_pib = energy_per_gdp
  )

ingei <- read.csv(file.path(RUTAS$datos, "ingei_energia_serie.csv"), stringsAsFactors = FALSE) |>
  transmute(anio = as.integer(anio),
            ingei_energia_kt = energia_total_ktco2eq,
            ingei_electricidad_kt = electricidad_publica_ktco2eq,
            ingei_transporte_kt = transporte_terrestre_ktco2eq,
            ingei_industria_kt = industrias_manuf_ktco2eq,
            ingei_fugitivas_kt = fugitivas_ktco2eq)

# PIB del Banco Mundial (USD constantes de 2015): serie completa 1960-2024,
# necesaria porque OWID/GCB no reporta PIB para 2023-2024.
wdi <- read.csv(file.path(RUTAS$brutos, "wdi_colombia.csv"), stringsAsFactors = FALSE) |>
  filter(indicador == "pib_usd2015") |>
  transmute(anio = as.integer(anio), pib_const2015 = valor)

upme <- read.csv(file.path(RUTAS$datos, "generacion_anual_por_grupo.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE) |>
  transmute(anio = as.integer(anio),
            gen_total_gwh = total_gwh,
            gen_hidro_gwh = `Hidroeléctrica (convencional + PAH)`,
            gen_biomasa_gwh = Biomasa, gen_solar_gwh = Solar, gen_eolica_gwh = Eólica,
            gen_fosil_gwh = Fósil, gen_otros_gwh = `Otros/no identificado`,
            pct_renovable = pct_renovable_total, pct_fosil, meses_reportados)

ember <- read.csv(file.path(RUTAS$datos, "ember_colombia_indicadores.csv"), stringsAsFactors = FALSE) |>
  mutate(anio = as.integer(anio), fuente = as.character(fuente))
ember_cap <- ember |>
  filter(fuente %in% c("Solar", "Wind", "Hydro", "Bioenergy")) |>
  mutate(fuente = recode(fuente, Solar = "solar", Wind = "eolica",
                         Hydro = "hidro", Bioenergy = "biomasa")) |>
  select(anio, fuente, capacidad_gw) |>
  tidyr::pivot_wider(names_from = fuente, values_from = capacidad_gw,
                     values_fill = 0, names_prefix = "cap_")
ember_sis <- ember |>
  filter(fuente == "Total generation") |>
  transmute(anio, fe_sistema_electrico = intensidad_gco2e_kwh / 1000,
            emisiones_elec_mt = emisiones_mtco2e)

fe_fuente <- ember |>
  filter(fuente %in% c("Coal", "Gas", "Fossil", "Bioenergy"), intensidad_gco2e_kwh > 0) |>
  mutate(fe = intensidad_gco2e_kwh / 1000) |>
  group_by(fuente) |>
  summarise(fe_mediana = median(fe, na.rm = TRUE), .groups = "drop")

cap_fncer <- proy |>
  count(year, wt = capacidad_mw, name = "cap_mw") |>
  arrange(year) |>
  mutate(cap_fncer_acum_mw = cumsum(cap_mw)) |>
  select(anio = year, cap_fncer_acum_mw)

panel <- expand_grid(anio = 1990:2026) |>
  left_join(owid, by = "anio") |>
  left_join(ingei, by = "anio") |>
  left_join(wdi, by = "anio") |>
  left_join(upme, by = "anio") |>
  left_join(ember_cap, by = "anio") |>
  left_join(ember_sis, by = "anio") |>
  left_join(cap_fncer, by = "anio") |>
  arrange(anio)

panel <- panel |>
  mutate(
    across(starts_with("cap_"), ~ tidyr::replace_na(.x, 0)),
    cap_fncer_acum_mw = tidyr::replace_na(cap_fncer_acum_mw, 0),
    gen_fncer_nc_gwh = tidyr::replace_na(gen_solar_gwh, 0) +
      tidyr::replace_na(gen_eolica_gwh, 0) + tidyr::replace_na(gen_biomasa_gwh, 0),
    gen_renovable_gwh = gen_hidro_gwh + gen_fncer_nc_gwh,
    pct_fncer_nc = 100 * gen_fncer_nc_gwh / gen_total_gwh,
    pct_hidro = 100 * gen_hidro_gwh / gen_total_gwh,
    emisiones_evitadas_mt = gen_fncer_nc_gwh * fe_sistema_electrico / 1000,
    fe_ingei_energia = ingei_energia_kt / gen_total_gwh,
    intensidad_carbono = co2_owid / pib,
    co2_lag1 = lag(co2_owid),
    co2_pct_change = 100 * (co2_owid / lag(co2_owid) - 1),
    ghg_per_capita = ghg_total_owid / (poblacion / 1e6)
  )

edgar_nac <- tibble::tibble(
  anio = anios_edgar,
  co2_edgar_mt = c(sum(grilla$co2_2015), sum(grilla$co2_2023), sum(grilla$co2_2025)) / 1e6
)
panel <- panel |> left_join(edgar_nac, by = "anio")

ct_nac <- ct |>
  group_by(grupo) |>
  summarise(mt = sum(emissionsQuantity, na.rm = TRUE) / 1e6, .groups = "drop") |>
  tidyr::pivot_wider(names_from = grupo, values_from = mt, names_prefix = "ct_nac_")
panel <- panel |> mutate(!!!as.list(ct_nac))

saveRDS(panel, file.path(RUTAS$procesados, "panel_nacional.rds"))
msg("panel_nacional.rds: ", nrow(panel), " años x ", ncol(panel), " variables")

# ---- Panel municipal (unidad cuasi-experimental) ---------------------------
paso("4. PANEL MUNICIPAL 2015 -> 2023")
municipios <- grilla |>
  st_drop_geometry() |>
  group_by(mun_id, municipio, depto, region) |>
  summarise(
    n_celdas = n(), area_km2 = sum(area_km2),
    co2_2015 = sum(co2_2015), co2_2023 = sum(co2_2023), co2_2025 = sum(co2_2025),
    poblacion_2015 = sum(poblacion_2015, na.rm = TRUE),
    poblacion_2020 = sum(poblacion_2020, na.rm = TRUE),
    cap_mw = first(cap_mw), n_proy = first(n_proy),
    tratado = first(tratado), dist_fncer_km = min(dist_fncer_km),
    .groups = "drop"
  ) |>
  mutate(
    dens_2015 = co2_2015 / area_km2, dens_2023 = co2_2023 / area_km2,
    d_co2 = co2_2023 - co2_2015,
    log_d_co2 = log1p(pmax(co2_2023, 0)) - log1p(pmax(co2_2015, 0)),
    crec_co2 = ifelse(co2_2015 > 0, (co2_2023 - co2_2015) / co2_2015, NA_real_),
    co2_pc_2015 = ifelse(poblacion_2015 > 0, co2_2015 / poblacion_2015, NA_real_),
    pop_dens = poblacion_2020 / area_km2,
    region = factor(region, levels = names(PAL_REGION)),
    log_pob = log(poblacion_2015 + 1),
    log_area = log(area_km2),
    log_dens_2015 = log1p(dens_2015)
  ) |>
  filter(area_km2 > 0)
saveRDS(municipios, file.path(RUTAS$procesados, "municipios.rds"))
msg("municipios.rds: ", nrow(municipios), " municipios | tratados: ", sum(municipios$tratado))

saveRDS(proy, file.path(RUTAS$procesados, "proyectos_fncer.rds"))
saveRDS(mun_trat, file.path(RUTAS$procesados, "municipios_tratados.rds"))
saveRDS(fe_fuente, file.path(RUTAS$procesados, "factores_emision.rds"))

auditoria <- tibble::tribble(
  ~elemento, ~valor, ~fuente,
  "Celdas del marco muestral (0,1°)", nrow(grilla), "EDGAR v2026 + GADM 4.1",
  "Área cubierta por el marco (km²)", round(sum(grilla$area_km2)), "GADM 4.1",
  "CO2 EDGAR 2015 (Mt)", round(sum(grilla$co2_2015) / 1e6, 2), "EDGAR v2026",
  "CO2 EDGAR 2023 (Mt)", round(sum(grilla$co2_2023) / 1e6, 2), "EDGAR v2026",
  "CO2 EDGAR 2025 (Mt)", round(sum(grilla$co2_2025) / 1e6, 2), "EDGAR v2026",
  "Cobertura del marco vs EDGAR nacional 2023 (%)",
  round(100 * sum(grilla$co2_2023) / extract(r_edgar$y2023, vect(tierra), sum, na.rm = TRUE)[1, 2], 1),
  "EDGAR v2026 (9.657 celdas vs polígono nacional)",
  "CO2 OWID 2023 (Mt)", round(panel$co2_owid[panel$anio == 2023], 2), "OWID / Global Carbon Budget",
  "Climate TRACE 2023 (Mt)", round(sum(ct$emissionsQuantity) / 1e6, 2), "Climate TRACE v5.7",
  "Municipios en el panel", nrow(municipios), "GADM 4.1",
  "Municipios tratados", sum(municipios$tratado), "MinEnergía (meta FNCER)",
  "Proyectos FNCER geolocalizados", sum(!is.na(proy$mun_idx)), "MinEnergía (meta FNCER)",
  "Población 2020 calibrada (M)", round(sum(grilla$poblacion_2020) / 1e6, 2), "WorldPop + OWID"
)
guardar_tabla(auditoria, "00_auditoria_datos")
print(auditoria, n = 20)

paso("PREPROCESAMIENTO COMPLETADO")
