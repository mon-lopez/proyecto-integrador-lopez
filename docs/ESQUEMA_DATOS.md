# Esquema de datos y convenciones del proyecto

Documento de referencia para todos los scripts y agentes. **Leer antes de escribir código.**

## Ejecución

- Raíz: `C:\Users\monic\OneDrive\Escritorio\co2-emisions-col`
- Intérprete (Rscript no está en PATH):
  `& "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" scripts\03_spatial.R`
- Antes de ejecutar: `$env:RENV_CONFIG_SYNCHRONIZED_CHECK="FALSE"` (evita un warning de renv).
- Todo script empieza con `source(file.path("scripts", "00_setup.R"))` y se ejecuta **desde la raíz del proyecto**.
- Paquetes instalados y disponibles: sf, terra, spdep, gstat, spsurvey, spatialreg, sfdep,
  randomForest, ranger, xgboost, caret, ggplot2, dplyr, tidyr, patchwork, viridis, scales, car,
  knitr, rmarkdown. **NO están disponibles**: glmnet, emmeans, agricolae, tidymodels, janitor.
  Para post-hoc usar `stats::TukeyHSD`; para modelos lineales regularizados usar
  `MASS::lm.ridge` o `stats::lm`; para ANOVA tipo III usar `car::Anova`.

## Objetos de `00_setup.R`

- `RAIZ`, `RUTAS` (`$datos`, `$brutos`, `$geometrias`, `$procesados`, `$salida`, `$figuras`, `$mapas`, `$tablas`, `$manuscrito`)
- `CONST`: `crs_geo=4326`, `crs_proj=9377` (MAGNA-SIRGAS / Origen Nacional), `anio_base=2015`,
  `anio_ref=2023`, `anio_horizonte=2030`, `bau_2030_mtco2e=345.8`, `meta_2030_mtco2e=169.44`, `semilla=2026`
- `PAL_REGION` (Caribe, Andina, Pacifica, Orinoquia, Amazonia), `PAL_TEC`, `PAL_ESC`
- `tema()` (tema ggplot), `msg()`, `paso()`, `fmt_num()`, `fmt_p()`, `ic95()`
- `guardar_fig(p, "nombre")` → `output/figures/nombre.png`
- `guardar_tabla(df, "nombre")` → `output/tables/nombre.csv`

## Datos procesados (`datos/procesados/`)

### `grilla_co2.rds` — sf, 9.657 celdas, CRS 9377
Marco muestral = celdas de 0,1° (~11 km, ~123 km²) de EDGAR que intersectan Colombia.
Columnas principales:
`celda_edgar`, `area_km2`, `cobertura` (fracción de Colombia en la celda), `lon`, `lat`,
`x_9377`, `y_9377`, `depto`, `municipio`, `mun_id` (GID_2 GADM), `region` (factor 5 niveles),
`co2_2015`, `co2_2023`, `co2_2025` (t CO2/año), `dens_2015`, `dens_2023`, `dens_2025` (t/km²),
`d_co2`, `crec_co2`, `log_co2_2023`, `log_dens_2023`, `poblacion_2015`, `poblacion_2020`
(calibrada a OWID), `pop_dens_2020`, `urbano` (0/1), `co2_pc_2023`,
`n_proy`, `cap_mw`, `n_solar`, `n_eolico`, `tratado` (0/1), `dist_fncer_km`,
`ct_power`, `ct_transporte`, `ct_edificaciones`, `ct_industria`, `ct_ffo`, `ct_lulucf`,
`ct_residuos`, `ct_otros` (t CO2 Climate TRACE 2023), `grupo_emision` (Baja/Media/Alta),
`estrato` (region_grupo).
Totales: 2015 = 90,0 Mt; 2023 = 102,1 Mt; 2025 = 94,3 Mt. 971 celdas con 0 emisiones.

### `panel_nacional.rds` — 37 filas (1990-2026) x 59 variables
Nacional-año. Variables clave: `co2_owid` (Mt CO2, GCB), `ghg_total_owid` (Mt CO2e),
`co2_carbon`, `co2_petroleo`, `co2_gas`, `co2_cemento`, `poblacion`, `pib`,
`gen_total_gwh`, `gen_hidro_gwh`, `gen_fosil_gwh`, `gen_solar_gwh`, `gen_eolica_gwh`,
`gen_biomasa_gwh`, `gen_fncer_nc_gwh`, `pct_fncer_nc`, `pct_hidro`, `pct_renovable`, `pct_fosil`,
`cap_solar`, `cap_eolica`, `cap_hidro`, `cap_biomasa` (GW), `fe_sistema_electrico` (t CO2e/MWh),
`emisiones_evitadas_mt`, `cap_fncer_acum_mw`, `ingei_energia_kt`, `ingei_electricidad_kt`,
`fe_ingei_energia`, `co2_edgar_mt` (solo 2015/2023/2025), `ct_nac_*` (Climate TRACE 2023).
Cobertura: OWID 1990-2024, UPME 2000-2026 (2026 parcial), Ember 2000-2025.

### `municipios.rds` — 1.008 municipios
Agregación de la grilla por municipio. `tratado` = 1 si el municipio tiene ≥1 proyecto
FNCER (solar/eólico) con entrada ≤ 2023 según MinEnergía (65 tratados).
Variables: `co2_2015`, `co2_2023`, `d_co2`, `log_d_co2` (resultado principal),
`crec_co2`, `dens_2015`, `log_dens_2015`, `co2_pc_2015`, `poblacion_2015`, `poblacion_2020`,
`log_pob`, `pop_dens`, `log_area`, `area_km2`, `cap_mw`, `n_proy`, `dist_fncer_km`,
`region`, `depto`, `n_celdas`.

### Otros
- `proyectos_fncer.rds`: 158 proyectos MinEnergía (tipo, year, capacidad_mw, departamento, municipio, mun_id, emisiones_co2_ton).
- `factores_emision.rds`: FE implícitos por fuente (t CO2e/MWh) de Ember.
- `municipios_tratados.rds`: 69 municipios con proyecto ≤2023.

## Resultados de cada módulo

Cada script debe guardar un `.rds` con la lista de resultados numéricos en
`datos/procesados/` (p. ej. `resultados_espaciales.rds`) **y** un resumen legible en
`output/tables/<modulo>_resumen.csv`. El paper se escribe después leyendo esos objetos:
por eso los nombres de los campos deben ser estables y estar documentados al inicio del script.

## Convenciones de figuras

- `ggplot2` con `tema()`, títulos y ejes **en español**, unidades explícitas.
- Mapas: `geom_sf()` con la grilla en EPSG:9377 (usar `x_9377`/`y_9377` solo para puntos).
- Guardar con `guardar_fig(p, "prefijo_nombre", ancho, alto)`; 300 dpi; fondo blanco.
- Evitar `viridis` para categóricas de región: usar `PAL_REGION`.
- Toda figura debe poder regenerarse ejecutando el script completo (sin estado manual).

## Reglas de calidad

1. Nada de datos simulados: **todo resultado debe provenir de los objetos anteriores**.
2. Reportar tamaños de muestra, unidades y valores p con 3 decimales.
3. Semilla fija (`CONST$semilla`) en todo procedimiento aleatorio.
4. Si un método no converge o da valores degenerados, **decirlo en el log y en el resumen**,
   no inventar un resultado alternativo.
