# Fuentes de datos de emisiones de CO₂ para estadística espacial

> Documento generado el **28 de octubre de 2026**. Todos los enlaces, APIs y
> tamaños de archivo de este documento fueron **probados y verificados** en esa
> fecha (no son referencias de segunda mano).

---

## 1. Recomendación corta

| Lo que necesitas | Fuente recomendada | Por qué |
|---|---|---|
| Panel país–año (regresión, series) | **EDGAR v2026** + **OWID** | 1970–2025, ISO3, licencia abierta |
| Rejilla global para estadística espacial | **EDGAR 0.1°** | malla regular global, 56 años, NetCDF/GeoTIFF |
| Datos de **punto** (plantas, minas, barcos) | **Climate TRACE** | lat/lon por activo; >500 000 activos |
| Alta resolución (1 km) / mensual | **ODIAC** o EDGAR 1 km | fósil, buena para kriging fino |
| Fronteras y polígonos | **Natural Earth / GADM / geoBoundaries** | ya incluidos en `datos/` |

**Combinación que yo usaría para el proyecto:** EDGAR 0.1° como capa primaria
(dato de área), Climate TRACE como capa de puntos, OWID para covariables
país–año y Natural Earth para los límites.

---

## 2. Tabla comparativa

| Fuente | Cobertura | Resolución | Unidad de observación | Años | Licencia |
|---|---|---|---|---|---|
| **EDGAR v2026** (JRC) | global | 0.1° (≈11 km) y 1 km | celda y país | 1970–2025 | CC BY 4.0 ¹ |
| **OWID / Global Carbon Budget** | global | — (país) | país–año | 1750–2024 | CC BY 4.0 |
| **Climate TRACE v5.7** | global | activo individual | punto / polígono | 2015–2025 ² | CC BY 4.0 |
| **Banco Mundial (WDI)** | global | — (país) | país–año | 1960–2024 | CC BY 4.0 |
| **ODIAC 2020b** | global | 1 km, mensual | celda | 2000–2021 | académica (registro) |
| **CAMS / Copernicus** | global | 0.75° / 0.1° | celda (atmósfera) | 2003–hoy | Copernicus |
| **CDIAC** (histórico) | global | país | país–año | 1751–2014 | dominio público |

¹ El CO₂ de combustión proviene de la AIE y está bajo **CC BY-NC-ND 4.0**: no
comercial y sin obras derivadas. Para un trabajo académico no hay problema; para
uso comercial hay que pedir permiso a la AIE (`compliance@iea.org`).
² Climate TRACE tiene años distintos según sector: 2015+ para país y continente,
2021+ para niveles subnacionales.

---

## 3. Fichas detalladas

### 3.1 EDGAR v2026 — la fuente principal

**Emissions Database for Global Atmospheric Research**, Comisión Europea / JRC.
Es la base de referencia mundial para emisiones con distribución espacial.

- Página: <https://edgar.jrc.ec.europa.eu/dataset_ghg2026>
- Repositorio de descarga directa (sin formularios):
  <https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/EDGAR_2026_GHG/>

**Qué contiene**

| Producto | Contenido | Formato |
|---|---|---|
| `IEA_EDGAR_CO2_1970_2025.zip` | CO₂ por país y sector, 1970–2025 | XLSX (hoja `TOTALS BY COUNTRY`) |
| `CO2/TOTALS/emi_nc/` | malla global 0.1°, un archivo por año | NetCDF, ~15 MB/año |
| `CO2/TOTALS/emi_txt/` | la misma malla en texto ESRI ASCII | ZIP, ~36 MB/año |
| `CO2/<SECTOR>/emi_nc/` | malla por sector (ENE, IND, REF_TRF, TNR_Ship…) | NetCDF |
| `EDGAR_2025_1km/` | CO₂ fósil a 1 km para 2024 | ZIP, 1.0 GB |
| `HTAP/`, `edgar_smod/`, `edgar_cities/` | productos derivados | varios |

**Estructura de la malla** (verificada con `ncdump`):

```
dimensions:  lat = 1800 ;  lon = 3600 ;
variables:   float emissions(lat, lon) ;
                 emissions:units       = "Tonnes" ;
                 emissions:global_total = "39.22Gt" ;
                 emissions:year        = "2023" ;
```

- Origen `(-180, 90)`, tamaño de píxel `0.1°`, norte arriba (listo para GDAL/R).
- **Unidad: toneladas de CO₂ por celda por año.**
- 6 480 000 celdas totales; **4.2 millones con emisión > 0**.

**Cómo filtrar**

```r
# por país
edgar_pais <- read.csv("datos/emisiones_co2_pais_edgar.csv")
subset(edgar_pais, iso3 == "COL" & anio == 2023)

# por coordenada: la celda es directa
#   col = floor((lon + 180) / 0.1) + 1 ; fila = floor((90 - lat) / 0.1) + 1
```

**Advertencias importantes** (verificadas empíricamente):

1. **Celdas atípicas en el mar.** Unas 500 celdas concentran ~3 Gt (8 % del
   total mundial) en aguas costeras. La peor: `113.65°E, 21.95°N` (estuario del
   río de las Perlas) con **430 Mt en una sola celda de 0,1°**, dentro del sector
   `ENE` (generación eléctrica). Es un artefacto de asignación espacial de EDGAR,
   no un error de lectura (se comprobó celda por celda contra el NetCDF original
   y aparece igual en todos los años). Consecuencia medida: al agregar la malla
   por país se obtienen **34,89 Gt frente a 39,22 Gt de la tabla oficial**, y la
   diferencia mediana por país es del **10,3 %**. Usa la tabla oficial para
   totales nacionales y la malla para patrones espaciales.
2. **No incluye LULUCF** (uso del suelo y silvicultura). Para eso usa OWID o
   Global Carbon Budget.
3. Los 2–3 últimos años son estimaciones *fast-track*, provisionales.
4. Los búnkeres marítimos y la aviación internacional van como filas `SEA` y
   `AIR` (660 y 553 Mt en 2023), no asignadas a países.

---

### 3.2 Our World in Data / Global Carbon Budget

La mejor fuente para **covariables país–año**: una sola tabla con 79 columnas.

- Descarga: <https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv>
  (14 MB) · también XLSX y JSON.
- Documentación: <https://github.com/owid/co2-data>
- Última actualización del archivo: **2026-06-01**.

**Verificado:** 50 411 filas × 79 columnas, años **1750–2024**, 254 entidades.

Columnas útiles: `co2`, `co2_per_capita`, `co2_per_gdp`, `co2_including_luc`,
`consumption_co2` (emisiones de consumo, no territoriales), `cumulative_co2`,
`coal_co2`, `oil_co2`, `gas_co2`, `cement_co2`, `flaring_co2`, `methane`,
`total_ghg`, `population`, `gdp`.

**Trampa importante:** la columna `country` mezcla países con agregados
regionales y de ingreso (`World`, `Africa`, `High-income countries`…). Filtra
siempre por `iso_code`:

```r
owid <- read.csv("datos/owid-co2-data.csv")
paises_reales <- subset(owid, nzchar(iso_code) & !grepl("^OWID_", iso_code))
```

---

### 3.3 Climate TRACE — la fuente de puntos

Inventario hecho con satélites e IA: cada planta, mina, pozo, barco y área de
deforestación es un **registro con coordenadas**.

- Web: <https://climatetrace.org/>
- **API REST** (sin llave, JSON): `https://api.climatetrace.org/v7/`
- Documentación interactiva: <https://api.climatetrace.org/v7/docs/index.html>
- Volcado completo en Zenodo: <https://doi.org/10.5281/zenodo.20830133>

**Endpoints verificados**

| Endpoint | Devuelve |
|---|---|
| `/v7/sources` | activos con `centroid.longitude` / `centroid.latitude` |
| `/v7/sources/emissions` | totales mensuales agregados por sector |
| `/v7/admins` | áreas administrativas GADM (nivel 0, 1, 2) |
| `/v7/definitions/countries` | catálogo de países |

```bash
# fuentes de Colombia en 2023
curl "https://api.climatetrace.org/v7/sources?gadmId=COL&year=2023&limit=100"

# total nacional
curl "https://api.climatetrace.org/v7/sources/emissions?gadmId=COL&year=2023&gas=co2"
```

Formato de cada registro:

```json
{ "id": 60951997,
  "name": "Colombia_OtherBasins_OtherResources",
  "sector": "fossil-fuel-operations",
  "subsector": "oil-and-gas-production",
  "country": "COL",
  "centroid": { "longitude": -74.290701, "latitude": 5.841388, "srid": 4326 },
  "gas": "co2e_100yr",
  "emissionsQuantity": 10947139.76,
  "year": 2023 }
```

**Volcado completo (Zenodo, 34.6 GB en 12 archivos):**

| Archivo | Tamaño |
|---|---|
| `geometries.gpkg.zip` | 652 MB — **geometrías poligonales**, muy útil para datos de área |
| `power.zip` | 304 MB |
| `fossil_fuel_operations.zip` | 126 MB |
| `manufacturing.zip` | 522 MB |
| `mineral_extraction.zip` | 32 MB |
| `agriculture.zip` | 14.4 GB |
| `buildings.zip` | 5.9 GB |
| `forestry_and_land_use.zip` | 5.4 GB |
| `transportation.zip` | 5.3 GB |
| `waste.zip` | 2.0 GB |

**Advertencia:** la API pagina de 100 en 100 y no acepta filtro por *bounding
box*; para filtrar por coordenada hay que traer los activos y hacer la unión
espacial localmente (o descargar el volcado de Zenodo).

---

### 3.4 Banco Mundial — World Development Indicators

Metadatos limpios y API estable; útil para covariables (PIB, población, energía).

- CO₂ total (excluye LULUCF), indicador `EN.GHG.CO2.MT.CE.AR5`:
  <https://api.worldbank.org/v2/country/COL/indicator/EN.GHG.CO2.MT.CE.AR5?format=json&per_page=1000>
- También `EN.ATM.CO2E.KT`, `EN.ATM.CO2E.PC` (per cápita), `EN.ATM.CO2E.GF.ZS`.
- Años 1960–2024, sin llave, JSON o XML.

---

### 3.5 ODIAC — 1 km mensual

**Open-source Data Inventory for Anthropogenic CO₂**: emisiones fósiles a 1 km
y resolución mensual, 2000–2021 (v2020b).

- Página: <https://www.odiac.org/>
- Distribución: <https://db.cger.nies.go.jp/dataset/ODIAC/>
- Requiere registro/aceptación de términos. Es la referencia para alta
  resolución antes de EDGAR 1 km.
- Formato: GeoTIFF y NetCDF, unidades t C (¡no t CO₂! multiplicar por 3.664).

---

### 3.6 Fuentes atmosféricas (si el enfoque es concentración, no emisión)

- **CAMS** (Copernicus Atmosphere Monitoring Service): reanálisis global de CO₂,
  <https://ads.atmosphere.copernicus.eu/> — requiere cuenta CDS.
- **OCO-2 / OCO-3** (NASA): columnas de XCO₂ satelitales,
  <https://ocov2.jpl.nasa.gov/>.
- **CarbonTracker** (NOAA): asimilación de CO₂ atmosférico.

Sirven si tu proyecto mide *concentraciones* o quiere validar inventarios; no son
inventarios de emisión y no se pueden sumar por país.

---

### 3.7 Fronteras y polígonos (lo que hace posible "filtrar por país")

| Fuente | Niveles | Licencia | Nota |
|---|---|---|---|
| **Natural Earth** | país (1:10m, 50m, 110m) | dominio público | ya incluido en `datos/` |
| **GADM 4.1** | país, provincia, municipio | no comercial | <https://gadm.org/> |
| **geoBoundaries** | ADM0–ADM4 | CC0 | <https://www.geoboundaries.org/> |
| **World Bank Official Boundaries** | país | CC BY 4.0 | |

**Comparación práctica de resoluciones de Natural Earth** (probada con 30 ciudades):

- `ne_10m` → 30/30 aciertos. **Úsala para asignar coordenadas.**
- `ne_50m` → falla en costas densas: Nueva York (Manhattan) queda fuera porque
  el Hudson y el East River se fusionan a esa escala.
- `ne_110m` → solo para mapas de continentes.

**Límite de precisión de Natural Earth 10m.** Es escala 1:10 000 000, así que
suaviza accidentes pequeños. Caso medido: el centro de **Leticia** (Amazonas,
Colombia, `-69.9406, -4.2153`) cae del lado brasileño porque el saliente de
~20 km desaparece en la generalización. Con **GADM nivel 0** de Colombia
(`gadm41_COL_0.json`) la misma coordenada devuelve `COL` correctamente.

Regla práctica:

- Fronteras nacionales "gruesas", mapas y análisis a escala continental →
  **Natural Earth 10m** (rápida, sin registro).
- Puntos individuales cerca de fronteras o costas, o análisis subnacional →
  **GADM** (`gadm41_<ISO3>_1.json`, `_2.json`) o **geoBoundaries**.

Descarga directa de GADM (sin registro):

```bash
# país, provincia y municipio de Colombia
curl -O https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_COL_0.json
curl -O https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_COL_1.json
curl -O https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_COL_2.json
```

`gadm41_COL_2.json` pesa 4,0 MB y trae **1119 municipios**: es la unidad natural
para un análisis espacial subnacional de Colombia.

---

## 4. Los tres patrones de filtrado

### Patrón A — filtrar por país (dato de panel)

```r
library(jsonlite)
edgar <- read.csv("datos/emisiones_co2_pais_edgar.csv")   # iso3, pais, anio, co2_gg
col   <- subset(edgar, iso3 == "COL")
col$co2_mt <- col$co2_gg / 1000
tail(col[, c("anio", "co2_mt")], 5)
```

### Patrón B — filtrar por coordenada (¿en qué país cae este punto?)

Con los polígonos: punto en polígono. Implementación sin dependencias en
`R/00_utilidades.R`, o con `sf`, que es la vía recomendada:

```r
library(sf)
paises <- st_read("datos/paises_ne_10m.geojson")
pts    <- st_as_sf(data.frame(lon = -74.07, lat = 4.71), coords = c("lon","lat"), crs = 4326)
paises[st_intersects(pts, paises)[[1]], "ISO_A3_EH"]   # -> "COL"
```

**Ambas implementaciones coinciden en 20 000 puntos de prueba (100 %).**

### Patrón C — país → subconjunto de la rejilla (dato de área)

```r
library(sf); library(terra)
paises <- st_read("datos/paises_ne_10m.geojson")
col    <- paises[paises$ISO_A3_EH == "COL", ]
r      <- rast("datos/co2_grid_edgar_2023.tif")

# celdas cuyo centroide cae en Colombia
v <- terra::extract(r, vect(st_transform(col, 4326)), cells = TRUE)
sum(v[, 2], na.rm = TRUE) / 1e6    # Mt
```

Para agregación **exacta por área** (repartiendo celdas de borde) usa
[`exactextractr`](https://github.com/isciences/exactextractr):

```r
exactextractr::exact_extract(r, col, "sum")
```

---

## 5. Notas para el diseño del proyecto de estadística espacial

**El tipo de dato define el método:**

| Estructura | Fuente | Métodos |
|---|---|---|
| Retícula regular (lattice) | EDGAR 0.1° | Moran's I, LISA, SAR/SEM, GWR |
| Puntos marcados | Climate TRACE | K de Ripley, kernel, procesos puntuales |
| Polígonos irregulares | países, provincias | contigüidad, pesos por longitud de frontera |
| Geoestadística continua | ODIAC, CAMS | variograma, kriging |

**Cosas que te van a morder si no las prevés:**

1. **Asimetría extrema.** La distribución de emisiones por celda es
   fuertemente asimétrica a la derecha (cuantiles de 2023: mediana 69 t, p99
   67 000 t, p99.99 11.8 Mt). Trabaja en escala **log** o usa métodos robustos;
   el variograma y los índices de Moran son sensibles a esto.
2. **MAUP** (problema de unidad de área modificable): los resultados de
   agregación cambian según el tamaño de celda. Haz análisis de sensibilidad con
   al menos dos resoluciones (0.1° y 0.5°).
3. **Ceros estructurales.** 2.3 millones de celdas tienen exactamente 0
   (océanos, desiertos). No son "datos faltantes": decide si los tratas como
   censura o los excluyes.
4. **Fronteras y costa.** Una celda de borde pertenece parcialmente a dos
   países; `exact_extract` con `weights` resuelve esto mejor que asignar por
   centroide.
5. **Búnkeres.** El transporte internacional (660 + 553 Mt) no pertenece a
   ningún país. Si comparas países sin decirlo, los totales no cuadran.
6. **Efectos escala.** Para GWR o modelos locales necesitas ~30+ unidades; con
   195 países es viable, con 33 provincias de un país también, pero con celdas
   de 0.1° (millones) el costo computacional explota.

---

## 6. Cómo citar

- **EDGAR:** Crippa, M. et al. (2026). *EDGAR v2026 GHG Emissions Database*.
  European Commission, JRC. <https://edgar.jrc.ec.europa.eu/dataset_ghg2026>
- **OWID:** Rosado, P., Ritchie, H., Roser, M. et al. *CO₂ and Greenhouse Gas
  Emissions*. Our World in Data. <https://ourworldindata.org/co2-and-other-greenhouse-gas-emissions>
- **Climate TRACE:** Climate TRACE (2025). *Climate TRACE Emissions Inventory
  V5.7.0*. Zenodo. <https://doi.org/10.5281/zenodo.20830133>
- **Global Carbon Budget:** Friedlingstein, P. et al. (2025). *Global Carbon
  Budget 2025*. ESSD. <https://globalcarbonbudget.org/>
