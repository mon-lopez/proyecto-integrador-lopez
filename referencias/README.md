# Emisiones de CO₂ para un proyecto de estadística espacial

Paquete de datos y código listo para usar. Todo fue **descargado y verificado**
el 28 de octubre de 2026.

- 📄 **[FUENTES_CO2_ESPACIAL.md](FUENTES_CO2_ESPACIAL.md)** — el inventario
  completo de fuentes: cobertura, resolución, licencias, endpoints probados,
  ventajas y trampas de cada una. **Empieza por ahí.**
- 💾 `datos/` — los datos ya descargados y convertidos.
- 🧑‍💻 `R/` — scripts que corren tal cual.

---

## Arranque rápido

```bash
cd emisiones_co2_espacial
Rscript R/01_ejemplo_pais.R          # serie y ranking por país
Rscript R/02_ejemplo_coordenadas.R   # coordenada -> país
Rscript R/03_analisis_espacial.R     # rejilla, agregación y Moran's I
Rscript R/04_climate_trace.R COL 2023   # activos georreferenciados (API)
```

Los scripts 01, 02 y 04 **no requieren instalar nada** (solo R base +
`jsonlite` y `curl`). El 03 necesita `terra` y `sf`, que **ya están instalados**
en la carpeta `.Rlib/` del proyecto (R 4.6, macOS arm64):

```r
# si alguna vez borras .Rlib/ o cambias de máquina:
install.packages(c("terra", "sf"))
```

La carpeta `.Rlib/` es solo una biblioteca local de R (38 MB); los scripts la
detectan solos y no afecta a tu instalación global.



---

## Qué hay en `datos/`

| Archivo | Contenido | Tamaño |
|---|---|---|
| `emisiones_co2_pais_edgar.csv` | **12 544 filas** — `iso3`, `pais`, `anio`, `co2_gg`<br>EDGAR v2026, 1970–2025, 224 países y agregados | 457 KB |
| `co2_grid_edgar_2015.tif` | malla global 0.1°, toneladas de CO₂ por celda | 15 MB |
| `co2_grid_edgar_2023.tif` | ídem 2023 | 15 MB |
| `co2_grid_edgar_2025.tif` | ídem 2025 | 15 MB |
| `owid-co2-data.csv` | OWID: 50 411 filas × 79 indicadores, 1750–2024 | 14 MB |
| `paises_ne_10m.geojson` | Natural Earth 1:10m, 258 entidades — **para asignar coordenadas** | 13 MB |
| `paises_ne_50m.geojson` | Natural Earth 1:50m, 242 entidades — más liviana, menos precisa | 2.9 MB |
| `paises_ne_110m.geojson` | Natural Earth 1:110m, 177 entidades — solo mapas | 819 KB |
| `gadm41_COL_0.json` | GADM Colombia nivel país (ejemplo de mayor precisión) | 211 KB |
| `climate_trace_COL_2023.csv` | **22 835 activos georreferenciados** de Colombia, 2023<br>con `lon`, `lat`, `sector`, `emissionsQuantity` | 5.7 MB |
| `fig_*.png` | figuras generadas por los scripts | — |

**Unidades:** EDGAR entrega **Gg** (= kilotoneladas); divide entre 1000 para
**megatoneladas** (la columna `co2_mt` que calcula `cargar_emisiones_pais()`).
Los GeoTIFF están en **toneladas por celda de 0.1°**.

Los GeoTIFF abren directamente en **QGIS, ArcGIS, R (`terra::rast`) y Python
(`rasterio`)** sin conversión.

---

## Las formas de filtrar (y un bonus de datos de punto)

### 1. Por país → tabla

```r
source("R/00_utilidades.R")
edgar <- cargar_emisiones_pais()
serie_pais("COL", edgar)                       # por ISO3
serie_pais("Colombia", edgar)                  # por nombre
```

### 2. Por coordenada → país

```r
source("R/00_utilidades.R")
paises <- leer_paises(); idx <- construir_indice(paises)

pais_de_coordenada(-74.07, 4.71, paises, idx)  # "COL"
```

Vectorizado y rápido: 200 000 puntos en ~4 s. Verificado contra `sf`
(20 000 puntos, coincidencia del 100 %).

### 3. País → celdas de la rejilla (dato de área)

```r
library(terra); library(sf)
r     <- rast("datos/co2_grid_edgar_2023.tif")
col   <- st_read("datos/paises_ne_10m.geojson") |> subset(ISO_A3_EH == "COL")
terra::extract(r, vect(col), fun = sum, na.rm = TRUE) / 1e6   # Mt
```

Para agregación exacta por área (repartiendo las celdas de borde), usa
`exactextractr::exact_extract(r, col, "sum")`.

### 4. Bonus: datos de punto reales (Climate TRACE)

`datos/climate_trace_COL_2023.csv` ya trae **22 835 activos** de Colombia con
`lon`/`lat` — ideal para análisis de patrones puntuales (K de Ripley, kernel):

```r
act <- read.csv("datos/climate_trace_COL_2023.csv")
act <- act[!is.na(act$lon), ]
sum(act$emissionsQuantity) / 1e6        # 95,05 Mt
sort(table(act$sector), decreasing = TRUE)
```

Para otro país: `Rscript R/04_climate_trace.R BRA 2023` (tarda ~3 min por la
paginación de la API, 100 registros por petición).

---

## Advertencia que conviene leer antes de analizar

Al agregar la malla de EDGAR por país salen **34,89 Gt** frente a **39,22 Gt** de
la tabla oficial (−11 %; diferencia mediana por país 10,3 %). La causa: unas 500
celdas oceánicas concentran ~3 Gt, incluida una celda aislada frente a Hong Kong
(`113.65°E, 21.95°N`) con **430 Mt en 0,1°** dentro del sector de generación
eléctrica. Es un artefacto de asignación espacial de EDGAR, no un error de
lectura: se comprobó celda por celda contra el NetCDF original y se repite en
todos los años.

→ Para **totales nacionales** usa la tabla oficial.
→ Para **patrones espaciales** usa la malla, pero mira las celdas extremas antes
de mapear (la escala logarítmica las esconde y la lineal las aplasta).

---

## Resumen de la recomendación

| Necesidad | Fuente | Estado |
|---|---|---|
| Panel país–año | **EDGAR v2026** + **OWID** | ✅ descargado |
| Rejilla espacial global | **EDGAR 0.1°** | ✅ descargado y convertido |
| Puntos georreferenciados | **Climate TRACE** (`api.climatetrace.org/v7`) | 🔗 listo para consultar |
| Alta resolución 1 km | **ODIAC** / EDGAR 1 km | 🔗 requiere registro |
| Fronteras | **Natural Earth** / **GADM** | ✅ descargado |

Detalle completo, con la justificación de cada elección y las trampas
metodológicas (MAUP, búnkeres, LULUCF, escala log), en
[FUENTES_CO2_ESPACIAL.md](FUENTES_CO2_ESPACIAL.md).
